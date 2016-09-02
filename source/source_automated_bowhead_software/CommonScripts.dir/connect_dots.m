function detect_out=connect_dots(connect,detect,T,F,Nfft,Fs,Idebug),
%%function detect_out=connect_dots(connect,detect,T,F,Nfft,Fs,Idebug),
%% Connect peak detections into segments
%% First version Oct. 31 stopped after finding a single contour
%%  Later Oct. 31 version finds all contours...
%%INPUTS:
%% connect.gap_t: parameter for maximum time jump allowed
%% connect.gap_f: parameter for maximum frequency jump allowed.
%% connect.total_duration.  minimum length of entire contour
%%   detect.matrix: previous detection matrix.
%% T,F, times and frequencies of corresponding specgram
%% Nfft, Fs, fft size and sampling frequency (Hz)

if nargin<7,
    Idebug=0;
end
detect_out=detect;
Igt=floor(connect.gap_t./(T(2)-T(1))); %Time gap in samples
Igf=round(connect.gap_f*Nfft/Fs); %Freq gap in samples
Ntduration=floor(connect.total_duration./(T(2)-T(1)));
Ntmaxt=floor(connect.max_duration./(T(2)-T(1)));
Icall=0;
status=zeros(size(detect.matrix));
status_test=2;
status_good=1;

for I=1:(length(T)-1),
    current_slice=detect.matrix(:,I);
    Igood=find(current_slice>1); %peak locations in this slice
    
    for If=1:length(Igood),
        branch_status=status_test;  %NEW
        while branch_status>status_good, %NEW
            Icurrent=current_slice(Igood(If)); %Frequency index to begin with
            
            Iseg=0;
            Itmark=I+1; %Current time
            Igap=0;
            Igtest=Igf; % Initialize freqeuncy box
            TT=[T(I)];
            FF=[F(Icurrent)];
            status(If,I)=status_test;
            %disp(sprintf('Current frequency bin to match is %i',Icurrent));
            while Igap<Igt&Itmark<length(T),  %While gap from end of previous segnment within bounds and not at end of time window.
                Ifuture_slice=detect.matrix(:,Itmark);
                % AARON change Nov 1, 2006
                Ipossible=find(abs(Icurrent-Ifuture_slice)<=Igtest);  %NEW
                if length(Ipossible)>1, %NEW
                       branch_status=branch_status+1;%NEW    
                end  %NEW
                errr=abs(Icurrent-Ifuture_slice(Ipossible));  %NEW
                
                %%errr=abs(Icurrent-Ifuture_slice); COmmented out recently
                [small_err,Ibest]=min(errr);%Frequency closest to current one, should also match to power..
                Ibest=Ipossible(Ibest);
                if small_err<=Igtest,
                    Igap=1;Igtest=Igf;
                    %status(Ibest,Itmark)=status_test;
                    status(Ibest,Itmark)=branch_status;  %NEW
                    if Ifuture_slice(Ibest)==0,
                        disp('Unexplained error');
                        Igap=Inf;
                    end
                    Icurrent=Ifuture_slice(Ibest); %New frequency bin
                  try, 
                    TT=[TT T(Itmark)];
                    FF=[FF F(Icurrent)];
	          catch,
		   disp('connect dots crash, Icurrent less than zero');
		    detect_out=[];
		    return;
		  end
                else
                    Igap=Igap+1; %We note that we have not found what we are looking for.
                    %Igtest=Igtest+Igf; %Widen tolerence in case signal has slope.
                    %keyboard;
                end
                Itmark=Itmark+1;  %New time window
            end %while Igap
             
            %Evaluate segment length
            %Istatus=find(status==status_test);
            
            Istatus=find(status>status_good); %NEW
            Nseg=length(Istatus);
            if Idebug>0,disp(sprintf('Current segment length: %i required length %i to %i',Nseg,Ntduration,Ntmaxt));end
            
      
            if Nseg>=Ntduration&Nseg<=Ntmaxt, %It has lasted long enough
                %keyboard;
                Icall=Icall+1;
                status(Istatus)=status_good;
                detect.matrix(Istatus)=-Inf; %Erase freqeuncies..
                detect_out.score(Icall)=Nseg;
                detect_out.T{Icall}=TT;
                detect_out.F{Icall}=FF;
                branch_status=1;  %NEW-Exit while loop
                %keyboard;
            else
                if branch_status>2,
                        %keyboard;
                    end
                Ibadd=find(status==branch_status);  %NEW, Eliminate last branch
                detect.matrix(Ibadd)=-100; %No point retracing this territory.  Mistake?
                %Istatus=find(status>status_good);
                status(Istatus)=0; %Return to pool   
                
                if branch_status==status_test, %%NEW, If no branches have occured
                    branch_status=1;  %NEW, exit while loop
                else %NEW
                    branch_status=status_test;%NEW
                end%NEW
            end  %if Nseg
        end % NEW-end while loop for branches...
    end %If, length(Igood)
end %I, length(T)
Ifix=find(detect_out.matrix<1);
detect_out.matrix(Ifix)=1;
%keyboard;
if Icall==0,
    detect_out.score=[];
    detect_out.T=[];
    detect_out.F=[];
end
