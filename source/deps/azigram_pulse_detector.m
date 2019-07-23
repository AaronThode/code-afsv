function output=azigram_pulse_detector(x,Fs,Nfft,ovlap,params)
%%%Inputs:x:  cell matrix with x{Idasar} being a N by 3 column matrix from
%%%         DASAR Idasar
%%%%Parameters
%params.tabs_start:  %%Start time of x in datenumber.
%params.threshold=5*4;  %Hz
%params.f_transition;  %Frequency at which DASAR p and v go out of phase.
%params.vertical_line=round(params.vertical_line/dF);
%params.min_TBW=round(params.min_TBW/(dT*dF));
%params.brefa=vector the same size as x;
%params.frange
%params.debug
%params.a_grid, da
%params.threshold
%params.time_chunk=10*60;  %How much data to read into memory at a
%           time--solves clock drift problem too

params_temp=params;
Nstation=length(x);  %How many DASARs are participating?

%

%%%Determine how big a chunk to import

params.time_chunk=min([params.time_chunk length(x{1})/Fs]);
%Nt=floor(params.time_chunk/dT);
%Npt=floor(1+((length(x{1})/Nfft)-1)/(1-ovlap));
Iwindow=unique([1:(params.time_chunk*Fs):length(x{1}) length(x{1})]);  %Index of chunks

if params.time_chunk==length(x{1})/Fs  %%%Chunk is whole window
    Iwindow=[1 length(x{1})+1];
end

%%%Index of detections
Icount=0;

dT=Nfft*(1-ovlap);
dF=Fs/Nfft;
params.vertical_line=round(params.vertical_line/dF);
params.min_TBW=round(params.min_TBW/(dT*dF));

SE=strel('line',params.vertical_line,90);  %%%closing object
    
for It=1:(length(Iwindow)-1)
    fprintf('Processing chunk %i of %i\n',It,length(Iwindow)-1);
    Iindex=Iwindow(It):(Iwindow(It+1)-1);
    
    
    for Ichan=1:Nstation
        params_temp.brefa=params.brefa(Ichan);
        [TT,FF,azi{Ichan},B{Ichan},~,Ix{Ichan},Iy{Ichan}]=compute_directional_metrics ...
            (x{Ichan}(:,Iindex),{'Directionality','Directionality','ItoERatio','ItoERatio'},Fs,Nfft,ovlap,params_temp, ...
            'gsi',[false true false true]);
        [~,Icut]=min(abs(FF-params.f_transition));
        
        %     [~,~,temp,~,~,Ixr,Iyr]=compute_directional_metrics ...
        %         (x{Ichan},'Directionality',Fs,Nfft,ovlap,params_temp, ...
        %         'gsi',true);
        
        %%%Correct for reactive/active intensity flip around 300 Hz...
        azi{Ichan}{1}((Icut+1):end,:)=azi{Ichan}{2}((Icut+1):end,:);
        Ix{Ichan}((Icut+1):end,:)=1i*Ix{Ichan}((Icut+1):end,:);  %Shift 90 degrees
        Iy{Ichan}((Icut+1):end,:)=1i*Iy{Ichan}((Icut+1):end,:);  %Shift 90 degrees
        ItoE{Ichan}=azi{Ichan}{3};  %Transport ratio
        ItoE{Ichan}((Icut+1):end,:)=azi{Ichan}{4}((Icut+1):end,:);  %%Correct for reactive error
        azi{Ichan}=azi{Ichan}{1};
    end
    
    
    %%%%Compute azigram difference
    %%%Convert parameter values into pixels
    dF=FF(2)-FF(1);
    dT=TT(2)-TT(1);
    
    
    
    %%%%If the image is very long, break up into smaller pieces for this loop
    %%%%Crude time synchronization....
    % Try to time-align spectrograms to best of ability..
    maxlag=1024;
    
    for Ichan=2:Nstation
        
        temp=zeros(1,2*maxlag+1);
        for Ir=1:size(B{Ichan},1) %For each row..
            [tmp,laggs]=xcov((B{1}(Ir,:)),(B{Ichan}(Ir,:)),maxlag);
            temp=temp+tmp;
        end
        clear tmp
        [~,Ioff]=max(temp);
        fprintf('Best-fit lag for channel %i is %i samples or %6.2f sec\n',Ichan,laggs(Ioff),laggs(Ioff)*dT);
        
        if laggs(Ioff)<0
            
            IIII=-laggs(Ioff)+1;
            delta_T(It,Ichan)=IIII*dT;
            azi{Ichan}=azi{Ichan}(:,IIII:end);
            ItoE{Ichan}=ItoE{Ichan}(:,IIII:end);
            B{Ichan}=B{Ichan}(:,IIII:end);
            
            for JJJ=1:Nstation
                if JJJ~=Ichan
                    IIII=(size(azi{JJJ},2)+laggs(Ioff));
                    azi{JJJ}=azi{JJJ}(:,1:IIII);  %Whatever you do to azi{1} you do to channels<Ichan
                    ItoE{JJJ}=ItoE{JJJ}(:,1:IIII);
                    B{JJJ}=B{JJJ}(:,1:IIII);
                end
            end
             
            TT=TT(1:size(azi{1},2));
            
        else
            % keyboard
            delta_T(It,Ichan)=laggs(Ioff)*dT;
            
            IIII=(size(azi{Ichan},2)-laggs(Ioff));
                    
            azi{Ichan}=azi{Ichan}(:,1:IIII);
            B{Ichan}=B{Ichan}(:,1:IIII);
            ItoE{Ichan}=ItoE{Ichan}(:,1:IIII);
            
            IIII=laggs(Ioff)+1;
            
            for JJJ=1:Nstation
                if JJJ~=Ichan
                    azi{JJJ}=azi{JJJ}(:,IIII:end);
                    B{JJJ}=B{JJJ}(:,IIII:end);
                    ItoE{JJJ}=ItoE{JJJ}(:,IIII:end);
                end
            end
             
            TT=TT(IIII:end);
            
        end
        fprintf('time shift is %6.2f\n',delta_T(It,Ichan));
        
        %%%Debug check
        temp=zeros(1,2*maxlag+1);
        for Ir=1:size(B{Ichan},1) %For each row..
            [tmp,laggs]=xcov((B{1}(Ir,:)),(B{Ichan}(Ir,:)),maxlag);
            temp=temp+tmp;
        end
        clear tmp
        [~,Ioff]=max(temp);
        fprintf('Best-fit lag for channel %i is %i samples or %6.2f sec\n',Ichan,laggs(Ioff),laggs(Ioff)*dT);
        
        
    end  %Ichan
    if params.debug.morph||params.debug.mask_construction
        figure(5);
        plot(temp,'k');grid on;
        plot_azigram_overview;
        keyboard;
        close(5);
    end
    %%%Create a series of binary images by azimuthal sector
    
    Iff=find(FF>=params.frange(1)&FF<=params.frange(2));
    for I=1:Nstation
        bw_all{I}=false(length(Iff),size(azi{1},2),length(params.a_grid));
        for Ia=1:length(params.a_grid)
            %Make binary images
            %fprintf('Searching between %6.2f and %6.2f degrees\n',params.a_grid(Ia)-params.da/2,params.a_grid(Ia)+params.da/2);
            bw_all{I}(:,:,Ia)=(azi{I}(Iff,:)>=(params.a_grid(Ia)-params.da/2))&(azi{I}(Iff,:)<=(params.a_grid(Ia)+params.da/2));
        end
    end
    
    
    %Ia=9;  %Demo for 10 seconds into DASAR A, 1/18 Tako City
    %Ia=13;  %Demo for 10 seconds into DASAR A, 1/18 Tako City
    err=zeros(length(params.a_grid),length(params.a_grid),length(TT),3,'uint8');
    disp('Starting a_grid');
    
    for Ia=1:length(params.a_grid)
        fprintf('Image processing azi %i degrees: %6.2f percent done\n',params.a_grid(Ia),100*Ia/length(params.a_grid));
        
        %%%Remove small objects and store original result
        bw{1}.raw=squeeze(bw_all{1}(:,:,Ia));
        bw{1}.clean=bwareaopen(bw{1}.raw,params.min_TBW);
        
        %%%keep clean image
        bw_all{1}(:,:,Ia)=bw{1}.clean;
        
        for Ia1=1:length(params.a_grid)
            %    for Ia1=2:6
            bw{2}.raw=squeeze(bw_all{2}(:,:,Ia1));
            bw{2}.clean=bwareaopen(bw{2}.raw,params.min_TBW);
            bw_all{2}(:,:,Ia1)=bw{2}.clean;
            
            comb.raw.noclose=bw{1}.raw&bw{2}.raw;
            comb.clean.noclose=bw_all{1}(:,:,Ia)&bw_all{2}(:,:,Ia1);
            
            %%%Try to connect gaps..
            %comb.raw.close=imclose(comb.raw.noclose,SE);
            %comb.clean.close=imclose(comb.clean.noclose,SE);
            %comb.clean.close=imerode(comb.clean.close,SE);
            %comb_close=bwareaopen(comb_clean,vertical_line);
            
            %%%Diagnostic figures
            err(Ia,Ia1,:,1)=sum(comb.raw.noclose);
            err(Ia,Ia1,:,2)=sum(comb.clean.noclose);
            %err(Ia,Ia1,:,3)=sum(comb.clean.close);
            if params.debug.morph
                if any(sum(comb.clean.noclose)>params.threshold/dF)
                    %Needs TT FF Iff azi azi{2} bw bw1 comb dBspread  B  B{2}  params
                    %Ia Ia 1
                    
                    plot_morphological_debug;
                    pause
                end
            end
            
        end  %Ia1
    end %Ia
    
    err=err*dF;  %Converted back to Hz (bandwidth)
    
    %Ngrid=length(params.a_grid);
    
    detect=(squeeze(err(:,:,:,2))); %cleaned surface, but not closed (comb.clean.noclose)
    %detection units (azimuth1, azimuth2, time)
    detect=[detect(end,:,:); detect; detect(1,:,:)]; %Takes care of cyclical permutation around north
    detect=[detect(:,end,:) detect detect(:,1,:)];
    
    %%Augment bw_all images by cyclical permuation
    for I=1:Nstation
        bw_all{I}(:,:,2:(end+1))=bw_all{I};
        bw_all{I}(:,:,1)=bw_all{I}(:,:,end);
        bw_all{I}(:,:,end+1)=bw_all{I}(:,:,2); %I had made a mistake here, June 27
    end
    
    % I1=2;I2=3;
    %  temp=(bw_all{1}(:,:,I1)&bw_all{2}(:,:,I2));
    % figure;plot(dF*sum(temp));hold on;plot(squeeze(detect(I1,I2,:)),'r');
    
    
    %%%%%%%%For every combination of angles, search for common
    %%%%%%%%  pulses.
    
    %for Ia=2:(length(params.a_grid)-1)
    for Ia=2:(length(params.a_grid)+1)  %The 2 is because we'd artifically bounded the detect function with
        % a wrap around
        
        for Ia1=2:(length(params.a_grid)+1)
            det=squeeze(detect(Ia,Ia1,:));  %Now just a time series with amplitude in Hz
            %fprintf('Angle 1: %6.2f Angle 2: %6.2f\n',params.a_grid(Ia-1),params.a_grid(Ia1-1));
            Itest=find(det>params.threshold);
            if isempty(Itest)  %%%If no pulses found, move on
                continue
            end
            
            %%%Identify distinct pulses by looking for time gaps
            Ibound=[0; find(diff(Itest)>2); length(Itest)];
            fprintf('There are %i pulses in this segment\n',length(Ibound));
            for Ipulse=1:(length(Ibound)-1)  %%%For each pulse...
                %figure(20);plot(det);grid on;title(sprintf('%i %i: %i segments',Ia,Ia1,length(Ibound)-1));
                if rem(Ipulse,100)==0,fprintf('%6.2f percent done\n',100*Ipulse/length(Ibound));end
                II=(Ibound(Ipulse)+1):(Ibound(Ipulse+1));
                index=Itest(II);  %%%indicies in det that define pulse..
                
                %%%Check that a higher score does not exist in an adjacent bin.
                %%%% I looked at both area of peak and height.  Problem with
                %%%% area is that a peak might split into two in an adjacent
                %%%% bin (index differs between bins).
                %%%  The peak height is a more robust check of "best"
                %%%  viewpoint.
                %score=sum(det(index));
                score=max(det(index));
                local_peak_flag=true;
                
                %%%Check if we are near the "best" combination of Ia and Ia1
                % for a particular peak
                for JJ=[-1 0 1]
                    for KK=[-1 0 1]
                        %temp=sum(detect(Ia+JJ,Ia1+KK,index),3);
                        temp=max(detect(Ia+JJ,Ia1+KK,index));
                        %fprintf('Ia: %i Ia1: %i, JJ:%i KK:%i, score: %6.2f\n',Ia, Ia1, JJ,KK,temp);
                        local_peak_flag=local_peak_flag&(score>=temp);
                    end
                end
                if ~local_peak_flag  %%Haven't hit the local maximum in error surface
                    % close(20)
                    continue
                end
                
                %Create a mask from adjacent azimuths.  Be careful to wrap around 0
                %azimuth
                index_all=index;index2=1:length(index_all);
                %index_all=1:length(TT);
                mask=false(size(bw_all{1}(:,index_all,1)));
                for JJ=[-1 0 1]
                    I1=Ia+JJ;
                    for KK=[-1 0 1]
                        I2=Ia1+KK;
                        temp=(bw_all{1}(:,index_all,I1)&bw_all{2}(:,index_all,I2));
                        
                        if params.debug.mask_construction
                            plot_mask_construction;
                        end
                        mask=(mask |temp); %%Grow mask with contributions from all sides
                        
                        %%%Clear detection function so won't trigger in future
                        %%% June 27, 2019
                        detect(I1,I2,index_all)=0;
                        
                    end
                end
                
                Icount=Icount+1;
                
                for I=1:Nstation
                    %%%%Weighted_median
                    
                    [output.azi.wm.med(I,Icount),output.azi.wm.iqr(I,Icount)]=get_weighted_median(mask(:,index2), ...
                        azi{I}(Iff,index),B{I}(Iff,index),params.da);
                    
                    %%%Compute total energy in two channels
                    output.azi.avg(I,Icount) = atan2d(sum(sum(mask(:,index2).*real(Ix{I}(Iff,index)))), ...
                        sum(sum(mask(:,index2).*real(Iy{I}(Iff,index)))));
                    output.azi.avg(I,Icount)=bnorm(params.brefa(I)+output.azi.avg(I,Icount));
                    
                end
                
                Ifr=(Iff(any(mask(:,index2),2)));
                output.freq.min(1,Icount)=FF(Ifr(1));
                output.freq.max(1,Icount)=FF(Ifr(end));
                
                temp=mask(:,index2).*(10.^(0.1*B{1}(Iff,index)));
                temp=temp(temp~=0);
                output.power.SEL(1,Icount)=10*log10(sum(temp));
                output.power.PSD.med(1,Icount)=10*log10(median(temp));
                output.power.PSD.iqr(1,Icount)=iqr(10*log10(temp));
                
                temp=mask(:,index2).*ItoE{1}(Iff,index);
                temp=temp(temp~=0);
                output.ItoE.med(1,Icount)=median(temp);
                output.ItoE.iqr(1,Icount)=iqr(temp);
                
                for Istation=1:Nstation
                    %output.tabs(1,Icount)=datenum(0,0,0,0,0,TT(index(1))+0.5*dT+Iindex(1)/Fs)+params.tabs_start;  %%%Need start of time chunk
                    output.tabs(Istation,Icount)=datenum(0,0,0,0,0,TT(index(1))+0.5*dT-delta_T(It,Istation)+Iindex(1)/Fs)+params.tabs_start;
                end
                
                output.duration.max(1,Icount)=dT*length(index);
                output.BW.max(1,Icount)=double(score);
                output.TWP.max(1,Icount)=dT*sum(double(det(index)));
                
                %output.BW.max_close(1,Icount)=dF*double(max(squeeze(err(Ia,Ia1,index,3))));
                %output.TWP.max_close(1,Icount)=dF*dT*sum(double(squeeze(err(Ia,Ia1,index,3))));
                
                %%%To speed up, zero out detections that are adjacent here...
                %detect(Ia+(-1:1),Ia1+(-1:1),index)=0;
            end %Ipulse
            
            
        end  %Ia1
    end %Ia
    
    
end  %It

output.delta_T=delta_T;




