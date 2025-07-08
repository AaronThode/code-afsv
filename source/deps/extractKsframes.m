function [Kstot,f,t,power]=extractKsframes(x,ovlap,Nfft,chann,frange,Fs,M,vector_only_flag,nowin,threshold,tiltdata,normalize_Ks_flag)

%%%%%%%%%%%%%%%%%%%%%%%extractKsframes.m%%%%%%%%%%%
%function [Kstot,f,t,pwr]=extractKsframes(x,ovlap,Nfft,chann,frange,Fs,M,vector_only_flag,nowin,threshold,tiltdata,normalize_Ks_flag);
%  Aaron Thode
%  July 3, 1996
%  Generates averaged cross-spectral outputs
%    from data input
%  CRUCIAL:  There is a line that flips pgoal so that the first
%		element represents the top phone.
%x=array of data, rows are time, columns are channels
%ovlap=overlap of time samples
%Nfft-number of point desired taken
%chann-vector containing element number to be used:
%frange-pairs of frequencies that define desired ranges If odd, last
%	element is bin spacing (i.e.) evaluate every other bin
% Fs: sampling rate in Hz
% M amount of dataused for the fft snapshot
%
% vector_only:  if 'true' only provide the sum of the squared FFT
%   amplitude, not the CSDM
% nowin-if exists, don't window the data before using fft....used for source signature estimates
% threshold-dB threshold of power (sum of energy across all frequencies and channels) to reject a snapshot.
%    Set to -Inf to ensure all snapshots used.
% tiltdata: tilt: estimated vertical tilt of array in degrees.
%           rd: element depths in m
% normalize_Ks_flag:  If true, divide the CSDM by its trace (keep only
%           phase information)
% Output: Kstot-CSDM (Nel,Nel,freq,time)
%          % t: times in second for each snapshot
%          % f: vector of frequencies
%         However, if 'vector_only_flag' is true
%           Kstot is (Nel, freq, time);
%          % power:  power of signal summed across all channels
%  April 1, 2004: normalize CSDM to have units of power spectral density:
%  |X(f)|^2/(Fs*Nfft), Power/Hz


MAXF=Inf;

if ~exist('normalize_Ks_flag','var')
    normalize_Ks_flag=false;
end

if ~exist('vector_only_flag','var')
    vector_only_flag=false;
end

if ~exist('threshold','var')||isempty(threshold)
    threshold = -Inf;
end


%%Fix tilt
if ~exist('tiltdata','var')||isempty(tiltdata)
    tiltdata.tilt=0; %
    tiltdata.rd=ones(length(chann),1);
end
if size(tiltdata.rd,2)>1
    tiltdata.rd=tiltdata.rd';
end
sintilt=sin(tiltdata.tilt*pi/180);
%%

if ~exist('nowin','var')
    nowin=0;
    disp('The signal will be windowed');
elseif nowin==1
    nowin=1;
    disp('The signal will NOT be windowed');
else
    nowin=0;
    disp('The signal will be windowed');
end

Nel=length(chann);
Ns=floor(((size(x,1)/M)-ovlap)/(1-ovlap));
fprintf('There will be %i snapshots\n',Ns);
%pause;
%Select appropriate frequency bin
%[f,findex]=makefaxis(Nfft,Fs,frange);
frange=round((Nfft/Fs)*frange);
findex=max([1 frange(1)]):min([Nfft/2 frange(2)]);
f=findex*(Fs/Nfft);

if length(findex)<MAXF & ~vector_only_flag
    Kstot=zeros(Nel,Nel,length(findex),Ns);
else
    disp('just making vectors, not matricies');
    Kstot=zeros(Nel,length(findex),Ns);

end
power=zeros(Ns,1);
Nf=length(findex);

Isnap=0:Ns-1; %average all
t=zeros(1,Ns);

if Ns<=0
    disp('Signal too short for one shapnot, will center pad for FFT:');
    %pause;
    Ns=1;Isnap=0;Nx=size(x,1);
    x0=zeros(Nfft,size(x,2));
    index=floor(Nfft/2-(Nx/2));
    index=(index:(index+Nx-1));
    x0(index,:)=x;
    x=x0;
    clear x0
    M=size(x,1);
end
if nowin==0
    win=kaiser(M,2.5);
else
    win=ones(M,1);
end

for I=Isnap
    index=round(I*M*(1-ovlap)+1);
    t(I+1)=index/Fs;
    xindex=(index:(index+M-1));
    xh=x(xindex,chann);
    for Ic=1:length(chann)
        xh(:,Ic)=xh(:,Ic)-mean(xh(:,Ic));

        xh(:,Ic)=xh(:,Ic).*win;
    end
    Xh=fft(xh,Nfft);
    pgoal=Xh(findex,:);
    %Make pgoal a vertical array
    %pgoal=pgoal(:,chann);   %Reject elements ;
    %pgoal=conj(pgoal);   %Test to see if conjugation is the problem
    pgoal=pgoal.';        %Columns are now single-frequency array snapshots
    %pgoal=flipud(pgoal);  %Puts topmost element first, according to lewis

    power(I+1)=sum(sum(abs(pgoal).^2));

    if normalize_Ks_flag
        pgoal=pgoal./sqrt(sum(abs(pgoal).^2));
    end

    if length(findex)<MAXF&&threshold<=10*log10(abs(power(I+1)))
        if vector_only_flag  %%Return only FFT vector
            for J=1:Nf
                %disp(f(J));

                tiltvec=exp(1i*2*pi*f(J)*sintilt*tiltdata.rd/1500);
                ptemp=pgoal(:,J).*tiltvec;
                Kstot(:,J,I+1)=ptemp'; %Top LH cornertop element autocor

            end

        else

            for J=1:Nf
                %disp(f(J));

                tiltvec=exp(1i*2*pi*f(J)*sintilt*tiltdata.rd/1500);
                ptemp=pgoal(:,J).*tiltvec;
                Kstot(:,:,J,I+1)=ptemp*ptemp'; %Top LH cornertop element autocor

            end


        end
    else
        power(I+1)=-1;
    end
end  %I=Isnap

%Normalize to have units of power spectral density...pwr per Hz

power=power./(Fs*Nfft);

if ~normalize_Ks_flag
    if vector_only_flag  %%Return only FFT vector
        Kstot=Kstot/sqrt(Fs*Nfft);
    else
        Kstot=Kstot/(Fs*Nfft);
    end
end


% if length(power)>2
%     figure
%     tt=Isnap*(1-ovlap)*Nfft/Fs;
%     plot(tt,10*log10(abs(power)));
%     grid on
%     if ~isinf(threshold)
%         hold on
%         line([min(tt) max(tt)],threshold*[1 1]);
%     end
%     xlabel('Time (s)');ylabel('dB power');
% end

if length(findex)>MAXF
    Kstot=[];
end




end
