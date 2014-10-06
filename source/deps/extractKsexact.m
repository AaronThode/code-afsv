function [Kstot,f,VV,EE_sort]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin,threshold,tiltdata)

%%%%%%%%%%%%%%%%%%%%%%%extractKsexact.m%%%%%%%%%%%
%function [Kstot,Ks_eig,Ns,f,power]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin);
%[Kstot,pgoal,Ns,f,pwr_est,SNRest,pphase]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin);
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
%Isnap %Select the Isnap window.  If negative, average
%M amount of dataused for the fft snapshot
% nowin-if exists, don't window the data before using fft.used for source signature estimates
% threshold-dB threshold of power (sum of energy across all frequencies and channels) to reject a snapshot.
%    Set to Inf to ensure all snapshots used.
% tiltdata: tilt: estimated vertical tilt of array in degrees.
%           rd: element depths in m
% Output: Kstot-CSDM
%          %V: eigenvectors of CSDM
%          %ev: eigenvalues of CSDM
%  April 1, 2004: normalize CSDM to have units of power spectral density:
%  |X(f)|^2/(Fs*Nfft), Power/Hz
VV=[];
EE_sort=[];
MAXF=Inf;
if ~exist('threshold','var')
    threshold = -Inf;
end

if ~exist('tilt','var')
    tiltdata.tilt=0; %
    tiltdata.rd=ones(length(chann),1);
end
if size(tiltdata.rd,2)>1
    tiltdata.rd=tiltdata.rd';
end
sintilt=sin(tiltdata.tilt*pi/180);
if ~exist('nowin','var')
    nowin=0;
    disp('The signal will be windowed');
elseif nowin==1,
    nowin=1;
    disp('The signal will NOT be windowed');
else
    nowin=0;
    disp('The signal will be windowed');
end

Nel=length(chann);
Ns=floor(((size(x,1)/M)-ovlap)/(1-ovlap));
disp(Ns);
%pause;
%Select appropriate frequency bin
%[f,findex]=makefaxis(Nfft,Fs,frange);
frange=round((Nfft/Fs)*frange);
findex=max([1 frange(1)]):min([Nfft/2 frange(2)]);
f=findex*(Fs/Nfft);

if length(findex)<MAXF
    Kstot=zeros(Nel,Nel,length(findex));
else
    disp('frange too long to make Kstot\n just making pgoal');
end
power=zeros(Ns,1);
Nf=length(findex);
if	~exist('Isnap', 'var')|| Isnap<0
    Isnap=0:Ns-1; %average all
end

if Ns<=0,
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
if nowin==0,
    win=kaiser(M,2.5);
else
    win=ones(M,1);
end

for I=Isnap
    index=round(I*M*(1-ovlap)+1);
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
    
    
    power(I+1)=(Fs/Nfft)*sum(sum(abs(pgoal).^2))/(Fs*Nfft);
    
    if length(findex)<MAXF&&threshold<=10*log10(abs(power(I+1)))
        for J=1:Nf
            %disp(f(J));
            
            tiltvec=exp(1i*2*pi*f(J)*sintilt*tiltdata.rd/1500);
            ptemp=pgoal(:,J).*tiltvec;
            Kstemp=ptemp*ptemp'; %Top LH cornertop element autocor
            Kstot(:,:,J)=Kstot(:,:,J)+Kstemp;
            
        end
    else
        power(I+1)=NaN;
    end
end  %I=Isnap

figure
tt=Isnap*(1-ovlap)*Nfft/Fs;
plot(tt,10*log10(abs(power)));
grid on
if ~isinf(threshold)
    hold on
    line([min(tt) max(tt)],threshold*[1 1]);
end
xlabel('Time (s)');ylabel('dB power');

if length(findex)>MAXF,
    Kstot=[];
end

%Normalize to have units of power spectral density...pwr per Hz
Kstot=Kstot/(Fs*Nfft);

%Compute eigenvalues of CSDM, and the ratio of the high-power to the
%low-power sound be like SNR


%SNR=zeros(1,Nf);
if Ns>1
    EE_sort=zeros(Nel,Nf);
    VV=zeros(Nel,Nel,Nf);
    for If=1:Nf
        [VV0,EE]=eig(Kstot(:,:,If));
        EE=diag(EE);
        [EE_sort(:,If),Isort]=sort(EE,1,'descend');
        VV(:,:,If)=VV0(:,Isort);
        
    end
end

% if Ns<=1
%     SNRest=-SNRest;
% end

end
