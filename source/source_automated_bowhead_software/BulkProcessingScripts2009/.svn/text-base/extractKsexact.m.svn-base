
%%%%%%%%%%%%%%%%%%%%%%%extractKsexact.m%%%%%%%%%%%
%function [Kstot,Ks_eig,Ns,f,power_est, SNRest,pphase]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin);
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
% nowin-if exists, don't window the data before using fft.used for source signatureestimates
% Output: Kstot-CSDM
%          Ks_eig: CSDM comprised only of primary eigenvalue
%          pwr_est: estimated power
%          SNRest: SNR estimate from ratio of first two eigenvalues of
%          Kstot
%  April 1, 2004: normalize CSDM to have units of power spectral density:
%  |X(f)|^2/(Fs*Nfft), Power/Hz

function [Kstot,Ks_eig,Ns,f,pwr_est,SNRest,pphase]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin)

MAXF=Inf;
if ~exist('nowin'),
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
Ks=zeros(Nel,Nel);
%pause;
%Select appropriate frequency bin
%[f,findex]=makefaxis(Nfft,Fs,frange);
frange=1+floor((Nfft/Fs)*frange);
findex=max([1 frange(1)]):min([Nfft/2 frange(2)]);
f=(findex-1)*(Fs/Nfft);

if length(findex)<MAXF,
    Kstot=zeros(Nel,Nel,length(findex));
else
    disp('frange too long to make Kstot\n just making pgoal');
end
power=zeros(Nel,length(findex));
Nf=length(findex);
if Isnap<0|~exist('Isnap'),
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

for I=Isnap,
    index=round(I*M*(1-ovlap)+1);
    xindex=(index:(index+M-1));
    xh=x(xindex,chann);
    for Ic=1:length(chann),
        xh(:,Ic)=xh(:,Ic)-mean(xh(:,Ic));
        xh(:,Ic)=xh(:,Ic).*win;
    end
    Xh=fft(xh,Nfft);
    pgoal=Xh(findex,:);
    %Make pgoal a vertical array
    %pgoal=pgoal(:,chann);   %Reject elements ;
    %pgoal=conj(pgoal);   %Test to see if conjugation is the problem
    pgoal=pgoal.';        %Rotate so vector is vertical, like KRAKEN
    %pgoal=flipud(pgoal);  %Puts topmost element first, according to lewis
    
    
    power=power+(abs(pgoal).^2);
    
    if length(findex)<MAXF,
        for J=1:length(findex),
            %disp(f(J));
            Kstemp=pgoal(:,J)*pgoal(:,J)'; %Top LH cornertop element autocor
            Kstot(:,:,J)=Kstot(:,:,J)+Kstemp;
        end
    end
end

if length(findex)>MAXF,
    Kstot=[];
end

%Normalize to have units of power spectral density...pwr per Hz
Kstot=Kstot/(Fs*Nfft);

SNRest=zeros(size(f));
pwr_est=zeros(size(f));
Ks_eig=zeros(length(chann),length(chann),length(f));

%Compute eigenvalues of CSDM, and the ratio of the high-power to the
%low-power sound be like SNR
for If=1:length(f),
    pvec(1:Nel,If)=diag(Kstot(:,:,If));
    pphase(1:Nel,If)=angle(squeeze(Kstot(:,1,If)))-angle(pvec(1,If));
    [VV,EE]=eig(Kstot(:,:,If));
    EE=diag(EE);
    [EE_sort,Isort]=sort(EE,1,'descend');
    SNRest(If)=10*log10(abs(EE_sort(1)./sum(EE_sort(2:end))));
    pwr_est(If)=10*log10(abs(EE_sort(1)));
    Eigvec(:,If)=VV(:,Isort(1));
    Ks_eig(:,:,If)=VV(:,Isort(1))*VV(:,Isort(1))';
    %keyboard;
end

if Ns<=1,
    SNRest=-SNRest;
end
%disp('done');pause
%for J=1:length(findex),
%	Ks=squeeze(Kstot(:,:,J));
%	[V,D,FLAG]=eigs(Ks,1);
%	disp(D)
%	Kstot(:,:,J)=Ks/D;
%	Kstot(:,:,J)=Ks/trace(Ks);
%end

%peaks=peakpick(power,.25*max(power));
%disp(sprintf('Peak frequency at %6.2f\n',f(peaks)))
