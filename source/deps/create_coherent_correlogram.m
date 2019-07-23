function [XC_eq,Trel,tt,pwr,pwr_tot,yscale]= create_coherent_correlogram(x,fs,param,flo,fhi)
%function [mean_corr,Trel,tt,pwr]= create_coherent_correlogram(x,fs,param,flo,fhi)
% Trel: relative time
% tt: lag times...

ovlap=round(param.ovlap*param.Nfft);
Nfft=param.Nfft;
teager=param.teager;  %Apply teager-kaiser operation to spectrogram before processing...
ici_range=param.ici_range;
noiseT=param.noiseT;
alpha=param.alpha;

frange=[0.8*flo flo fhi fhi+0.2*flo];
[N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.05],fs);
Bfilt = firpm(N,Fo,Ao,W);
y=filtfilt(Bfilt,1,x-mean(x));
%Adaptive filter option

if teager==1
    y=y(2:end-1).^2-(y(1:end-2).*y(3:end));

end

index1=1:round((Nfft-ovlap)):length(y);
Ncoll=length(index1)-1;
XC=zeros(Nfft/2,Ncoll);
pwrr=zeros(1,Ncoll);
for I=1:Ncoll
    try
        Indexx=index1(I)+(0:(Nfft/2-1));
        tmp=(xcorr(y(Indexx),'coef'));  %Delphine, note that this should be windowed...
        tmp=fftshift(tmp);
        XC(:,I)=tmp(1:(Nfft)/2);
        pwrr(I)=tmp(end);
        %%Mulitpath can generate strong correlations.  Surface-reflected paths should have negative
        %%correlations.  Therefore, only keep levels greater than zero.
        %Iplus=find(XC(:,I)<=0);
        %XC(Iplus,I)=0;
    catch
        disp('Problem');
    end
end
tt=(1:(Nfft/2))/fs;
dT=ovlap/fs;
%Tabs=twant(Itimes)+datenum(0,0,0,0,0,index1/fs);
Trel=index1/fs;

%%Trim away portions of autocorrelation that will not have creaks (ici_range)
Iwant=find(tt>ici_range(1)&tt<ici_range(2));

%Remove times not of interest, normalize zero lag to one.
XC=XC(Iwant,:);
%XC0=XC0(Iwant,:);
%signn=signn(Iwant,:);
tt=tt(Iwant);
%Quick click removal...

Inoise=round(fs*noiseT);  %y-axis bins to estimate autocorrelation level


%%Adaptively equalize along y-axis

XC_eq=zeros(size(XC));
eq=mean(XC(2:Inoise,:));  %Original equalization estimate

for I=1:size(XC,1)
    eq=(1-alpha).*eq+alpha.*XC(I,:);
    XC_eq(I,:)=XC(I,:)-eq;
end



end