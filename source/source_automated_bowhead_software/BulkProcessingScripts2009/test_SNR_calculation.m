
clear all;close all
path(path,'../CommonScripts.dir');
Fs=1000;
Nfft=128;
freq=linspace(0,Fs,Nfft);

%%Parameters for noise
Tlen=3;
x=randn(Tlen*Fs,1);  %i.e. we will pad with

param.image.threshold.eq_time=1;
x_rms{1}=sqrt(mean(x(1:(param.image.threshold.eq_time*Fs)).^2));

%%%%%%%%%%RMS derived from power spectral density%%%%%%%%
[S,F,T,Pxx_spec]=spectrogram(x,Nfft,0,Nfft,Fs);
Imedian=find(T<=param.image.threshold.eq_time);
Pxx_spec.noise=mean(Pxx_spec(:,Imedian)');
x_rms{2}=sqrt(trapz(Pxx_spec.noise)*(Fs/Nfft))

%%%%%%%%%%%




% %Create signal...
SNR=1000;
fc=100; %Hz
param.sig_time=2;
t=(1:round(param.sig_time*Fs))/Fs;
xs=sqrt(SNR)*cos(2*pi*fc*t).*hanning(length(t)).';
%
t=(1:round(Tlen*Fs))/Fs;
index=round(param.image.threshold.eq_time*Fs)+(1:length(xs));
x(index)=x(index)+xs.';
xs_rms{1}=sqrt(mean(x(index).^2))
disp(sprintf('SNR is %6.2f dB',20*log10(xs_rms{1}/x_rms{1})));

%[S,F,T,Pxx_spec]=spectrogram(x,Nfft,round(0.9*Nfft),Nfft,Fs);
[S,F,T,PP]=spectrogram(x,Nfft,0,Nfft,Fs);

Imedian=find(T>param.image.threshold.eq_time&T<param.image.threshold.eq_time+param.sig_time);
Pxx_spec.signal=mean(PP(:,Imedian)');
xs_rms{2}=sqrt(trapz(Pxx_spec.signal)*(Fs/Nfft))

%Now try the morph approach..
B=10*log10(abs(PP)./x_rms{1});
Imedian2=find(B>0);
xs_rms{3}=sqrt((Fs/Nfft)*sum(PP(Imedian2))/length(Imedian))
param=TOC_params('April21_2008_initial_both');
param.image.threshold.eq_time=1;
[stats,BWfinal,Bmean]=contour_postprocessor_morph(x,0,param,1,[]);