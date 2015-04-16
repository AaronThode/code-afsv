%function [y,ysignal,ynoise]=simulate_FMsweep_with_noise(Fs,Nfft,ovlap,fstart,fend,total_window_time,tduration,tstart,SNRdB,SNRchc, no_noise);
% Create an FM sweep in white noise
% Input:
%   Fs: sampling rate in Hz
%   Nfft: FFT size used for computing spectrogram (not really needed).
%   fstart: start frequency of sweep in Hz
%   fend: end freqency of sweep in Hz
%   total_wind_time: total length of time series in seconds
%   tduration: duration of sweep in seconds
%   tstart: time at which sweep begins in output time series
%   SNRdB: desired SEL SNR in dB
%   SNRchc: 'SEL' or 'PSD'- how should SNR be defined?
%   no_noise:  if exists, don't add any noise to signal regardless of value
%       of SNRdB.
% Output:
%   y: export time series, sampled at Fs Hz

function [y,ysignal,ynoise, noise_var]=simulate_FMsweep_with_noise(Fs,Nfft,ovlap,fstart,fend,total_window_time,tduration,tstart,SNRdB,SNRchc,no_noise)

ysignal=[];ynoise=[];y=[];

slope=(fend-fstart)/tduration;
t=linspace(0,tduration,round( Fs*tduration));
dt=t(2)-t(1);
y0=cos(2*pi*(fstart*t+0.5*slope*t.^2));
y=zeros(round(total_window_time*Fs),1);
y(round(tstart* Fs)+(1:length(t)))=y0;
ysignal=y;

SEL=trapz(y0.^2)*dt;
Power=SEL/tduration;

disp(sprintf('Narrowband SEL in time domain: %12.8f, Narrowband Power in time domain: %12.8f', ...
    SEL,Power));

if exist('no_noise'),
    disp('No noise added in simulate_FMsweep_with noise.');
    y=ysignal;
    return
end
SNR=10^(SNRdB/10);

if strcmp(SNRchc,'SEL')
    noise_var=0.5*(Fs/abs(fend-fstart))*(SEL/(SNR*tduration));
else
    noise_var=0.5*(Fs/abs(fend-fstart))*(Power/SNR);
end
ynoise=sqrt(noise_var)*randn(length(y),1);

Power_noise=dt*trapz(ynoise.^2)/total_window_time;
Power_noise=Power_noise*(abs(fend-fstart)/(0.5*Fs));

disp(sprintf('Narrowband Noise power in time domain: %12.8f, SNR: %12.8f', ...
    Power_noise,Power/Power_noise));


y=ysignal+ynoise;


end
