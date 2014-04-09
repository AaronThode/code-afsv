%function [SEL,power]=extract_SEL_n_PSD(x,Fs,frange)
% Input: 
%  x: signal of interest
%  Nfft: FFT window length to apply to x
%  Fs: sampling rate in Hz
%  frange: two-element vector of upper and lower frequency desired in Hz
%  
%  Output:
% SEL.broadband=dt*trapz(x.^2);
% SEL.broadband_fft=df*sum(abs(Xa).^2);
% SEL.narrowband=2*df*sum(abs(Xa(Ifrange)).^2);  %Note factor of two for two-sided psd
%
% power.broadband=SEL.broadband/Tduration;
% power.narrowband=SEL.narrowband/Tduration;
function [SEL,power]=extract_SEL_n_PSD(x,Nfft,Fs,frange,Idebug)

dt=1/Fs;
Tduration=dt*length(x);
%Nfft=2^(ceil(log10(length(x))/log10(2)));
df=Fs/Nfft;

%%IMPORTANT!  DEMEAN THE SIGNAL WHEN APPENDING ZEROS!
x=x-mean(x);
Xa=dt*fft(x,Nfft);  %Analytic Fourier transform

% figure;
% plot(linspace(0,Fs,Nfft),20*log10(abs(Xa)));

Ifrange=round(frange*Nfft/Fs);
if Ifrange(1)==0,Ifrange(1)=1;end
Ifrange=Ifrange(1):Ifrange(2);

SEL.broadband=dt*trapz(x.^2);
SEL.broadband_fft=df*sum(abs(Xa).^2);
SEL.narrowband=2*df*sum(abs(Xa(Ifrange)).^2);  %Note factor of two for two-sided psd

power.broadband=SEL.broadband/Tduration;
power.narrowband=SEL.narrowband/Tduration;

if exist('Idebug')
    plotstr='brk';
    figure(12);
    F=linspace(0,Fs,Nfft);
    %subplot(2,1,1);
    plot(F,20*log10(abs(Xa)),plotstr(Idebug));grid on;hold on
    ylim([-50 50]);
    line(frange(1)*[1 1],[-50 50],'Color',[0 0 0],'linewidth',3);
    line(frange(2)*[1 1],[-50 50],'Color',[0 0 0],'linewidth',3);
    
    %title('SEL')
    %subplot(2,1,2);
    %plot(F,10*log10(abs(SEL.narrowband)),plotstr(Idebug));grid on;hold on
    
    disp(sprintf('Broadband SEL in time domain: %12.8f, Broadband SEL in freq domain: %12.8f, narrowband SEL: %12.8f', ...
        SEL.broadband,SEL.broadband_fft,SEL.narrowband));

    disp(sprintf('Broadband power in time domain: %12.8f, Narrowband power in freq domain: %12.8f', ...
        power.broadband,power.narrowband));
    
   
end



%
%
% noise_SEL=trapz(noise_y(1:length(t)).^2)*dt;
% noise_Power=noise_SEL/max(t);
%
% %ovlap=0;
% [S,F,T,PSD]=spectrogram(y, Nfft,round(ovlap*Nfft), Nfft, Fs,'yaxis');
%
% Igood=find(T<=max(t));  %Take length of time equal to signal
%
% if (max(t)>tstart),
%     disp('Warning, your freq-domain noise esstimate is contaminated by signal--increase tstart');
% end
% %%Get frequency range of signal.
% I1=floor(fstart*Nfft/Fs);
% I2=floor(fend*Nfft/Fs);
% %noise_SEL2=sum(sum(PSD(I1:I2,Igood)));
%
% noise_SEL2=sum(PSD(I1:I2,Igood(1)))+(1-ovlap)*sum(sum(PSD(I1:I2,Igood(2:end))));
% noise_Power2=(Fs/Nfft)*noise_SEL2/length(Igood);
%
% disp('******************');
% disp(sprintf('%\n SEL of signal, time-domain noise SEL, freq-domain noise SEL'));
% disp(' The first and last quantities should be equal');
% [SEL noise_SEL noise_SEL2]
%
% disp('******************');
% disp(sprintf('%\n SEL of signal, time-domain noise SEL, freq-domain noise SEL'));
% disp(' The first and last quantities should be equal');
% [Power noise_Power noise_Power2]
%
%
% imagesc(T,F,10*log10(abs(PSD)));colorbar
