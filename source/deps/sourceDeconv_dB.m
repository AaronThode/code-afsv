function [x_deconv,source] = sourceDeconv_dB(x,Fs,NFFT,Nwindow,f_min,f_max,Npoint)

N=length(x);
% plot the tf response to select the points and interpolate the 1st mode
sig_whale=abs(tfrstft(x,1:N,NFFT,hamming(Nwindow)));

% Selecting points
disp('Beginning source deconvolution')
fig=imagescFun((1:N)/Fs,Fs*(1:NFFT)/NFFT,10*log10(sig_whale),[f_min f_max],'Choose the points to interpolate the 1st mode');
caxis([prctile(prctile(10*log10(sig_whale),15),15) prctile(prctile(10*log10(sig_whale),95),95)])
[t_points,f_points]=ginput(Npoint);

[~,sortIndex]=sort(t_points);
[t_points,f_points]=deal(t_points(sortIndex),f_points(sortIndex));
close(fig)

% Linear interpolation between the selected points

t_interp=linspace(t_points(1),t_points(end),round(1+Fs*(t_points(end)-t_points(1))));
iflaw=interp1(t_points,f_points,t_interp,'linear')';

% Deconvolution (source phase cancelation)
x_f=fft(x,N);
source=fmodany(iflaw/Fs);
source_phase=angle(fft(source,N));

x_deconv=ifft(x_f.*exp(-1i*source_phase));

disp('End of source deconvolution')

end

