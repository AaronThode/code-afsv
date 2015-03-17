function [x_sel,y_min,y_max] = sourceSelect_dB(x,N,NFFT,Nwindow)

% Time selection
tfr=tfrstft(x,1:N,NFFT,hamming(Nwindow));
TFR=abs(tfr);

fig=imagescFun(1:N,(1:NFFT)/NFFT,10*log10(TFR),[0 1/2],...
    'Click twice on the spectrogram to zoom on the signal');

disp('Selecting time-frequency window')
lims=ginput(2);
close(fig)

[y_min,y_max]=deal(min(lims(:,2)),max(lims(:,2))); % min and max are used so that the order of the clicks doesn't matter
[i_min,i_max]=deal(max(1,floor(y_min*N)),min(N,floor(y_max*N)));
x_f=fft(x,N);
x_f(1:i_min)=0;
x_f(i_max:end)=0;
x=ifft(x_f,N);

[j_min,j_max]=deal(max(1,round(min(lims(:,1)))),min(N,round(max(lims(:,1))))); % min and max are used so that the order of the clicks doesn't matter
x_sel=x(j_min:j_max);

disp('End of time-frequency selection')
end