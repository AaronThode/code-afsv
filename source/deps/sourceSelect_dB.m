function [x_sel,f_min,f_max,j_max] = sourceSelect_dB(x,Fs,NFFT,Nwindow,flims)
% x_sel: snippet
% y_min,y_max, min and max frequencies of bounding box
% Time selection

N=length(x);
% Time selection
tfr=tfrstft(x,1:N,NFFT,hamming(Nwindow));
TFR=abs(tfr);

fig=imagescFun((1:N)/Fs,Fs*(1:NFFT)/NFFT,10*log10(TFR),flims,...
    'Click twice on the spectrogram to zoom on the signal');
caxis([prctile(prctile(10*log10(TFR),15),15) prctile(prctile(10*log10(TFR),95),95)])
disp('Selecting time-frequency window')
lims=ginput(2);
close(fig)

[f_min,f_max]=deal(min(lims(:,2)),max(lims(:,2))); % min and max are used so that the order of the clicks doesn't matter
[x,~]=quick_filter(x,Fs,f_min,f_max);


[j_min,j_max]=deal(max(1,round(Fs*min(lims(:,1)))),min(N,round(Fs*max(lims(:,1))))); % min and max are used so that the order of the clicks doesn't matter
x_sel=x(j_min:j_max);

disp('End of time-frequency selection')
end