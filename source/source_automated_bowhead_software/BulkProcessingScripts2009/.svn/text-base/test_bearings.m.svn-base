%%%test_bearings.m%%
function test_bearings

param.Nfft=128;
param.ovlap=7/8;

filedir='../CommonScripts.dir/';
fnames=dir([filedir 'sound*mat']);

for I=1:length(fnames),
    name{I}=fnames(I).name;
end
Ichc=menu('Select a file:',name);

data=load([filedir  name{Ichc}]);
Fs=data.Fs;

plot_spectrogram(data.x);
disp('Select fmin and fmax:')
tmp=ginput(2);
fmin=tmp(1,2);
fmax=tmp(2,2);

[thet,x]=extract_bearings(data.x,0,param.Nfft,Fs,fmin,fmax,'CSDM');
disp(sprintf('uncorrected bearing is %6.2f',thet));

function plot_spectrogram(xx)
        figure(1)
        for Ishow=1:3,
            [S,F,T,PP1]=spectrogram(xx(:,Ishow),param.Nfft,round(param.Nfft*param.ovlap),param.Nfft,Fs,'yaxis');
            set(gcf,'units','norm','pos',[6.666e-02     4.858333e-01     2.91666e-01     3.5000e-01]);
            subplot(3,1,Ishow);
            imagesc(T,F,10*log10(abs(PP1)));axis('xy');
        end
        %pause
end

end