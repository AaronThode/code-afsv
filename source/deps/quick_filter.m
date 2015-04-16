%function [y,B]=quick_filter(x,Fs,minfreq,maxfreq)
%%minfreq and maxfreq in Hz
function [y,B]=quick_filter(x,Fs,minfreq,maxfreq,B)
df=maxfreq-minfreq;
if ~exist('B','var')
    frange=[max([0 minfreq-0.2*df]) minfreq maxfreq min([maxfreq+0.2*df Fs/2])];
    [N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.1],Fs);
    B = firpm(N,Fo,Ao,W);
end
y=filtfilt(B,1,x-mean(x));

%freqz(B,1,linspace(0,Fs/2,128),Fs);