function [] = Whale_multi_hydro(filt,x,Fs,N_window,NFFT,r,c)

% Signal selection
N=length(x);
x_f=fft(x,N);
x_f(1:N*filt.lims(3))=0;
x_f(N*filt.lims(4):end)=0;
x=ifft(x_f,N);
x=x(filt.lims(1):filt.lims(2));

tfr=tfrstft(x,1:N,NFFT,hamming(N_window));
TFR=abs(tfr);
imagescFun(1:N,(1:NFFT)/NFFT,10*log10(TFR),[0 1/2],'Selection');

%Decimation
decim_fact=ceil(1/2.1/filt.lims(4));
x=decimate(x,decim_fact);
N_dec=length(x);

%Deconvolution
source_phase=angle(fft(filt.source,N_dec));
x=ifft(fft(x,N_dec).*exp(-1i*source_phase));

sig_whale=abs(tfrstft(x,1:N_dec,NFFT,hamming(N_window)));
imagescFun(1:N_dec,(1:NFFT)/NFFT,10*log10(sig_whale),filt.lims(3:4),'Signal after source deconvolution');


