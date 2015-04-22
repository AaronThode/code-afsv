function [x_filt] = Whale_multi_hydro(filt,x)

% Signal selection
x=x(filt.xlims(1):filt.xlims(2));

%Deconvolution
Nx=length(x);
source_phase=angle(fft(filt.source,Nx));
x_filt=ifft(fft(x,Nx).*exp(-1i*source_phase));