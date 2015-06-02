function [x_filt] = Whale_multi_hydro(filt,x)

Nx=length(x);
Fsource=fft(filt.source,Nx);
Fsource(Nx/2+2:end)=conj(Fsource(Nx/2:-1:2));
source_phase=angle(Fsource);

x_filt=ifft(fft(x,Nx).*exp(-1i*source_phase+1i*2*pi*filt.xmin*(0:Nx-1)'/Nx));