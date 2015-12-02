%%%%%%%%%%plot_time_series.m%%%%%%%%%%%
% given TL, compute time series and plot results..
%
% Inputs:
%       tilt: horizontal offset between shallowest and deepest element in meters
%
% function [x,Tmin]=plot_time_series(TL,rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,xsource,tilt)

function [x,Tmin]=plot_time_series_all(TL,rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,xsource,tilt,tlimm)

rd=zplot;
tiltt=tilt*(rd-min(rd))./(max(rd)-min(rd));
kk=2*pi*freq/1495;
index=round(fs*tlimm(1)):round(fs*tlimm(2));  

x=zeros(length(zplot),length(rplot),length(index));
for Iazi=1:length(azi)  %Azimuthal
    for Ir=1:length(rplot)
        for Iz=1:length(zplot)
            Xfft=zeros(1,Nsamp);
            extrct=squeeze(TL(Iz,Ir,:,Iazi));
            
            x0=zeros(1,Nsamp);
            
            %%%Create tilt matrix...
            tilt_matrix=exp(1i*tiltt(Iz)*kk);
            if size(extrct,1)==nfreq
                extrct=extrct.';
            end
            Xfft(flo:fhi)=extrct(1:nfreq).*tilt_matrix.*exp(-1i*2*pi*freq*Tmin(Ir));
            x0(1,:)=(2/Nsamp)*real(fft(Xfft,Nsamp));%Factor of two to account for negative frequencies
            
            %Trim
            %%%Include source signal
            if exist('xsource','var')&&~isempty(xsource)
                x(Iz,Ir,:)=conv(x0(1,index),xsource','full');
            else
                x(Iz,Ir,:)=x0(1,index);
                
            end
            
            
           
        end
        %Correct for trimmed time series
        Tmin(Ir)=Tmin(Ir)+tlimm(1);
        
    end
end


end %Iazi