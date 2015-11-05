%%%%%%%%%%plot_time_series.m%%%%%%%%%%%
% given TL, compute time series and plot results..
%
% Inputs:
%       tilt: horizontal offset between shallowest and deepest element in meters
%
% function x=plot_time_series(TL,sd,rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,range_normalization,xsource,tilt)

function x=plot_time_series(TL,sd,rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,range_normalization,xsource,tilt)
%
%%First, check if have a range/depth grid and force selection of a row or
%%column...



if length(rplot)>1&&length(zplot)>1
    choice=menu('Choose plot:','range slice (fixed depth)','depth slice (fixed range)');
    if choice==1, %range slice
        for II=1:length(zplot),itemlist{II}=zplot(II);end
        Iz_slice=menu('Choose depth:',itemlist);
        Ir_slice=1:length(rplot);
    else %depth slice
        for II=1:length(rplot),itemlist{II}=rplot(II);end
        Ir_slice=menu('Choose range:',itemlist);
        Iz_slice=1:length(zplot);
        
    end
elseif length(zplot)==1,
    Iz_slice=1;Ir_slice=1:length(rplot);choice=1;
else
    Ir_slice=1;Iz_slice=1:length(zplot);choice=2;
end


for Iazi=1:length(azi)
    ro=rplot(Ir_slice);
    rd=zplot(Iz_slice);
    Xfft=zeros(1,Nsamp);
    %How to plot from TL matrix?  depth example...
    if choice==2 %multiple receiver depths and only source depth
        extrct=squeeze(TL(Iz_slice,Ir_slice,:,Iazi));
        rfactor=ro(end)*ones(1,length(rd)); %geometric removal factor
        %if length(ro)>1,disp('Can''t plot multiple range and depths for time series, using value for max range');end
        nplots=length(rd);
    else choice==1 %multiple ranges
        extrct=squeeze(TL(Iz_slice,Ir_slice,:,Iazi));
        rfactor=ro;
        nplots=length(ro);
        
    end
    
    x0=zeros(nplots,Nsamp);
    %x=zeros(nplots,Nsamp);
    Nplots_total=5;
    Iplots=0;
    
    if length(rd)>1
        tiltt=tilt*(rd-min(rd))./(max(rd)-min(rd));
    else
        tiltt=0; %tilt has no meaning for a single phone...
    end
    kk=2*pi*freq/1495;
    for Ir=1:nplots
        Iplots=Iplots+1;
        if Iplots>Nplots_total
            figure
            Iplots=1;
        end
        subplot(Nplots_total,1,Iplots);
        
        
        if choice==1,  %range slice
            t_axis=Tmin(Ir)+(1:Nsamp)/fs;
            tlimm=[Tmin(Ir) max(t_axis)];
            Tminn=Tmin(Ir);
        else  %depth slice
            Tminn=Tmin(Ir_slice);
            t_axis=Tminn+(1:Nsamp)/fs;
            tlimm=[Tminn max(t_axis)];
            
        end
        %subplot(nplots,1,Ir)
        
        %%%Create tilt matrix...
        tilt_matrix=exp(1i*tiltt(Ir)*kk);
        %extrct is [nel nfreq]
        if size(extrct,1)==nfreq
           extrct=extrct.'; 
        end
        Xfft(flo:fhi)=extrct(Ir,1:nfreq).*tilt_matrix.*exp(-1i*2*pi*freq*Tminn);
        x0(Ir,1:Nsamp)=(2/Nsamp)*real(fft(Xfft,Nsamp));%Factor of two to account for negative frequencies
        
        %%%Include source signal
        
        if exist('xsource')&&~isempty(xsource)
            x(Ir,:)=conv(x0(Ir,1:Nsamp),xsource','full');
        else
            x(Ir,1:Nsamp)=x0(Ir,1:Nsamp);
            
        end
        %x(Ir,1:Nsamp)=fliplr(x(Ir,1:Nsamp));
        if range_normalization==1,
            plot(t_axis,rfactor(Ir)*x(Ir,:),'k');  %Spherical spreading normalization
            disp('pressures are normalized by range to account for spherical spreading');
        else
            %plot(t_axis,x(Ir,:)*1e-6,'k');
            %plot(t_axis,x(Ir,:),'k');
             Nfft=4*1024;
            [S,FF,TT,B] = spectrogram(x(Ir,:),hanning(Nfft/4),round(.9*Nfft/4),Nfft,fs);
            %figure('name','Pulse spectrogram');
            B=10*log10(B);
            imagesc(TT,FF/1000,B-max(max(B)));%
            grid on
            axis('xy');
            %spectrogram(x(Ir,:),hanning(2048),round(0.9*2048),2048,fs,'yaxis')
            ylim([min(freq) max(freq)]/1000);
            ylim([0 .25]);
            colormap(jet);caxis([-20 0]);xlim([4.5 5.5]);colorbar
        end
        ylabel('Hz','fontweight','bold');xlabel('Time (s)','fontweight','bold');
        set(gca,'fontweight','bold');
        if length(rd)>1 %multiple receiver depths...
            title(sprintf('Range: %6.2f km Receiver depth: %6.2f m Source depth: %6.2f m azimuth: %6.2f', ...
                ro/1000,rd(Ir),sd,azi(Iazi)));
        elseif length(ro)>1,
            title(sprintf('Range: %6.2f km Receiver depth: %6.2f m Source depth: %6.2f m azimuth: %6.2f', ...
                ro(Ir)/1000,rd,sd,azi(Iazi)));
        end
        %xlim(tlimm);
        %ylim([-200 200]);  %%%%%%%%%%%Educational note...this gives the same result...
        
        
    end
    
    %     ylimm=input('Enter frequency range if desired:');
    %     if ~isempty(ylimm),
    %         for Ir=1:nplots,
    %             subplot(nplots,1,Ir)
    %             ylim(ylimm);
    %         end
    %     end
    %
    %     xlimm=input('Enter time range if desired:');
    %     if ~isempty(xlimm),
    %         for Ir=1:nplots,
    %             subplot(nplots,1,Ir)
    %             xlim(xlimm);
    %         end
    %     end
    
end %Iazi