%%%%%%%%%%%%%%%%%%%%%%%%plot_TL.m%%%%%%%%%
% plot transmission loss 
   
for Iazi=1:length(azi),
figure;
    subplot(2,1,1);
    plot(Xrot{Iazi}(:),Yrot{Iazi}(:),'o','markersize',10,'markerfacecolor',[1 0 0]);
    title(sprintf('Geometry of %s array %i azimuth',case_array,azi(Iazi)));
    xlabel('Xrot: propagation plane (m)','fontweight','bold');
    ylabel('Yrot: perpendicular plan (m)','fontweight','bold');
    grid on;axis('equal');
    %Just plot TL of single element
    
     %subplot(length(azi),1,Iazi);
     subplot(2,1,2);
     TLdB=-20*log10(abs(TL(:,:,:,Iazi))+eps);
    
    imagesc(rplot,zplot,TLdB)
    if length(rplot)~=size(TLdB,2),
        disp('ranges dont match');
    end
    if length(zplot)~=size(TLdB,1),
        disp('depths dont match')
    end
    colormap(flipud(jet));
    xlabel('range(m)','fontweight','bold');ylabel('depth (m)','fontweight','bold');
    set(gca,'fontweight','bold')
    title(num2str(rd));
    %caxis([40 120]);
     caxis([20 90]);
     colorbar;
   
    %title(sprintf('Average of %i frequencies between %i and %i Hz',length(freq),freq(1),freq(end)));
    title(sprintf('Source depth: %6.2f m, frequency range: %i-%i Hz, bottom type: %s, azimuth %i model:%s', ...
        sd,min(freq),max(freq),case_bottom', azi(Iazi),code_chc));
    
    if ~isempty(bath),
        line(bath(:,1),bath(:,2));    
    end
end
   % gtext(plottitle);

   %MOV(Iazi)=getframe;
   
   
