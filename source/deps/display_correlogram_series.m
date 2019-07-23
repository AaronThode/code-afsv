function display_correlogram_series(handles,x,Fs,hdr,ovlap,Nfft)

fmax=1000*str2double(get(handles.edit_fmax,'String'));
fmin=1000*str2double(get(handles.edit_fmin,'String'));
param.ovlap=ovlap;
param.Nfft=Nfft;
prompt1={'ICI range (s)','Correlation sample time (s)','Teager-Kaiser treatment?','Incoherent=1, Coherent=0 ?'};
dlgTitle1='Parameters for correlogram...';
def1={'[0.01 .25]', '0.25','0','1'};
answer=inputdlg(prompt1,dlgTitle1,1,def1);
param.ici_range=eval(answer{1});
param.time_sample=str2double(answer{2});
param.teager=str2double(answer{3});
alg_chc=str2double(answer{4});

if alg_chc==1
    contents=get(handles.popupmenu_WindowSize,'String');
    Nfft_window=str2double(contents{get(handles.popupmenu_WindowSize,'Value')});
    
    [S,FF,TT,B] = spectrogram(x(:,1),hanning(Nfft_window),round(ovlap*Nfft_window),Nfft,Fs);
    [mean_corr_org,tindex,TT_plot,pwr,pwr_tot]= create_incoherent_correlogram(TT,FF,B,param,fmin,fmax);
    
    
else
    param.Nfft=round(param.time_sample*Fs);
    [mean_corr_eq,mean_corr_org,tindex,TT_plot,pwr]= create_coherent_correlogram(x(:,1),Fs,param,fmin,fmax);
    %
    
end

dX=tindex(2)-tindex(1);  %X axis of new correlation image
Ntime=length(tindex)-1;
imagesc(tindex(1:(end-1)),TT_plot,mean_corr_org);caxis([0 0.3]);

%colorbar('east','ycolor','w');
axis('xy')
% title(sprintf('mean correlation value extracted between %6.2f and %6.2f Hz, %6.2f overlap',freq(Ifreq),freq(Ifreq+1),param.ovlap));
climm(1)=str2double(get(handles.edit_mindB,'String'));
climm(2)=climm(1)+str2double(get(handles.edit_dBspread,'String'));

if climm(1)>1
    climm(1)=0;
    set(handles.edit_mindB,'String','0');
end
if climm(2)>1
    climm(2)=1;
    set(handles.edit_dBspread,'String','1');
end
caxis(climm);

if alg_chc==0  %coherent Correlogram
    caxis([-20 0]);
    set(handles.edit_mindB,'String','-20');
    set(handles.edit_dBspread,'String','20');
end

if get(handles.checkbox_grayscale,'Value')==1,
    colormap(flipud(gray));
else
    colormap(jet);
end
colorbar;
% set(gcf,'pos',[30   322  1229   426])
set(gca,'fontweight','bold','fontsize',14);
xlabel('Time (sec)');ylabel('Correlation lag (sec)');

end