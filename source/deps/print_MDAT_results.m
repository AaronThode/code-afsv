
mychc=menu('Which print option?','Single phone','all phones');

if mychc==2  %%Multiphone data...
    
    figure;
    tdate_start=handles.tdate_start;
    tlen=handles.tlen;
    contents=get(handles.popupmenu_Nfft,'String');
    Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});
    
    contents=get(handles.popupmenu_ovlap,'String');
    ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
    ovlap=min([1-1/Nfft ovlap]);
    
    mydir=pwd;
    figure(1)
    set(gcf,'pos',[291         628        1513         991]);
    Iplot=0;
    [x,~,Fs,~,~,head]=load_data(handles.filetype, tdate_start,tlen,'all',handles);
    
    for Ichan=1:head.Nchan
        
        if Iplot<8
            Iplot=Iplot+1;
        else
            figure(2)
            set(gcf,'pos',[291         628        1513         991]);
            Iplot=1;
            
        end
        subplot(4,2,Iplot);
        
        contents=get(handles.popupmenu_WindowSize,'String');
        Nfft_window=str2double(contents{get(handles.popupmenu_WindowSize,'Value')});
        
        [~,FF,TT,B] = spectrogram(x(Ichan,:),hanning(Nfft_window),round(ovlap*Nfft_window),Nfft,Fs);
        %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
        %axes(handles.axes1);
        imagesc(TT,FF,10*log10(B));%
        %xlimm=str2double(get(handles.))
        axis('xy')
        fmax=str2double(get(handles.edit_fmax,'String'));
        fmin=str2double(get(handles.edit_fmin,'String'));
        if fmax==0
            ylim([0 Fs/2]);
            set(handles.edit_fmax,'String',num2str(Fs/2));
        else
            ylim([fmin fmax]*1000);
        end
        %ylim([0 1]);axis('xy')
        climm(1)=str2double(get(handles.edit_mindB,'String'));
        climm(2)=climm(1)+str2double(get(handles.edit_dBspread,'String'));
        caxis(climm);
        if get(handles.checkbox_grayscale,'Value')==1
            colormap(flipud(gray));
        else
            colormap(jet);
        end
        % set(gcf,'pos',[30   322  1229   426])
        set(gca,'fontweight','bold','fontsize',14);
        ylabel('Hz');
        if rem(Ichan,2)==0
            ylabel('');
        end
        if Iplot<7
            set(gca,'xticklabel',[]);
        else
            xlabel('Time (sec)');
        end
        poss=get(gca,'pos');
        poss(4)=.2;poss(3)=.4;
        set(gca,'pos',poss);
        set(gca,'yminorgrid','on','xminorgrid','on');
        try
            text(0.1,0.9,num2str(head.geom.rd(Ichan),3),'color','y','units','norm','fontsize',14,'fontweight','bold');
        catch
            disp('No elements depths provided');
        end
        hh=colorbar('East');
        set(hh,'fontweight','bold','fontsize',14)
        
        %             if Ichan==1
        %                 xall=zeros(8,length(x));
        %             end
        %xall(Ichan,:)=x;
    end
    if ~isempty(strfind(lower(computer),'mac'))
        Islash=strfind(handles.mydir,'/')+1;
    else
        Islash=strfind(handles.mydir,'\')+1;
    end
    save_name=sprintf('soundNfft%i_%s.%s_%s_%s',Nfft,handles.mydir((Islash(end-1)):(end-1)),handles.myfile,datestr(tdate_start,30),'all_channels');
    pause
    
    for I=1:2
        if any(get(0,'child')==I)
            save_name1=sprintf('%s_%i',save_name,I);
            save_path	=	fullfile(handles.outputdir, save_name1);
            disp(['Printing %s ...' save_name1]);
            figure(I)
            orient landscape
            print(I,'-djpeg','-r500',[save_path '.jpg']);
            save([save_path '.mat'],'x','Fs');
            close(I);
        end
        
    end
    
    return
end  %if mychc==2