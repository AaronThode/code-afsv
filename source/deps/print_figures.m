if strcmpi(handles.filetype,'MDAT')
    print_MDAT_results;
end  %MDAT

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');
OwnFig=handles.togglebutton_NewFigure.Value;

if ~OwnFig
    figchc=[];
    axes(handles.axes1);
else
    figure(gcf);
    chcc=get(0,'child');
    chcc=[chcc.Number];
    if isempty(chcc)
        return
    end
    Igoodd	=	(chcc-round(chcc)==0);
    chcc=chcc(Igoodd);
    if isempty(chcc)
        return
    end
    for III=1:length(chcc) %#ok<*FORPF>
        tmp{III}=chcc(III);
    end
    figchc=menu('Select a figure number:',tmp);
    figure(chcc(figchc));
end
tdate_start=handles.tdate_start;
tlen=handles.tlen;
chan=get(handles.edit_chan,'String');
Idot=strfind(handles.myfile,'.')-1;

if OwnFig
    strr='NewFig';
else
    strr='Screen';
end
save_name=sprintf('soundsamp_%s_%s_%s_%s_%s', ...
    handles.myfile(1:Idot),datestr(tdate_start,30),chan,handles.display_view,strr);

save_path	=	fullfile(handles.outputdir, save_name);
disp(['Printing %s ...' save_name]);

%orient landscape
%print(gcf,'-djpeg','-r300',sprintf('%s.jpg',save_path));
if ~OwnFig
    exportapp(gcf,sprintf('%s.jpg',save_path))
else
    set(gcf,'pos',[136         231        1255         766]);
    print(gcf,'-djpeg','-r300',sprintf('%s.jpg',save_path));
end