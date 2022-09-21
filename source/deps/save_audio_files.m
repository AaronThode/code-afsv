tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;
try
    chan=eval(get(handles.edit_chan,'String'));
    if strcmpi(handles.filetype,'GSI')
        choice=menu('Store all channels?','Yes','No');
        if choice==1
            chan=1:3;
        else
            chan=-1;
        end
    end
    if chan<0
        Ichan=chan;
    else
        Ichan='all';
    end
catch
    Ichan='all';
end


%[x,t,Fs,tstart,junk,hdr]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,Ichan,handles);
[x,~,Fs,~,~,hdr]=load_data(handles.filetype,tdate_start,tlen,Ichan,handles);

Islash=strfind(handles.mydir,filesep);
chan=get(handles.edit_chan,'String');

%save_name=sprintf('soundsamp%s_%s',handles.myfile((Islash(end)+1):end),datestr(tdate_start,30));
save_name_wav=sprintf('soundsamp%s_%s_%3.1fx',handles.myfile,datestr(tdate_start,30),handles.audio_scale);
save_name_mat=sprintf('soundsamp%s_%s',handles.myfile,datestr(tdate_start,30));

disp(['Saving ...' save_name_wav]);

save_path	=	fullfile(handles.outputdir, save_name_wav );  %AARON: save to local directory, not server
save_path_mat=fullfile(handles.outputdir, save_name_mat);
try
    if size(x,2)>size(x,1)
        x=x';
    end
    for Ichan=1:size(x,2)
        xfilt(:,Ichan)=filter(handles.b,1,x(:,Ichan));
    end
    audiowrite(save_path,xfilt/(1.1*max(max(abs(xfilt)))),round(Fs));
    
catch
    disp('No filtering desired... exists; saving raw acoustic data to WAV');
    xfilt=[];
    matVers=version('-release');
    
    if(str2double(matVers(1:4))>2012 || strcmp(matVers,'2012b'))
        audiowrite([save_path '.wav'],(x/(1.1*max(max(abs(x))))),round(Fs*handles.audio_scale));
        % audiowrite([save_path '_32d.wav'],(x/(1.1*max(max(abs(x))))),round(Fs/32));
        
        %audiowrite([save_path '.wav'],(x/(1.1*max(max(abs(x))))),round(Fs));
    else
        %        wavwrite(x/(1.1*max(max(abs(x)))),round(Fs),[save_path '.wav']);
        audiowrite([save_path '.wav'], x/(1.1*max(max(abs(x)))),round(Fs));
    end
    
end
save([save_path_mat '.mat'],'x','xfilt','Fs','tdate_start','hdr');