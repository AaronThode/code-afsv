function handles = load_and_display_spectrogram(app, handles)

cla;

%What type of Display?
handles.old_display_view=[];
if isfield(handles,'display_view')
    handles.old_display_view=handles.display_view;
end

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');


%%Read off parameters from figure...

%tdate_start		=	handles.tdate_start;
if isfield(handles,'tlen')
    tlen	=	handles.tlen;
else
    tlen=[];
end
%if tlen<1e-2 %somehow window length is way too short
%    tlen=1;
%   set(handles.edit_winlen,'String',num2str(tlen));
%end
if isempty(tlen)
    tlen    =   str2double(get(handles.edit_winlen,'String'));
    
    handles.tlen=tlen;
end
mydir	=	pwd;

% [~,~,~,~,~,hdr]=load_data(handles.filetype, ...
%        handles.tdate_start,tlen,Ichan,handles);

want_directionality=strcmp(handles.display_view,'Directionality')||strcmpi(handles.display_view,'KEtoPERatio');
want_directionality=want_directionality||strcmpi(handles.display_view,'ItoERatio')||strcmpi(handles.display_view,'IntensityPhase');
%want_directionality=want_directionality || handles
if want_directionality&strcmpi(handles.filetype,'gsi')
    Ichan='all';
elseif want_directionality&contains(handles.myfile,'drifter')
    Ichan=9:11;
else
    Ichan	=	eval(get(handles.edit_chan,'String'));
end

try
    % [x,t,Fs,tmin,tmax,head]	=	...
    %load_data(filetype,tdate_start,tlen,Ichan,handles)
    
    if ~exist(fullfile(handles.mydir,handles.myfile),'file')
        fprintf('%s does not exist! \n',fullfile(handles.mydir,handles.myfile));
        return
    end
    [x,t,Fs,tstart,junk,hdr]=load_data(handles.filetype, ...
        handles.tdate_start,tlen,Ichan,handles);
    handles.file_flags.multichannel=hdr.multichannel;
    handles.file_flags.linked=hdr.linked;
    handles.file_flags.vector_sensor=hdr.vector_sensor;
    
    %%Change file display if a transformation of a basic file has
    %%  occurred...
    if isfield(hdr,'myfile')  %file has been processed (e.g. an accelerometer file)
        handles.myfile=hdr.myfile;
        additional_text=sprintf(hdr.transformation.description{1},hdr.transformation.value{1});
    else
        additional_text=[];
    end
    set(handles.text_filename,'String',[fullfile(handles.mydir, handles.myfile) ' ' additional_text]);
    
catch
    uiwait(errordlg('Cannot load spectrogram: perhaps event or time desired too close to edge'));
    return
end

%%Imagesc PSD file

if strcmpi(handles.filetype,'psd')
    FF=(0:hdr.Nmax)*hdr.Fs/hdr.Nfft;
    set(handles.edit_maxfreq,'String',num2str(hdr.Fs/2));
    set(handles.edit_fmax,'String',num2str(hdr.Fs/2000));
    
    contents=get(handles.popupmenu_Nfft,'String');
    Value=find(strcmp(contents,int2str(hdr.Nfft))>0);
    set(handles.popupmenu_Nfft,'Value',Value);
    
    ovlap=100*(1-hdr.dn/hdr.Nfft);
    Nfft=hdr.Nfft;
    contents=get(handles.popupmenu_ovlap,'String');
    Value=find(strcmp(contents,int2str(ovlap))>0);
    set(handles.popupmenu_ovlap,'Value',Value);
    
    set(handles.edit_chan,'String','1');
    
else
    
    if isempty(x)
        errordlg('Cannot load spectrogram: perhaps event or time desired too close to edge');
        return
    end
    if max(t)<=tlen
        tlen	=	max(t);
        handles.tlen	=	max(t);
        set(handles.edit_winlen,'String',num2str(tlen));
    end
    
    if size(x,2)>1
        x=x';
    end
    
    Fs=round(Fs);
    
    mymaxfreq=str2double(get(handles.edit_maxfreq,'String'));
    
    if mymaxfreq==0||mymaxfreq>Fs/2
        handles.filter.f_max	=	Fs/2;
        set(handles.edit_maxfreq,'String',num2str(Fs/2));
    end
    %disp(sprintf('Fs=%i',Fs));
    contents=get(handles.popupmenu_Nfft,'String');
    Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});
    
    contents=get(handles.popupmenu_WindowSize,'String');
    Nfft_window=str2double(contents{get(handles.popupmenu_WindowSize,'Value')});
    
    contents=get(handles.popupmenu_ovlap,'String');
    ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
    ovlap=min([1-1/Nfft ovlap]);
    
    
end  %if PSD

%%%%Select different views
if ~handles.togglebutton_NewFigure.Value
    axes(handles.axes1);
else
    figure;
end


if want_directionality
    
    if isfield(handles,'azigram')
        azigram_param=handles.azigram;
    end
    
    if strcmpi(handles.filetype,'gsi')
        azigram_param.brefa=hdr.brefa;
    else
        azigram_param.brefa=0;
    end
    reactive_flag=handles.checkbox_reactive.Value;
    
    %%AARON test if sample delay between pressure and P.V r
    %responsible for pressure/velocity phase shift with frequency
    
    %Ishift=0;
    %x(1,:)=[x(1,(1+Ishift):end) zeros(1,Ishift)];
    %x(2:3,:)=[x(2:3,(1+Ishift):end) zeros(2,Ishift)];
    
    
    
    %%%%Wavelet test%%%%%%%%%%%%%
    
    if get(handles.checkbox_wavelet,'Value')==1
        use_wavelets=true;
    else
        use_wavelets=false;
    end
    if use_wavelets
        [TT,FF,output_array,PdB,azigram_param]=compute_directional_metrics_wavelet ...
            (x,handles.display_view,Fs,Nfft,ovlap,azigram_param, ...
            handles.filetype,reactive_flag);
        if strcmpi(handles.filetype,'gsi')&&~reactive_flag
            params.f_transition=300;
            [~,Icut]=min(abs(FF-params.f_transition));
            [~,~,temp]=compute_directional_metrics_wavelet ...
                (x,handles.display_view,Fs,Nfft,ovlap,azigram_param, ...
                'gsi',true);
            output_array{1}((Icut+1):end,:)=temp{1}((Icut+1):end,:);
        end
    else
        %%%%Here is directional metric calculation...
        %handles.display_view='PhaseSpeed';
        [TT,FF,output_array,PdB,azigram_param]=compute_directional_metrics ...
            (x,handles.display_view,Fs,Nfft,ovlap,azigram_param, ...
            handles.filetype,reactive_flag);
        %if strcmpi(handles.filetype,'gsi')&&~reactive_flag
        %params.f_transition=300;
        %[~,Icut]=min(abs(FF-params.f_transition));
        %[~,~,temp]=compute_directional_metrics ...
        %    (x,handles.display_view,Fs,Nfft,ovlap,azigram_param, ...
        %    'gsi',true);
        %output_array{1}((Icut+1):end,:)=temp{1}((Icut+1):end,:);
        %end
    end
    
    if ~contains(handles.display_view,'Directionality')
        myfig=gcf;
        myax=gca;
        switch handles.display_view
            case 'ItoERatio'
                ylimm=[0 1];
                fignum=3;
            case 'KEtoPERatio'
                ylimm=[0 10];
                fignum=4;
            case 'IntensityPhase'
                ylimm=[0 90];
                fignum=2;
            case 'PhaseSpeed'
                ylimm=[0 300];
                fignum=5;
            otherwise
                ylimm=[0 90];
                figunum=10;
        end
        figure(fignum);hold on
        hh=plot(FF,prctile(output_array{1}',50));
        %plot(FF,prctile(output_array{1}',25),'--',FF,prctile(output_array{1}',75),'--');
        %plot(FF,mean(output_array{1}'),'r',FF,mean(output_array{1}')+std(output_array{1}'),'r--',FF,mean(output_array{1}')-std(output_array{1}'),'r--');grid on;hold on
        xlabel('Frequency (Hz)');ylabel(handles.display_view);grid on
        ylim(ylimm)
        
        
        figure(10);hold on
        plot(FF,prctile(PdB',50),FF,prctile(PdB',25),'--',FF,prctile(PdB',75),'--');
        xlabel('Frequency (Hz)');ylabel('PSD dB re 1 uPa^2/Hz');grid on
        
        figure(myfig);
        axes(gca);
    else
        myfig=gcf;
        myax=gca;
        figure(1);hold on
        plot(FF,bnorm(180*circ_median(output_array{1}'*pi/180)/pi),'k.', ...
            FF,bnorm((180/pi)*(circ_mean(output_array{1}'*pi/180)+circ_std(output_array{1}'*pi/180))),'r.', ...
            FF,bnorm((180/pi)*(circ_mean(output_array{1}'*pi/180)-circ_std(output_array{1}'*pi/180))),'g.');grid on;
        xlabel('Frequency (Hz)');ylabel(handles.display_view);
        %ylim([0 360])
        figure(myfig);
        axes(gca);
    end
    
    plot_directional_metric(TT,FF,output_array{1},handles,azigram_param,PdB,use_wavelets);
    % To recover matrix use handles.axes1.Children.CData;
    %handles.azigram.azi=azi;
    handles.azigram=azigram_param;
    handles.azigram.TT=TT;
    handles.azigram.FF=FF;
    
    
    
    %%%%Take difference of azigrams
    %%FOR LUDOVIC TENORIO
    if str2num(handles.edit_chan.String)>1&strcmpi(handles.filetype,'gsi')
        
        azigram_param.sec_avg='0';
        handles1=handles;
        DASAR_ltr(1)=handles1.mydir(end-2);
        DASAR_ltr(2)=char(double(DASAR_ltr(1))+1);
        handles1.mydir(end-2)=DASAR_ltr(2);  %%%Increment DASAR letter
        handles1.myfile(5)=DASAR_ltr(2);
        
        [xall{2},~,~,~,~,hdr1]=load_data('gsi', handles1.tdate_start,tlen,'all',handles1);
        
        %Package everything together
        xall{2}=xall{2}';
        xall{1}=x;
        %azigram_param.brefa1=hdr1.brefa;
        azigram_param.brefa(2)=hdr1.brefa;
        azigram_param.goodName{1}=handles.myfile;
        azigram_param.goodName{2}=handles1.myfile;
        
        clear handles1
        
        %%%%Parameters for the automated processing
        azigram_param.tabs_start=handles.tdate_start;
        
        azigram_param.da=60;  %width of a sector to search azigram in degrees
        
        %Azimuthal grid with 50%% overlap
        azigram_param.a_grid=azigram_param.da/2:(azigram_param.da/2):(360-azigram_param.da/2);  %Azimuthal grid
        
        %%Frequency range to serach for pulses
        azigram_param.frange=[100 450];
        
        %%Minimum time/bandwidth product (e.g. pixel size) that is allowed
        %%to remain in binary azigram
        azigram_param.min_TBW=2.0480;
        
        %%%Maximum bandwidth gap permitted between pulse components
        azigram_param.vertical_line=12;  %Hz
        %params.threshold=5*4;  %Hz
        
        %%%Minimum total bandwidth required to be registered as a pulse
        %%%(doesn't have to be a single connected oval
        azigram_param.threshold=15*4;  %Hz
        
        %%%This parameter fixes a probelm with DASAR sensor.  Normally
        %%%don't adjust
        azigram_param.f_transition=300;  %Frequency at which pressure/velocity go out of phase
        
        %%%Set these to true to see the morphological image processing and
        %%%mas construction
        azigram_param.debug.morph=false;
        azigram_param.debug.mask_construction=false;
        
        dBspread=str2double(handles.edit_mindB.String) + [0 str2double(handles.edit_dBspread.String)];
        azigram_param.debug.dBspread=dBspread;
        
        %%%%%%How much data to read in at a time to RAM--solves clock drift
        %%%%%%problem
        params.time_chunk=10*60;
        
        output=azigram_pulse_detector(xall,Fs,Nfft,ovlap,azigram_param);
        fprintf('Finished processing, starting localization...\n');
        
        %%%Compute positions
        bearings=output.azi.wm.med;
        
        location='Eel';
        if contains(lower(location),'eel')
            Ikeep=double(DASAR_ltr)-double('X')+1;
            pos_DASAR0=[0 -14.69; 0 15.73; 0 0];  %Tako City
            rot=-20;  %Optional rotation of frame size we don't have alignment exact.
            
        elseif contains(lower(location),'tako')
            
            Ikeep=double(DASAR_ltr)-double('A')+1;
            pos_DASAR0=[4 -1; 0 0; 4+3 -1+4.4];  %Tako City
            rot=0;  %Optional rotation of frame size we don't have alignment exact.
            
        end
        pos_DASAR=([cosd(rot) sind(rot);-sind(rot) cosd(rot)]*pos_DASAR0')';
        
        output.loc=plot_pulse_positions(bearings,pos_DASAR,Ikeep);
        
        bearings=output.azi.avg;
        plot_pulse_positions(bearings,pos_DASAR,Ikeep);
        
        %%%Restrict what detections/locations to keep and write to
        %%%annotation file.
        
        filter_params.keyword='all';  %%For now pass all annotation
        %station_position.easting=pos_DASAR(:,1);station_position.northing=pos_DASAR(:,2);
        Igood=output.azi.wm.med(1,:)>=0&output.azi.wm.med(1,:)<=360;
        
        output_write=filter_structure(output,Igood);
        output_write.goodName=azigram_param.goodName;
        
        success_flag=convert_automated_PulseTracker_into_annotations(output_write,filter_params,Igood);
        
    end  %For ludovic tenorio
elseif strcmp(handles.display_view,'Spectrogram')
    
    %if strcmp(handles.display_view,'Spectrogram')
    %    axes(handles.axes1);
    %else
    %    figure;
    %end
    
    if length(x(:,1))<Nfft/2
        return
    end
    if ~(strcmp(handles.filetype,'PSD'))
        
        %%% KEY SPECTROGRAM COMMAND
        [S,FF,TT,B] = spectrogram(x(:,1),hanning(Nfft_window),round(ovlap*Nfft_window),Nfft,Fs);
        %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
        %     For real signals, variable 'B'
        %     returns the one-sided modified periodogram estimate of the PSD of each
        %     segment; for complex signals and in the case when a vector of
        %     frequencies is specified, it returns the two-sided PSD.
        handles.sgram.T		=	TT;
        handles.sgram.F		=	FF;
        % handles.sgram.B		=	B;
        handles.sgram.Nfft	=	Nfft;
        handles.sgram.ovlap	=	ovlap;
        handles.sgram.Fs	=	Fs;
        
        
        
        %%Add spectral calibration curve, if present
        if isfield(hdr,'calcurv')
            Xp_cal_fin=polyval(hdr.calcurv,FF/Fs);
            
            imagesc(TT,FF/1000,10*log10(B)+Xp_cal_fin*ones(1,length(TT)));
            %elseif isfield(hdr,'cable_factor')
            %    Xp_cal_fin=20*log10(1+hdr.cable_factor*FF);  %Unit resistance 140 ohm, capacitance 110 nF
            %    imagesc(TT,FF/1000,10*log10(B)+Xp_cal_fin*ones(1,length(TT)));
        else
            %%%KEY SPECTROGRAM IMAGE COMMAND
            imagesc(TT,FF/1000,10*log10(B));%
        end
    else
        
        ppsd=10*log10(x);
        imagesc(t,FF/1000,ppsd);
        handles.sgram.T		=	t;
        handles.sgram.F		=	FF;
        %handles.sgram.B		=	ppsd;
        handles.sgram.Nfft	=	hdr.Nfft;
        handles.sgram.ovlap	=	ovlap;
        handles.sgram.Fs	=	Fs;
    end  %if PSD
    grid on
    axis('xy')
    fmax=str2double(get(handles.edit_fmax,'String'));
    fmin=str2double(get(handles.edit_fmin,'String'));
    if fmax==0,
        ylim([0 Fs/2000]);
        set(handles.edit_fmax,'String',num2str(Fs/2000));
    else
        ylim([fmin fmax]);
    end
    %ylim([0 1]);axis('xy')
    climm(1)=str2double(get(handles.edit_mindB,'String'));
    climm(2)=climm(1)+str2double(get(handles.edit_dBspread,'String'));
    %%If switching from correlogram, reset to suggested values
    if climm(1)==0&&climm(2)<1
        climm=[40 70];
        set(handles.edit_mindB,'String',num2str(climm(1)));
        set(handles.edit_dBspread,'String',num2str(climm(2)));
        climm(2)=sum(climm);
    end
    
    caxis(climm);
    if get(handles.checkbox_grayscale,'Value')==1
        colormap(flipud(gray));
    else
        colormap(jet);
    end
    colorbar;
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
    xlim([0 str2num(handles.edit_winlen.String)]);
    if ~strcmp(handles.display_view,'Spectrogram')
        title(get(handles.text_filename,'String'));
    end
elseif strcmp(handles.display_view,'Time Series') %%Time series
    display_time_series(handles,x,Fs,hdr);
    
    
elseif strcmp(handles.display_view,'Correlogram') %%Correlogram
    display_correlogram_series(handles,x,Fs,hdr,ovlap,Nfft)
    
    
end  %Spectrogram, new fig

handles.x	=	x;
handles.Fs	=	Fs;
%tmp=ginput(2)


%Update all control buttons, depending on loaded data

%set(handles.pushbutton_GSIbearing,'vis','on');
%set(handles.pushbutton_GSI_localization,'vis','on');
%set(handles.pushbutton_next_linked_annotation,'vis','on');
%set(handles.pushbutton_previous_linked_annotation,'vis','on');
%set(handles.pushbutton_next_linked_annotation,'enable','on');
%set(handles.pushbutton_previous_linked_annotation,'enable','on');

if strcmpi(handles.filetype,'psd')
    set(handles.pushbutton_binary,'vis','off');
    set(handles.pushbutton_pausesound,'vis','off');
    set(handles.pushbutton_playsound,'vis','off');
    set(handles.pushbutton_save,'vis','off');
    set(handles.radiobutton_correlogram,'vis','off');
    set(handles.radiobutton_timeseries,'vis','off');
    
else
    set(handles.pushbutton_binary,'vis','on');
    set(handles.pushbutton_pausesound,'vis','on');
    set(handles.pushbutton_playsound,'vis','on');
    set(handles.pushbutton_save,'vis','on');
    set(handles.radiobutton_correlogram,'vis','on');
    set(handles.radiobutton_timeseries,'vis','on');
    %set(handles.radiobutton_directionality,'vis','on');
    
end

if handles.file_flags.multichannel
    status='on';
else
    status='off';
end
for II=1:length(handles.buttongroup.array)
    set(handles.buttongroup.array(II),'vis',status);
    set(handles.buttongroup.array(II),'enable',status);
end


%%%Are these vector sensor files?
[~,~,~,~,~,hdr]=load_data(handles.filetype,-1,10,1,handles);
if hdr.vector_sensor||strcmpi(handles.filetype,'gsi')||~isempty(strfind(handles.myfile,'DIFAR'))
    status='on';
else
    status='off';
end

for II=1:length(handles.buttongroup.linked)
    set(handles.buttongroup.linked(II),'vis',status);
    set(handles.buttongroup.linked(II),'enable',status);
end

for II=1:length(handles.buttongroup.GSI)
    set(handles.buttongroup.GSI(II),'vis',status);
    set(handles.buttongroup.GSI(II),'enable',status);
end

if strcmpi(lower(handles.filetype),'gsi')||handles.file_flags.multichannel
    status='on';
else
    status='off';
end

for II=1:length(handles.buttongroup.multichan)
    set(handles.buttongroup.multichan(II),'vis',status);
    set(handles.buttongroup.multichan(II),'enable',status);
end


%if strcmpi(handles.filetype,'sio')
%    set(handles.pushbutton_CSDM,'vis','on');
%end
end  %load_and_display_spectrogram
