

function varargout = AllFile_specgram_viewer(varargin)
%Test comment
% AllFile_SPECGRAM_VIEWER M-file for AllFile_specgram_viewer.fig
%      AllFile_SPECGRAM_VIEWER, by itself, creates a new AllFile_SPECGRAM_VIEWER or raises the existing
%      singleton*.
%
%      H = AllFile_SPECGRAM_VIEWER returns the handle to a new
%      AllFile_SPECGRAM_VIEWER or the handle to
%      the existing singleton*.
%
%      AllFile_SPECGRAM_VIEWER('CALLBACK',hObject,eventData,handles,...)
%      calls
%      the local
%      function named CALLBACK in AllFile_SPECGRAM_VIEWER.M with the given input arguments.
%
%      AllFile_SPECGRAM_VIEWER('Property','Value',...) creates a new AllFile_SPECGRAM_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AllFile_specgram_viewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AllFile_specgram_viewer_OpeningFcn
%      via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AllFile_specgram_viewer

% Last Modified by GUIDE v2.5 07-Mar-2013 16:48:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AllFile_specgram_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @AllFile_specgram_viewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

end

% End initialization code - DO NOT EDIT

% --- Executes just before AllFile_specgram_viewer is made visible.
function AllFile_specgram_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AllFile_specgram_viewer (see VARARGIN)


%%%%Load directory locations from setup file
% specificies three file names:
%   (1) Default directory location
%   (2) Function handles file name and location
%   (3) Multipath results file name and location

startupinfo=gui_startup_information;


%%%%Set up times based on available times in input directory...%%%%
%handles.mydir='/Users/thode/Projects/Insta-array/Sperm_Sitka/May9_Kendall_Cobra/Lucy';
handles.mydir=startupinfo.default_directory;
handles.inputdir=startupinfo.default_inputfiledir;
handles.annotation_file=startupinfo.annotation_file;
try
    handles.calibration_DASAR2007_dir=startup_info.calibration_DASAR2007_dir;
catch
    disp('2007 calibration directory not defined');
end
set(handles.text_filename,'String',[handles.mydir ]);

mydir=pwd;
try
    cd(handles.mydir);
end
handles.filetype='MT';
[handles,errorflag]=set_slider_controls(handles,handles.filetype);
%Make GSIbearing button and CSDM invisible
set(handles.pushbutton_GSIbearing,'Vis','off');
set(handles.pushbutton_GSI_localization,'Vis','off');
set(handles.pushbutton_CSDM,'Vis','off');
set(handles.pushbutton_Mode,'Vis','off');
set(handles.pushbutton_tilt,'Vis','off');
set(handles.pushbutton_modalfiltering,'Vis','off');

set(handles.edit_maxfreq,'String','0');
set(handles.edit_minfreq,'String','0');
cd(mydir);
%%%Create list of function handles to be used in other specialized
%%%applications



% Choose default command line output for AllFile_specgram_viewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using AllFile_specgram_viewer.
if strcmp(get(hObject,'Visible'),'off')
    %A = imread('/Users/thode/Personal/meriah_n_me.jpg','jpg');
    %image(A);set(gca,'vis','off');
    %plot(rand(5));
end
end

% UIWAIT makes AllFile_specgram_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AllFile_specgram_viewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


end
% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


load_and_display_spectrogram(hObject,eventdata,handles);



end

function load_and_display_spectrogram(hObject,eventdata,handles)

ha	=	handles.axes1;
cla(ha);
axes(ha);

tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
mydir=pwd;
Ichan=str2num(get(handles.edit_chan,'String'));  %Hardwire first channel

[x,t,Fs,tstart,junk,hdr]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,Ichan,handles);
if max(t)<tlen
    set(handles.edit_windowlength,'String',num2str(max(t)));
end

if size(x,2)>1
    x=x';
end

Fs=round(Fs);

mymaxfreq=str2double(get(handles.edit_maxfreq,'String'));

if mymaxfreq==0||mymaxfreq>0.8*Fs/2
    set(handles.edit_maxfreq,'String',num2str(0.8*Fs/2));
end
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if strcmp(handles.display_view,'Spectrogram')||strcmp(handles.display_view,'New Fig')
    %[B,FF,TT]=specgram(x(:,1),Nfft,Fs,hanning(Nfft),round(ovlap*Nfft));
    [S,FF,TT,B] = spectrogram(x(:,1),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
    %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
    
	%	Jit: save spectrogram for ship detection metrics
	handles.Specgram.PSD	=	B;
	handles.Specgram.T		=	TT;
	handles.Specgram.F		=	FF;
	handles.Specgram.Fs		=	Fs;
	handles.Specgram.Nfft	=	Nfft;
	guidata(hObject,handles);
	
	
    if strcmp(handles.display_view,'Spectrogram')
        axes(handles.axes1);
    else
        figure;
    end
    
    %%Add spectral calibration curve, if present
    if isfield(hdr,'calcurv')
        Xp_cal_fin=polyval(hdr.calcurv,FF/Fs);
        
        imagesc(TT,FF/1000,10*log10(B)+Xp_cal_fin*ones(1,length(TT)));
    %elseif isfield(hdr,'cable_factor')
   %    Xp_cal_fin=20*log10(1+hdr.cable_factor*FF);  %Unit resistance 140 ohm, capacitance 110 nF
    %    imagesc(TT,FF/1000,10*log10(B)+Xp_cal_fin*ones(1,length(TT)));
    else
        
        imagesc(TT,FF/1000,10*log10(B));%
    end
    grid on
    axis('xy')
    fmax=str2num(get(handles.edit_fmax,'String'));
    fmin=str2num(get(handles.edit_fmin,'String'));
    if fmax==0,
        ylim([0 Fs/2000]);
        set(handles.edit_fmax,'String',num2str(Fs/2000));
    else
        ylim([fmin fmax]);
    end
    %ylim([0 1]);axis('xy')
    climm(1)=str2num(get(handles.edit_mindB,'String'));
    climm(2)=climm(1)+str2num(get(handles.edit_dBspread,'String'));
    %%If switching from correlogram, reset to suggested values
    if climm(1)==0&climm(2)<1
        climm=[40 70];
        set(handles.edit_mindB,'String',num2str(climm(1)));
        set(handles.edit_dBspread,'String',num2str(climm(2)));
        climm(2)=sum(climm);
    end
    
    caxis(climm);
    if get(handles.checkbox_grayscale,'Value')==1,
        colormap(flipud(gray));
    else
        colormap(jet);
    end
    colorbar;
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
    handles.T=TT;
    handles.F=FF;
    if ~strcmp(handles.display_view,'Spectrogram')
        title(get(handles.text_filename,'String'));
    end
elseif strcmp(handles.display_view,'Time Series') %%Time series
    
    %%Check that we are looking at acoustic data
    if isempty(findstr(handles.myfile,'Press'))
        disp('Enter min and max frequency:');
        tmp=ginput(2);
        if ~isempty(tmp)&&size(tmp,1)==2
            freq=sort(tmp(:,2))*1000;
            minfreq=freq(1);maxfreq=freq(2);
            %y=quick_filter(x(:,1),Fs,freq(1),freq(2))
            frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
            [N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.1],Fs);
            B = firpm(N,Fo,Ao,W);
            
            y=filter(B,1,x(:,1)-mean(x(:,1)));
        else
            y=x(:,1)-mean(x(:,1));
        end
    else
        y=x(:,1)/1000;
    end
    
    t=(1:length(x(:,1)))/Fs;
    plot(t,y);grid on;
    xlabel('Time (sec)');ylabel('Amplitude');
elseif strcmp(handles.display_view,'Correlogram') %%Correlogram
    fmax=1000*str2num(get(handles.edit_fmax,'String'));
    fmin=1000*str2num(get(handles.edit_fmin,'String'));
    param.ovlap=ovlap;
    param.Nfft=Nfft;
    prompt1={'ICI range (s)','Correlation sample time (s)','Teager-Kaiser treatment?','Incoherent=1, Coherent=0 ?'};
    dlgTitle1='Parameters for correlogram...';
    def1={'[0.01 .25]', '0.25','0','1'};
    answer=inputdlg(prompt1,dlgTitle1,1,def1);
    param.ici_range=eval(answer{1});
    param.time_sample=str2num(answer{2});
    param.teager=str2num(answer{3});
    alg_chc=str2num(answer{4});
    
    if alg_chc==1
        [S,FF,TT,B] = spectrogram(x(:,1),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
        [mean_corr_org,tindex,TT_plot,pwr,pwr_tot,yscale]= create_incoherent_correlogram(TT,FF,B,param,fmin,fmax);
        
        
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
    climm(1)=str2num(get(handles.edit_mindB,'String'));
    climm(2)=climm(1)+str2num(get(handles.edit_dBspread,'String'));
    
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

handles.x=x;
handles.Fs=Fs;
%tmp=ginput(2)

if strcmp(lower(handles.filetype),'gsi')
    set(handles.pushbutton_GSIbearing,'vis','on');
    set(handles.pushbutton_GSI_localization,'vis','on');
else
    set(handles.pushbutton_GSIbearing,'vis','off');
    set(handles.pushbutton_GSI_localization,'vis','off');
end

if strcmp(lower(handles.filetype),'mdat')||strcmp(lower(handles.filetype),'wav')||strcmp(lower(handles.filetype),'simulation')
    set(handles.pushbutton_CSDM,'vis','on');
    set(handles.pushbutton_Mode,'vis','on');
    set(handles.pushbutton_tilt,'vis','on');
    set(handles.pushbutton_modalfiltering,'vis','on');
else
    set(handles.pushbutton_CSDM,'vis','off');
    set(handles.pushbutton_Mode,'vis','off');
    set(handles.pushbutton_modalfiltering,'vis','off');
end

guidata(hObject, handles);
end

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end
% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if isempty(handles.mydir),
%    handles.mydir = uigetdir('/Users/thode/Projects/Insta-array');
%else
%    handles.mydir = uigetdir(handles.mydir);
%end


mydir=pwd;
handles.filetype=get(get(handles.uipanel_type,'SelectedObject'),'String');
if ~strcmp(handles.filetype,'Simulation')
    menustr={['*.' lower(handles.filetype)];['*.' (handles.filetype)];'*.*';['*1004*.' lower(handles.filetype)]};
else
    menustr={'*.mat'};
end
if ~isempty(handles.mydir),
    disp(handles.mydir);
    try
        cd(handles.mydir);
    end
end
[filename, pathname, filterindex] = uigetfile(menustr, 'Select a file:');

handles.mydir = pathname;
handles.myfile=filename;
%yes_wav=get(handles.togglebutton_getwav,'value');
[handles,errorflag]=set_slider_controls(handles,handles.filetype);
set(handles.text_filename,'String',[handles.mydir '/' handles.myfile]);
cd(mydir);

guidata(hObject, handles);
end
% Update handles structure

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
end
% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

end

% --- Executes during object creation, after setting all properties.
function edit_datestr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_datestr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

end

function edit_datestr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_datestr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myval=get(handles.slider_datestr,'Value');

%tmin=datenum(get(handles.text_mintime,'String'));
%tmax=datenum(get(handles.text_maxtime,'string'));
tmin=handles.min_time;
tmax=handles.max_time;

olddate=tmin+myval*(tmax-tmin);

try
    handles.tdate_start=datenum(get(hObject,'String'));
    set(handles.text_datestr_demo,'String',get(hObject,'String'));
    newval=(handles.tdate_start-tmin)/(tmax-tmin);
    set(handles.slider_datestr,'Value',newval);
catch
    errordlg('Incorrect datestr');
    set(hObject,'String',get(handles.text_datestr_demo,'String'));
    set(handles.slider_datestr,'Value',myval);
end
guidata(hObject, handles);

end

% Hints: get(hObject,'String') returns contents of edit_datestr as text
%        str2double(get(hObject,'String')) returns contents of edit_datestr as a double


% --- Executes during object creation, after setting all properties.
function edit_windowlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

end

function edit_windowlength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowlength as text
%        str2double(get(hObject,'String')) returns contents of edit_windowlength as a double
tlen=str2num(get(hObject,'String'));
if tlen<=240
    maxx=get(handles.slider_datestr,'max');
    minn=get(handles.slider_datestr,'min');
    slider_step=get(handles.slider_datestr,'sliderstep');
    slider_step(1)=datenum(0,0,0,0,0,tlen)/(maxx-minn);
    set(handles.slider_datestr,'sliderstep',slider_step)
else
    errordlg('window length greater than 240 sec');
end



guidata(hObject, handles);

end
% --- Executes during object creation, after setting all properties.
function popupmenu_Nfft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
contents=get(hObject,'String');
set(hObject,'Value',4);

end
% --- Executes on selection change in popupmenu_Nfft.
function popupmenu_Nfft_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_Nfft contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Nfft

end
% --- Executes during object creation, after setting all properties.
function popupmenu_ovlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ovlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
contents=get(hObject,'String');
set(hObject,'Value',4);
end

% --- Executes on selection change in popupmenu_ovlap.
function popupmenu_ovlap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ovlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_ovlap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ovlap

end
% --- Executes during object creation, after setting all properties.
function slider_datestr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_datestr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

end
% --- Executes on slider movement.
function slider_datestr_Callback(hObject, eventdata, handles)
% hObject    handle to slider_datestr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
%keyboard;
myval=get(hObject,'Value');

%%Following commented out because tmax fails to read in WAV data correctly.
%tmin=datenum(get(handles.text_mintime,'String'));
%tmax=datenum(get(handles.text_maxtime,'string'));

tmin=handles.min_time;
tmax=handles.max_time;

dT_min=(tmax-tmin)*24*60;
set(hObject,'sliderstep',[(20/60)/dT_min 1/dT_min]); % 20sec and 1min increments

datenumm=tmin+myval*(tmax-tmin);
set(handles.edit_datestr,'String',datestr(datenumm,0));
set(handles.text_datestr_demo,'String',datestr(datenumm,0));


end
% --- Executes during object creation, after setting all properties.
function edit_mindB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mindB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

end

function edit_mindB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mindB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mindB as text
%        str2double(get(hObject,'String')) returns contents of edit_mindB as a double
climm=get(handles.axes1,'clim');
climm=str2num(get(hObject,'String'))+[0 diff(climm)];
caxis(climm);colorbar

end
% --- Executes during object creation, after setting all properties.
function edit_dBspread_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dBspread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

end

function edit_dBspread_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dBspread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dBspread as text
%        str2double(get(hObject,'String')) returns contents of edit_dBspread as a double
climm=get(handles.axes1,'clim');
climm(2)=str2num(get(hObject,'String'))+climm(1);
caxis(climm);colorbar

end
% --- Executes during object creation, after setting all properties.
function edit_minfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
handles.bflag=0;
guidata(hObject, handles);


end
function edit_minfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_minfreq as a double

handles.bflag=0;  %Has a filter been created?
guidata(hObject, handles);

end


% --- Executes during object creation, after setting all properties.
function edit_maxfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

end




function edit_maxfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_maxfreq as a double

handles.bflag=0;  %Has a filter been created?
guidata(hObject, handles);

end

% --- Executes on button press in pushbutton_playufsound.
function pushbutton_playufsound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playufsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.Fs>48000
    n=ceil(handles.Fs/48000);
    disp('Resampling sound');
    
    soundsc(downsample(handles.x,n),handles.Fs/n);
else
    soundsc(handles.x,handles.Fs);
end
end

% --- Executes on button press in pushbutton_recomputeFilter.
function pushbutton_recomputeFilter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_recomputeFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keep_trying=1;
while keep_trying
    maxfreq=str2num(get(handles.edit_maxfreq,'String'));
    minfreq=str2num(get(handles.edit_minfreq,'String'));
    
    try
        frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
        [n,fo,mo,w] = remezord(frange, [0 1 0], [0.1 0.01 0.1], handles.Fs );
        handles.b = remez(n,fo,mo,w);
        handles.bflag=1;
        guidata(hObject, handles);
        keep_trying=0;
    catch
        prompt={'Min (Hz)','Max (Hz)'};
        name='Your selections for frequency range are not valid';
        numlines=1;
        
        if isfield(handles,'Fs')
            defaultanswer={'20',num2str(0.8*handles.Fs/2)};
        else
            defaultanswer={'20',num2str(0.8*1000/2)};
        
        end
        
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        set(handles.edit_minfreq,'String',answer{1});
        set(handles.edit_maxfreq,'String',answer{2});
    end
end
end


% --- Executes on button press in pushbutton_playsound.
function pushbutton_playsound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.bflag==0
    keep_trying=1;
    while keep_trying
        maxfreq=str2num(get(handles.edit_maxfreq,'String'));
        minfreq=str2num(get(handles.edit_minfreq,'String'));
        
        try
            frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
            [n,fo,mo,w] = remezord(frange, [0 1 0], [0.1 0.01 0.1], handles.Fs );
            handles.b = remez(n,fo,mo,w);
            handles.bflag=1;
            guidata(hObject, handles);
            keep_trying=0;
        catch
            prompt={'Min (Hz)','Max (Hz)'};
            name='Your selections for frequency range are not valid';
            numlines=1;
            if isfield(handles,'Fs')
                defaultanswer={'20',num2str(0.8*handles.Fs/2)};
            else
                defaultanswer={'20',num2str(0.8*1000/2)};
                
            end
            answer=inputdlg(prompt,name,numlines,defaultanswer);
            set(handles.edit_minfreq,'String',answer{1});
            set(handles.edit_maxfreq,'String',answer{2});
        end
    end
    
end
if handles.Fs>48000
    disp('Resampling sound');
    n=ceil(handles.Fs/48000);
    x=filter(handles.b,1,handles.x);
    soundsc(downsample(x,n),handles.Fs/(1*n));
else
    
    soundsc(filter(handles.b,1,handles.x),handles.Fs/1);
end

end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;
Ichan='all';
[x,t,Fs,tstart,junk,hdr]=load_data(handles.filetype,handles.min_time,tdate_start,tlen,Ichan,handles);

if ~isempty(findstr(lower(computer),'mac'))
    Islash=findstr(handles.mydir,'/');
else
    Islash=findstr(handles.mydir,'\');
end
chan=get(handles.edit_chan,'String');

save_name=sprintf('soundsamp%s_%s',handles.mydir((Islash(end)+1):end),datestr(tdate_start,30));
disp(sprintf('Saving %s ...',save_name));

try
    if size(x,2)>size(x,1)
        x=x';
    end
    for Ichan=1:size(x,2)
        xfilt(:,Ichan)=filter(handles.b,1,x(:,Ichan));
    end
    wavwrite(xfilt/(1.1*max(max(abs(xfilt)))),Fs,save_name);

catch
    disp('No filtering desired... exists; saving raw acoustic data to WAV');
    xfilt=[];
    wavwrite(x/(1.1*max(max(abs(x)))),Fs,save_name);

end
save(save_name,'x','xfilt','Fs','tdate_start','hdr');

end

function [handles,errorflag]=set_slider_controls(handles,filetype)
%%Set min, max and other times associated with slider controls
mydir=pwd;
errorflag=0;
handles.min_time=-1;
handles.max_time=-1;

if isempty(handles.mydir)
    handles.mydir=pwd;
end
%cd(handles.mydir);

try
    [x,t,Fs,minn,maxx]=load_data(filetype,-1,-1,1,1,handles);
catch %no file selected
    %errordlg(sprintf('No %s file selected',filetype));
    errorflag=1;
    cd(mydir);
    return
end


set(handles.edit_datestr,'String',datestr(minn,0));
%set(handles.slider_datestr,'Max',maxx);
set(handles.slider_datestr,'Max',1);

%set(handles.slider_datestr,'Min',minn);
set(handles.slider_datestr,'Min',0);

set(handles.text_mintime,'String',datestr(minn,0));

set(handles.slider_datestr,'Value',0.5);
set(handles.edit_datestr,'String',datestr(0.5*(minn+maxx),0));

set(handles.text_maxtime,'String',datestr(maxx,0));

handles.min_time=minn;
handles.max_time=maxx;

set(handles.edit_fmax,'String',Fs/2000);
set(handles.edit_fmin,'String',0);
%slider_step(1) = datenum(0,0,0,0,0,str2num(get(handles.edit_windowlength,'String')))/(maxx-minn);
%slider_step(2) = min([datenum(0,0,0,0,5,0) 0.1*(maxx-minn)])/(maxx-minn);
%set(handles.slider_datestr,'sliderstep',slider_step)
% keyboard
cd(mydir);

end

% --- Executes on button press in pushbutton_print.
function pushbutton_print_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.filetype,'MDAT')
    mychc=menu('Which print option?','Single phone','all phones');
    
    if mychc==2  %%Multiphone data...
        
        figure;
        tdate_start=datenum(get(handles.edit_datestr,'String'));
        tlen=str2num(get(handles.edit_windowlength,'String'));
        contents=get(handles.popupmenu_Nfft,'String');
        Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});
        
        contents=get(handles.popupmenu_ovlap,'String');
        ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;
        
        mydir=pwd;
        figure(1)
        set(gcf,'pos',[291         628        1513         991]);
        Iplot=0;
        %[x,t,Fs,tstart,junk1,head]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,'all',handles);
        [x,t,Fs,tstart,junk1,head]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,'all',handles);
        
        for Ichan=1:head.Nchan
            
            if Iplot<8
                Iplot=Iplot+1;
            else
                figure(2)
                set(gcf,'pos',[291         628        1513         991]);
                Iplot=1;
                
            end
            subplot(4,2,Iplot);
            
            [S,FF,TT,B] = spectrogram(x(Ichan,:),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
            %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
            %axes(handles.axes1);
            imagesc(TT,FF,10*log10(B));%
            axis('xy')
            fmax=str2num(get(handles.edit_fmax,'String'));
            fmin=str2num(get(handles.edit_fmin,'String'));
            if fmax==0,
                ylim([0 Fs/2]);
                set(handles.edit_fmax,'String',num2str(Fs/2));
            else
                ylim([fmin fmax]*1000);
            end
            %ylim([0 1]);axis('xy')
            climm(1)=str2num(get(handles.edit_mindB,'String'));
            climm(2)=climm(1)+str2num(get(handles.edit_dBspread,'String'));
            caxis(climm);
            if get(handles.checkbox_grayscale,'Value')==1,
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
                text(0.1,0.9,num2str(head.rd(Ichan),3),'color','y','units','norm','fontsize',14,'fontweight','bold');
            catch
                disp('No elements depths provided');
            end
            hh=colorbar('East');
            set(hh,'fontweight','bold','fontsize',14,'ycolor','y')
            
            %             if Ichan==1
            %                 xall=zeros(8,length(x));
            %             end
            %xall(Ichan,:)=x;
        end
        if ~isempty(findstr(lower(computer),'mac'))
            Islash=findstr(handles.mydir,'/')+1;
        else
            Islash=findstr(handles.mydir,'\')+1;
        end
        save_name=sprintf('soundNfft%i_%s.%s_%s_%s',Nfft,handles.mydir((Islash(end-1)):(end-1)),handles.myfile,datestr(tdate_start,30),'all_channels');
        pause
        
        for I=1:2
            if any(get(0,'child')==I)
                save_name1=sprintf('%s_%i',save_name,I);
                disp(sprintf('Printing %s ...',save_name1));
                figure(I)
                orient landscape
                print(I,'-djpeg',[save_name1 '.jpg']);
                save([save_name1 '.mat'],'x','Fs');
                close(I);
            end
            
        end
        
        return
    end
end

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if ~strcmp(handles.display_view,'New Fig')
    figchc=[];
    axes(handles.axes1);
else
    figure(gcf);
    chcc=get(0,'child');
    Igoodd=find(chcc-round(chcc)==0);
    chcc=chcc(Igoodd);
    for III=1:length(chcc)
        tmp{III}=chcc(III);
    end
    figchc=menu('Select a figure number:',tmp);
    figure(chcc(figchc));
end
tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
chan=get(handles.edit_chan,'String');
Idot=findstr(handles.myfile,'.')-1;
save_name=sprintf('soundsamp_%s_%s_%s',handles.myfile(1:Idot),datestr(tdate_start,30),chan);
disp(sprintf('Printing %s ...',save_name));

orient landscape
print(figchc,'-djpeg',save_name);
print(figchc,'-dtiff',save_name);

%print('-depsc',save_name);

end
% --- Executes on button press in pushbutton_annotate.
function pushbutton_annotate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_annotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Click on two points defining a rectangle: ');
tmp=ginput(2);

start_time=datenum(get(handles.edit_datestr,'String'))+datenum(0,0,0,0,0,min(tmp(:,1)));
min_freq=num2str(1000*min(tmp(:,2)));
max_freq=num2str(1000*max(tmp(:,2)));
duration=num2str(abs(tmp(2,1)-tmp(1,1)));

call_types={'Pulsive','FM modulated'};
Icall=menu('Signal type?',call_types);
call_type=call_types{Icall};

prompt={'File name','pulse or FM?','Call Type','Min Freq(Hz)','Max Freq(Hz)','Duration(sec)','Number_pulses','Number_harmonics','modulation (Hz)'};
if Icall==1  %Pulsive
    def={handles.annotation_file,'pulse','S1',min_freq,max_freq,duration,'10','-1','0'};
else
    def={handles.annotation_file,'FM','moan',min_freq,max_freq,duration,'-1','1','0'};
end
dlgTitle=sprintf('Annotation for event at %s',datestr(start_time));
lineNo=1;
answer1=inputdlg(prompt,dlgTitle,lineNo,def);


fid=fopen(answer1{1},'a');
if ftell(fid)==0
    fprintf(fid,'%s','Start Date and Time,');
    for I=2:length(prompt)
        fprintf(fid,'%s,',prompt{I});
    end
    fprintf(fid,'\n');
end

fprintf(fid,'%s,',datestr(start_time,0));
for I=2:length(prompt)
    fprintf(fid,'%s, ',answer1{I});
end
fprintf(fid,'\n');
fclose(fid);

end

% --- Executes on button press in pushbutton_selectpoints.
function pushbutton_selectpoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt={'Number of selections'};
def={'3'};
start_time=datenum(get(handles.edit_datestr,'String'));
dlgTitle=sprintf('Selecting points at %s',datestr(start_time));
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
tmp=ginput(str2num(answer{1}));
disp('Absolute times:')
disp(datestr(start_time+datenum(0,0,0,0,0,tmp(:,1))));
disp('Relative times and frequencies:');
disp(tmp(:,1)');
disp(tmp(:,2)');
disp('Time elaspsed from first point:');
disp(tmp(:,1)'-tmp(1,1));
end

% --- Executes on button press in pushbutton_fileread.
function pushbutton_fileread_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fileread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Assumes datenumber is first part of line
tline=fgetl(handles.fid);
%Remove leading blanks;
tline=fliplr(deblank(fliplr(tline)));
disp(tline);
Idash=findstr('-',tline);
Icolon=findstr(':',tline);
newtime=tline(1:(max(Icolon)+2));
if tline~=-1,
    try
        handles.tdate_start=datenum(newtime);
        set(handles.edit_datestr,'String',newtime);
        set(handles.text_datestr_demo,'String',newtime);
        set(handles.slider_datestr,'Value',handles.tdate_start);
    catch
        errordlg('Incorrect datestr');
        set(hObject,'String',get(handles.text_datestr_demo,'String'));
    end
    
    guidata(hObject, handles);
else
    disp('End of file reached');
    fclose(handles.fid);
end

end
% --------------------------------------------------------------------
function OpenMenuItem_inputfile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_inputfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.mydir = uigetdir(handles.inputdir);
mydir=pwd;
cd(handles.inputdir);
handles.inputfilename=uigetfile('*','Select input file:');
handles.fid=fopen(handles.inputfilename);
cd(mydir);
guidata(hObject, handles);

end
% --- Executes on button press in checkbox_restart.
function checkbox_restart_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_restart
if get(hObject,'Value')==1,
    try,
        fclose(handles.fid);
        mydir=pwd;
        cd(handles.inputdir);
        handles.fid=fopen(handles.inputfilename);
        cd(mydir);
        guidata(hObject, handles);
    end
end
end

% --- Executes on button press in checkbox_grayscale.
function checkbox_grayscale_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_grayscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_grayscale

if get(hObject,'Value')==1,
    colormap(flipud(gray));
else
    colormap(jet);
end
end

% --- Executes during object creation, after setting all properties.
function edit_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'String',num2str(0));

end
function edit_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fmax as text
%        str2double(get(hObject,'String')) returns contents of edit_fmax as a double
fmin=str2num(get(handles.edit_fmin,'String'));
fmax=str2num(get(hObject,'String'));
ylim([fmin fmax]);
end

function edit_fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'String',num2str(0));

end
function edit_fmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fmax as text
%        str2double(get(hObject,'String')) returns contents of edit_fmax as a double

fmin=str2num(get(hObject,'String'));
fmax=str2num(get(handles.edit_fmax,'String'));
ylim([fmin fmax]);

%%%%%%%%%gui_startup_information.m%%%%%%%%%%%%%%%%
% Aaron Thode
% November 22, 2004
% Contains directory locations and names of various files
%   used by the Schultz GUI system

end
function [startup_info]=gui_startup_information

%%%%Directory where GSI_specgram_viewer.m and function_handles files will be
%%%%stored....
[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);
startup_info.calibration_DASAR2007_dir=[];

if strcmp(local_machine,'macmussel.ucsd.edu')
    
    startup_info.base_directory='/Users/Shared/MATLAB/AllFile_specgram_viewer';
    
    
    startup_info.default_directory='/Volumes/';
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    
    %%Default annotation file
    startup_info.annotation_file='annotated.txt';
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
    startup_info.calibration_DASAR2007_dir='';
elseif findstr(local_machine,'Jan-Straleys')
    
    startup_info.base_directory='/Users/janstraley/SEASWAP_2011';
    
    
    startup_info.default_directory='/Users/janstraley/SEASWAP_2011';
    startup_info.default_inputfiledir='/Users/janstraley/SEASWAP_2011';
    
    %%Default annotation file
    startup_info.annotation_file='annotated.txt';
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
    startup_info.calibration_DASAR2007_dir='';
elseif findstr(local_machine,'Janice')
    
    startup_info.base_directory='~/Desktop/MATLAB/AllFile_specgram_viewer';
    
    
    startup_info.default_directory='~/Desktop';
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    
    %%Default annotation file
    startup_info.annotation_file='annotated.txt';
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
    startup_info.calibration_DASAR2007_dir='';
    
elseif strcmp(local_machine,'KatsMacPro.local')  %KatMacGSI
    startup_info.base_directory='/Users/thode/MATLAB/AllFile_specgram_viewer';
    
    
    startup_info.default_directory='/Volumes';
    startup_info.default_inputfiledir='/Volumes';
    
    startup_info.annotation_file='annotated.txt';
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
    
    
elseif ~isempty(findstr('dgrebner',local_machine))
    
    startup_info.base_directory='/Users/dgrebner/Desktop/';
    startup_info.default_directory=startup_info.base_directory;
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    startup_info.annotation_file='annotated.txt';
    startup_info.calibration_DASAR2007_dir='';
    
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
elseif ~isempty(findstr('thode',local_machine))
    
    startup_info.base_directory='/Users/thode/Projects/Arctic_2010/Data/Bottom_Unit';
    startup_info.default_directory=startup_info.base_directory;
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    startup_info.annotation_file='annotated.txt';
    
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    startup_info.calibration_DASAR2007_dir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/RawData';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
else
    startup_info.base_directory='~/Desktop/MATLAB';
    startup_info.default_directory='~/Desktop';
    startup_info.default_inputfiledir='~/Desktop';
    startup_info.annotation_file='annotated.txt';
    
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
    
end

end

% --- Executes on button press in pushbutton_binary.
function pushbutton_binary_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
Ichan=str2num(get(handles.edit_chan,'String'));
[x,t,Fs,tstart]=load_data(handles.filetype,handles.min_time,tdate_start,tlen,Ichan,handles);

Nmax=(2^16)-1;
Fs_want=125000;  %Actual playback rate

if Fs<Fs_want
    [q,p]=rat(Fs/Fs_want,0.0001);
    x=resample(x,p,q);
end

y_scale=(x-min(x))*(Nmax/(max(x)-min(x)));
Istart=input('Enter an integer for filename: ');

fid=fopen(['PLAY.' int2str(Istart)],'w','ieee-be');
COUNT = fwrite(fid,y_scale,'uint16');
fclose(fid)

end

function tstart=convert_date(data,delimiter)
%%%%%%convert_date.m%%%%%
% Extract a datenumber from a filename string...
% function tstart=convert_date(data,delimiter),
%%get start date from file name.
%%Two types of filenames, one "MT" style the other MATLAB 'T" style
%%Start intelligent search for a timestamp in string name...


mt_style_flag=~isempty(findstr('Sound',data));
fname_bounds=findstr(data,delimiter);

if length(fname_bounds)==1,
    fname_bounds(2)=length(data)+1;
end
Idot=findstr(data,'.');
fname_bounds=sort([fname_bounds Idot(1)]);

nm=[];
if mt_style_flag==1,
    Isound=findstr(data,'Sound')+6;
    yr=str2num(data(Isound:(Isound+3)));
    mo=str2num(data((Isound+5):(Isound+6)));
    dy=str2num(data((Isound+7):(Isound+8)));
    hr=str2num(data((Isound+9):(Isound+10)));
    mn=str2num(data((Isound+11):(Isound+12)));
    %sec=str2num(data((Isound+13):(Isound+14)));
    tstart=datenum(yr,mo,dy,hr,mn,0);
    
else
    try
        
        for I=1:(length(fname_bounds)-1),
            test=data((fname_bounds(I)+1):(fname_bounds(I+1)-1));
            if length(test)==15,
                if strcmp(test(9),'T'),
                    nm=test;
                end
            end
        end
        if ~isempty(nm)
            yr=str2num(nm(1:4));
            mo=str2num(nm(5:6));
            dy=str2num(nm(7:8));
            hr=str2num(nm(10:11));
            mn=str2num(nm(12:13));
            sec=str2num(nm(14:15));
            
            tstart=datenum(yr,mo,dy,hr,mn,sec);
        end
    catch
        error('String does not contain valid datestr');
    end
end
end

% --- Executes during object creation, after setting all properties.
function pushbutton_update_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
% --- Executes during object creation, after setting all properties.
function pushbutton_annotate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_annotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
% --- Executes during object creation, after setting all properties.
function pushbutton_selectpoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
% --- Executes during object creation, after setting all properties.
function pushbutton_playsound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_playsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
function uipanel_type_SelectionChangeFcn(hObject,eventdata,handles)


switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_MT'
        disp('MT selected')
        
    case 'radiobutton_WAV'
        disp('WAV selected');
        
    case 'radiobutton_GSI'
        %  set(hObject,'filetype','GSI');
        disp('GSI selected');
    case 'radiobutton_DAT'
        
        % Continue with more cases as necessary.
    case 'radiobutton_MDAT'
        disp('MDAT selected');
    case 'radiobutton_ADI'
        disp('ADI selected');
        % Continue with more cases as necessary.
    otherwise
        % Code for when there is no match.
end


end
% --- Executes during object creation, after setting all properties.
function uipanel_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'SelectionChangeFcn',@uipanel_type_SelectionChangeFcn);


end

function [x,t,Fs,tmin,tmax,head]=load_data(filetype,tstart_min,tdate_start,tlen,Ichan,handles)
persistent fs keyword
mydir=handles.mydir;
myfile=handles.myfile;
teager=get(handles.checkbox_teager,'Value');

x=[];
t=[];
tmin=[];
tmax=[];
head=[];
switch filetype
    case 'Simulation'
        
        simulated=load([mydir '/' myfile]);
        Fs=simulated.fs;
        
        x=simulated.x_sweep';
        if strcmp(Ichan,'all')
            Ichan=1:size(x,2);
        end
        if max(Ichan)>max(size(x,2)),
            disp(['Channel too high, restricting to ' int2str(max(size(x,2)))]);
            Ichan=max(size(x,2));
        end
        
        x=x(:,Ichan);
        
        
        if tdate_start==-1
            tmin=datenum(0);
        elseif tdate_start>0 %We are trimming beginning time
            tmin=tdate_start;
            tmp=datevec(tdate_start);
            dn=1+floor(Fs*tmp(end));
            x=x(dn:end,:);
            t=t(dn:end);
        else
            tmin=tdate_start;
        end
        
        if tlen*Fs<size(x,1)
            x=x(1:floor(tlen*Fs),:);
        end
        t=(1:length(x))/Fs;
        
        tmax=tmin+datenum(0,0,0,0,0,max(t));
        head.rd=simulated.rd;
        head.Nchan=size(x,2);
        x=x';
    case 'GSI'
        if strcmp(Ichan,'all')
            Ichan=1:3;
        end
        if max(Ichan)>3,
            disp('Channel too high, restricting to 3');
            Ichan=3;
        end
        [x,t,head]=readGSIfile([mydir '/' myfile],tdate_start,tlen,Ichan,'datenum','nocalibrate');
        if isempty(keyword)
            keyword=input('Enter a keyword for GSI calibration [DASARC]:','s');
            if isempty(keyword)
                keyword='DASARC';
            end
        end
        x=calibrate_GSI_signal(x, keyword);
        
        Fs=head.Fs;
        
        tmin=datenum(1970,1,1,0,0,head.ctbc);
        tmax=tmin+datenum(0,0,1,0,0,0);
        head.Nchan=length(Ichan);
    case 'MT'
        %[x,t,Fs]=load_mt_mult(handles.mydir,tdate_start,tlen);
        head=read_mt_header([mydir '/' myfile]);
        
        tmin=head.tstart;
        tmax=head.tend;
        Fs=head.Fs;
        if Ichan>1
            disp('WARNING: load_data: MT can only have one channel');
            Ichan=1;
        end
        if tdate_start>0
            tdate_vec=datevec(tdate_start-tmin);
            nsec=tdate_vec(6)+60*tdate_vec(5)+3600*tdate_vec(4);
            [x,t]=load_mt([mydir '/' myfile],nsec,tlen);
        end
        head.Nchan=1;
    case 'DAT'
        [x,tmin,tmax,fs]=read_dat_file([mydir '/' myfile],[],-1,tlen,0);
        
        if isempty(fs)
            fs=input('Enter sampling rate in Hz:');
        end
        Fs=fs;
        [x,tmin,tmax]=read_dat_file([mydir '/' myfile],Fs,tdate_start,tlen,0); %output in uPa
        %[x,tmin,tmax]=read_dat_file([mydir '/' myfile],Fs,tdate_start,tlen,1);  %Voltage output
        
        head.Nchan=1;
        switch fs
            case 50000
                head.calcurv=[
                    1.179288464673746e+06
                    -3.417289147406752e+06
                    3.972100408634462e+06
                    -2.459193259685826e+06
                    8.904700994689314e+05
                    -1.924134277822444e+05
                    2.476608484423531e+04
                    -2.235739303825218e+03
                    2.904887584919255e+02
                    -5.381149759460806e+00
                    7.841554559708414e-03
                    ];
            case 6250
                head.calcurv=[
                    -9.864342626384007e+06
                    2.675183405132254e+07
                    -3.072255757018830e+07
                    1.946983114345214e+07
                    -7.445224085881455e+06
                    1.766054734429601e+06
                    -2.570588847834060e+05
                    2.188411119767746e+04
                    -9.803725367146685e+02
                    1.959124505642275e+01
                    -2.811936435415921e-01];
            case 12500
                head.calcurv=[
                    9.262441626302190e+07
                    -2.151487191990283e+08
                    2.069375942078056e+08
                    -1.063702102525421e+08
                    3.153159716612202e+07
                    -5.458152352141772e+06
                    5.382152297627985e+05
                    -2.765563363629215e+04
                    6.113088605208859e+02
                    -1.301582987521525e+00
                    -1.634557871607174e-01];
        end
        
    case 'ADI'
        [x,tmin,tmax,fs]=read_adi_file(mydir,myfile,[],0,tlen,0);
        
        if isempty(fs)
            fs=input('Enter sampling rate in Hz:');
        end
        Fs=fs;
        [x,tmin,tmax]=read_adi_file(mydir,myfile,Fs,tdate_start,tlen,0); %Output in uPa
        %[x,tmin,tmax]=read_adi_file(mydir,myfile,Fs,tdate_start,tlen,1);  %Output in voltage
        
        head.Nchan=1;
        
    case 'MDAT'
        
        if findstr(mydir,'Arctic_2010')  %%If we have two arrays to be merged together
            merge_two_mdat_files;
        else
            if strcmp(Ichan,'all')
                Ichan=1:8;
            end
            if max(Ichan)>8,
                disp('Channel too high, restricting to 8');
                Ichan=8;
            end
            [x,fparms]=read_mdat_file(mydir,myfile,tdate_start,tlen,1);
            Fs=fparms.fs;
            tmin=fparms.tfs;
            tmax=fparms.tfe;
            try
                x=x(Ichan,:);
            end
        end
        head.Nchan=size(x,1);
        head.calcurv=[
            -3.826740371096994e+07
            8.955964717194067e+07
            -9.041409824956748e+07
            5.195683344725087e+07
            -1.873830725900567e+07
            4.336813690670102e+06
            -6.251710636164501e+05
            5.339062841084976e+04
            -2.452316303554637e+03
            5.178423592184026e+01
            -1.229223450332553e+00];
        
    case 'WAV'
        sizz=wavread([mydir '/' myfile],'size');
        [~,Fs]=wavread([mydir '/' myfile],1,'native');
        handles.Fs=Fs;
        
        %%Can we calibrate the data?
        %%  load_wav often normalizes the data so the peak value is 1.
        %sens0=157; %dB re 1 unit of wav entry
        %sens=input('Enter sensitivity of entire system (flat-spectrum calibration) [180 dB re 1 unit]:');
        %if isempty(sens)
        %sens=sens0;
        %end
        %sens=10^(sens/20);
        
        %%Calibration for ADAT 24 attached to Sonotech hydrophone array
        %//ADAT HD24 has 6.9 V Rms max input for 24 bit data
        %//Factor of 2 from differential inputs...
        %//Hydrophone gain set to Sonotech array -157 dB
        %// Don't have cable attenuation  here...
        sens=(6.9*sqrt(2)/16777215)*0.5*10^(157.0 / 20.0);
        %freq_cal[i]=(1.0+2*Math.PI*f*(110e-9)*140.0);
        
        % Nov 15, 2011
        %         Hi Aaron,
        %
        % I thumbed through my notes for a few minutes to refresh my memory on this project.  The attenuation problem can be hugely simplified by eliminating many of the parameters right away.  The preamp output impedance is ordinarily very small, and the cable insulation resistance and receiver input impedance will be very large.  All of the above parameters can generally be ignored.
        %
        % This really only leaves the cable resistance and capacitance and so becomes a fairly simple voltage divider problem.  From my notes, I measured the cable (the entire 800+ length) resistance to be ~70 ohms and the capacitance to be ~110 nF.  The model would look like a series resistance with a shunt capacitance, so:
        %
        % R = 140 ohms (round trip)
        % X = 1/(2*pi*f*C)
        %
        % So, the 6 dB down frequency will be when R= X, so:
        % f = 1/(2*pi*C*R) = ~10 kHz.
        %
        % And the attenuation at any frequency:
        % X/(R + X)  or 1/(1+ 2*pi*f*C*R)
        %
        % The cable attenuation seen by the vector sensor will essentially look like a single pole low pass filter with Fc = 10 kHz:
        % 5   kHz: -3 dB
        % 10 kHz: -6 dB
        % 20 kHz: -9 dB
        %
        % This is consistent with my notes where I measured the cable attenuation to roll off by 3 dB/octave starting at 4-5 kHz.  The above is for the case of the vector sensor driving the entire length of the array.  For the other hydrophones driving shorter lengths of cable, the calculations will be the same but you will need to use different R and C values in the equation 20 log (1/(1+ 2*pi*f*C*R)) to compute the new attenuation values vs frequency.
        %
        % Jeff
        
        head.cable_factor=2*pi*(110e-9)*140.0;  %Unit resistance 140 ohm, capacitance 110 nF
        
        if tstart_min<0
            try
                tmin=convert_date(myfile,'_');
                tmax=convert_date(myfile,'_')+datenum(0,0,0,0,0,sizz(1)/Fs);
            catch
                disp([myfile ': convert_date failure']);
                minn=input('Enter start date in format [yr mo day hr min sec]: ');
                if isempty(minn)
                    minn=zeros(1,6);
                end
                tmin=datenum(minn);
                tmax=tmin+datenum(0,0,0,0,0,sizz(1)/Fs);
            end
        else
            tmin=tstart_min;
            tmax=tmin+datenum(0,0,0,0,0,sizz(1)/Fs);
            tdate_vec=datevec(tdate_start-tmin);
            nsec=tdate_vec(6)+60*tdate_vec(5)+3600*tdate_vec(4);
            N1=1+round(nsec*handles.Fs);
            N2=N1+round(tlen*handles.Fs);
            [x,Fs]=wavread([mydir '/' myfile],[N1 N2],'native');
            if ~strcmp(Ichan,'all')
                x=x(:,Ichan);
            end
            
            t=(1:length(x))/Fs;
        end
        x=double(x)*sens;
        head.Nchan=size(x,2);
        
end


%%%Optional Teager-Kaiser filtering...
if teager
    %%Assume that x is in form [ channel time]
    x=x(:,2:end-1).^2-x(:,1:end-2).*x(:,3:end);
end
    function merge_two_mdat_files
        if ~findstr(mydir,'Bottom_Unit')
            errordlg('Need to select a bottom unit');
            return
        end
        
        
        thisdir=pwd;
        cd(mydir);
        [xbot,fparms]=read_mdat_file('.',myfile,tdate_start,tlen,1);
        Fs=fparms.fs;
        tmin=fparms.tfs;
        tmax=fparms.tfe;
        cd('../Top_Unit');
        [head.rd,head.D,Iorder,torigin,toffset,tdrift,Iwant]=get_VA_2010_depths_and_offset(2010,myfile);
        
        if tdate_start<0
            tdate_start=torigin;
        end
        
        elapsed_time=tdate_start-torigin;
        if elapsed_time>0
            etime=datevec(elapsed_time);
            nsec=3600*etime(4)+60*etime(5)+etime(6);
        else
            etime=datevec(-elapsed_time);
            nsec=3600*etime(4)+60*etime(5)+etime(6);
            nsec=-nsec;
        end
        
        toffset=toffset+tdrift*nsec/(3600*1000);
        
        [xtop,fparms_top]=read_mdat_file('.',myfile,tdate_start+datenum(0,0,0,0,0,toffset),tlen,1);
        
        x=[xbot; xtop];  %Start with bottom element, move to top.
        
        x=x(Iorder,:);
        
        if strcmp(Ichan,'all')
            Ichan=Iwant;
        end
        if max(Ichan)>16
            disp('Channel too high, restricting to 16');
            Ichan=16;
        end
        
        x=x(Ichan,:);
        x=flipud(x);  %Shallowest element now first element, matching head.rd
        cd(thisdir);
        disp('Hello');
        
        %%Short-circuit data load with synthesized file...
        
    end

end  %function load_data




function [y,t,head]=readGSIfile(rawfile,cbegin,tlen,nchan,formatt,calibrate)
%function [y,t,head]=readGSIfile(rawfile,cbegin,tlen,nchan,formatt,calibrate)
% Input Parameters:
% rawfile = Include extension;
% cbegin = Start time;
% tlen = Length of sample to load;
% nchan = Index of desired channel, 1 for sound;
% formatt = Describes time input, string 'ctime' or 'datenum';
% calibrate = String 'calibrate' to convert Volts to microP;

if findstr(rawfile,'.sio')
    
    [y,t,head]=readsiof(rawfile,cbegin,tlen,formatt);
    y=y-0.5*(2^16);  %Remove DC bias in A/D converter
    
elseif findstr(rawfile,'.gsi')
    [y,t,head]=readgsi(rawfile,cbegin,tlen,formatt);
    if isempty(y)  %request time is befine file start
        dt=cbegin-head.ctbc;
        tlen=tlen+dt;
        [y,t,head]=readgsi(rawfile,0,tlen,formatt);
        
    end
    y=y(nchan,:);
    if strcmp(calibrate,'calibrate'),
        y=y-0.5*(2^16);  %Remove DC bias in A/D converter
        
        y=y*(2.5/65535)*(10^(150/20));
    end
    y=y.';
    
    if (abs(size(y,1)-floor(head.Fs*tlen))>2),
        disp('End of file reached, setting y to empty');
        y=[];
    end
end

end


% readsgi
% matlab script to read and display sample of gsi file.
% whose name is in fn
%function [x,t,head]=readgsi(fn,ctstart,tlen,formatt);
%  fn-file name, including extension
%  ctstart: c-time of start time...
%  nsec: number of seconds to read in...
%  formatt: 'datenum' or 'ctime'


function [omi,t,head]=readgsi(fn,ctstart,tlen,formatt)

t=[];omi=[];
head=readgsif_header(fn);
if nargin==1,
    return
end
if nargin<4,
    formatt='ctime';
end

if strcmp(formatt,'datenum')&ctstart>0,
    %tfile_start=datenum(1970,1,1,0,0,head.ctbc);
    ctstart=86400*(ctstart-datenum(1970,1,1,0,0,0));
    
end
%fid = fopen(char(fn),'r','ieee-le');
fs=head.Fs*(1+head.tdrift/86400);
fid = fopen(char(fn),'r','ieee-be');

if head.ctbc>ctstart&ctstart>0,
    disp('Start time less than file start');
    return;
end

if ctstart<=0,
    disp('cstart less than or equal to zero; interpret as seconds offset from file start');
    tsec=abs(ctstart);
else
    tsec=ctstart-head.ctbc;
end

fseek(fid,512+head.nc*floor(tsec*fs)*2,-1);

omi = fread(fid,[head.nc floor(head.Fs*tlen)],'uint16');
fclose(fid);

t=(1:size(omi,2))/head.Fs;

%N = 128;
%spectrogram(omi,hanning(N),N/2,N,head.Fs,'yaxis') % gram
%end

if 1==0,
    fn='S508A0T20080819T072522.gsi';
    tstart=0;
    %tstart=datenum(2008,8,19,12,0,0);
    tsec=3;
    [x,t,head]=readgsi(fn,tstart,tsec,'datenum');
    for I=1:size(x,1),
        subplot(size(x,1),1,I);
        spectrogram(x(I,:),256,128,256,1000,'yaxis')
        
        
    end
    
end
end

% readgsif_header
% matlab script to read and display sample of 2008 gsi file.
% whose name is in fn
% head=readsiof_header(fn);
function head=readgsif_header(fn)

%fid = fopen(char(fn),'r','ieee-le');
fid = fopen(char(fn),'r','ieee-be'); %for s1sT
head.dlabel(1:10) = char(fread(fid,10,'uchar'));
%head.nc=fread(fid,1,'uchar');
head.contents=char(fread(fid,4,'uchar'));

head.nc=length(find(double(head.contents)>32));
%head.nc=4;
%head.nc=length(deblank(head.contents'));

fillx(1:50) = fread(fid,50,'uint8');

a = fread(fid,9,'double');
head.ctbc=a(1);
head.ctec=a(2);
head.tdrift=a(3);
head.Fs=a(4);
head.UTMX=a(5);
head.UTMY=a(6);
head.depth=a(7);
head.UTMZone=a(8);
head.brefa=a(9);
head.ts = fread(fid,1,'uchar');
head.tabs_start=datenum(1970,1,1,0,0,head.ctbc);

fclose(fid);
end
%%%%%%%%%%%%%load_mt_mult.m%%%%%%%%%%%%
%[x,t,Fs]=load_mt_mult(folders,tdate_start,tlength,tshift,ftype);
%Load time series from two recorders, adjusting for initial clock
%% drift between them, then make coherence function...
%%Basic data....
% folders=strvcat('Linus','Lucy'), etc.;
% tdate_start: datenumber of absolute time at which to load_data.  Instead the
%   program uses this time to select which file to start processing.
% tlength:   how many seconds of data desired
% tshift=[0.7034 0]; %Inital time shift applied to channels, relative to
% tdate_start, in sec.  Index order is same as folder rows
%  ftype: if exists, specifys type of file to look for...
%Output:
%   x: matrix of data, rows are time, columns are channels corresponding to
%   folder order

% Sept 10, 2004:  Added option to load auxiliary mt files as well as sound,
% eliminated tdrift_rate, since often computed on outside
%
% Sept. 21, 2007:  Corrected an error that assumes that the mt files
%  are contiguous...now have 'samples expected' vs. 'samples got'.

function [x,t,Fs]=load_mt_mult(folders,tdate_start,tlength,tshift,ftype);

Nel=size(folders,1);
goodel=1:Nel;
if ~exist('ftype'),
    ftype='Sound';
end

%%%% if tshift a variable, adjust the start_time for loading data%%%

%%Make sure there are no negative time lags, otherwise load_mt wont work
% if tshift(1)<0,
%     tshift=[0 -tshift(1)];
% end
%
% if tdrift_rate<0,
%     tdrift_rate=[0 -tdrift_rate(1)];
% end

myloc=pwd;
for Idir=1:Nel,
    cd(deblank(folders(Idir,:)));
    if nargin>3,
        tdate_start_file=(tdate_start+datenum([0,0,0,0,0,tshift(Idir)]));
    else
        tdate_start_file=tdate_start;
    end
    tdate_end_file=(tdate_start_file+datenum([0,0,0,0,0,tlength]));
    %dirs_all{Idir}=dir('*Sound*');
    dirs_all{Idir}=dir(['*' ftype '*.mt']);
    
    Nfile=0;
    starttime=0;
    endtime=Inf;
    readflag=0;
    xx=[];
    %%AARON bug fix..in case each directory has different file names
    head=read_mt_header(dirs_all{Idir}(1).name);
    nx=floor(tlength*head.Fs);
    
    for Ifile=1:length(dirs_all{Idir}),
        head=read_mt_header(dirs_all{Idir}(Ifile).name);
        %%Identify data file where data start
        if tdate_start_file>=head.tstart&tdate_start_file<head.tend,  %All in one location...
            %dirs{Idir}.startfile=dirs_all{Idir}(Ifile);
            telapsed=datevec(tdate_start_file-head.tstart);
            starttime=3600*telapsed(4)+60*telapsed(5)+telapsed(6);
            readflag=1;
            disp(sprintf('Element:%s start file: %s start time in file:%6.2f',folders(Idir,:),dirs_all{Idir}(Ifile).name,starttime));
        end
        %%Identify data file where data end
        if tdate_end_file>=head.tstart&tdate_end_file<=head.tend,
            %dirs{Idir}.endfile=dirs_all{Idir}(Ifile);
            telapsed=datevec(tdate_end_file-head.tstart);
            endtime=3600*telapsed(4)+60*telapsed(5)+telapsed(6);
            readflag=1;
            disp(sprintf('   end file: %s start time in file:%6.2f',dirs_all{Idir}(Ifile).name,endtime));
        end
        if head.tstart>tdate_start_file&head.tend<tdate_end_file,
            readflag=1;
            disp(sprintf('   mid file: %s ',dirs_all{Idir}(Ifile).name));
        end
        
        
        if readflag==1,
            [xxx,tttt,head]=load_mt(dirs_all{Idir}(Ifile).name,starttime,endtime-starttime); %Note initial time offset
            %Safety check--is the size of xxx what we expect it to be?
            tleft=datevec(head.tend-tdate_start);
            samples_expected=round(head.Fs*(tleft(:,6)+60*tleft(:,5)+3600*tleft(:,4)))
            samples_got=length(xxx)
            if isinf(endtime)&&(samples_expected>samples_got)
                disp('Trying to read to end of file, but samples are missing...');
                xx=[xx; xxx; zeros(samples_expected-samples_got,1)];
            else
                xx=[xx; xxx];
            end
            readflag=0;
            starttime=0;
            endtime=Inf;
            %keyboard
        end  %Ifile
    end
    x(:,Idir)=xx;
    
    cd(myloc);
end
Fs=head.Fs;
t=(1:size(x,1))/Fs;
end

function [x,t,head,tdate]=load_mt(fname,tstart,tsec,uncal_chc)
%function [x,t,head,tdate]=load_mt(fname,tstart,tsec,uncal_chc),
% fname: file name
% tstart: start of file in seconds
% tsec: length of segment desired in seconds
% uncal_chc; If exists, return uncalibrated values.
% All times in seconds

if isempty(findstr(fname,'.mt')),
    fname=[fname '.mt'];
end
head=read_mt_header(fname);
Fs=head.Fs;
if strcmp(head.signing,'S'),
    dattype=['int' num2str(8*head.wordsize)];
else
    dattype=['uint' num2str(8*head.wordsize)];
end
if strcmp(head.swapping,'U'),
    fid=fopen(fname,'r','ieee-be');
else
    fid=fopen(fname,'r','ieee-le');
end

if fid<0,
    
    keyboard;
end
fseek(fid,512,-1);  %Skip the 512 bytes

offset=head.wordsize*floor(Fs*tstart);
status=fseek(fid,offset,0);
x=fread(fid,floor(tsec*Fs)+1,dattype);
fclose(fid);
%Optional calibration...
N=(2.^head.samplebits)-1;
if ~exist('uncal_chc'),
    switch head.signing,
        case 'S',
            calmean=0.5*(head.calmin+head.calmax);
            x=calmean+x*(head.calmax-head.calmin)/N;
        case 'U',
            calmean=head.calmin;
            x=calmean+x*(head.calmax-head.calmin)/N;
    end
    x=x*1000; %Convert from milliPa to microPa
    
end
t=(0:(length(x)-1))/Fs;
tdate=datenum(head.year,head.month,head.day,head.hours,head.minutes,tstart+t+head.seconds+head.msec/1000);
end
%%%%%%%%%%%read_mt_header.m%%%%%%%
% Read in mt header information
%function head=read_mt_header(fname),

% Nov 21, 2004: Changed head.Fs so that result is rounded to nearest
% value...
function head=read_mt_header(fname)

if isempty(findstr(fname,'.mt')),
    fname=[fname '.mt'];
end
%fname='aa_Sound_2003_08120945.mt';
fid=fopen(fname,'r');

if fid>0,
    head.magicstring=read_char(fid,8);
    head.totalhdrs=read_char(fid,3);
    head.abbrev=read_char(fid,8);
    head.stationcode=read_char(fid,3);
    head.title=read_char(fid,82);
    head.month=read_num(fid,3);
    head.day=read_num(fid,3);
    head.year=read_num(fid,5);
    head.hours=read_num(fid,3);
    head.minutes=read_num(fid,3);
    head.seconds=read_num(fid,3);
    head.msec=read_num(fid,4);
    head.sampling_period=read_char(fid,15);
    head.Fs=round(1./str2num(head.sampling_period));
    head.samplebits=read_num(fid,3);
    head.wordsize=read_num(fid,2);
    head.typemark=read_char(fid,1);
    head.swapping=read_char(fid,1);
    head.signing=read_char(fid,1);
    head.caltype=read_char(fid,1);
    head.calmin=read_num(fid,15);
    head.calmax=read_num(fid,15);
    head.calunits=read_char(fid,40);
    head.recordsize=read_num(fid,6);
    head.sourcevers=read_char(fid,9);
    head.sourcesn=read_char(fid,16);
    head.tstart=datenum(head.year,head.month,head.day,head.hours, ...
        head.minutes,head.seconds+head.msec/1000);
    
    fseek(fid,0,1);
    byte_length=ftell(fid);
    sec_length=byte_length/(head.wordsize*head.Fs);
    head.tend=datenum(head.year,head.month,head.day,head.hours, ...
        head.minutes,head.seconds+head.msec/1000+sec_length);
else
    disp(sprintf('%s could not be opened',fname));
end
fclose(fid);
end
function nhout=read_num(fid,N)
nhout=str2num(char(fread(fid,N,'char')'));
end
function chout=read_char(fid,N)
chout=char(fread(fid,N,'char')');
end


function [x, tfs, tfe, fs]    =...
        read_dat_file(file_name, fs, tstart, ns, units_voltage, sens, raw)
%READ_DAT_FILE     Read in a DAT file
%
%   WARNING!  LOG file must be in same location.
%
% Function form:
%  function [x, tfs, tfe, fs]    =...
%        read_dat_file(file_name, fs, tstart, ns, units_voltage, sens, raw)
%
% Input Parameters:
%   fname           -   name of *DAT file.  MUST include DAT extension and full pathname!
%   raw             -   t/f, weather output should be raw data or scaled
%   fs              -   sampling frequency in Hz.
%   tstart          -   datenumber of desired start time.
%                           If zero, data are read from start of file
%   nsec            -   number of seconds to be pulled.
%   units_voltage   -   if 1, then output in terms of voltage,otherwise
%                            output in terms of uPa. Default acoustic uPa
%   sens            -   Relative scaling factor for acoustic units
%
% Output Parameters:
%   x               -   Scaled data vector
%   tfs             -   datenumber of the start of the file
%   tfe             -   datenumber of the end of the file
%   fs              -   Sampling frequency (in case overridden internally)
%
%
% Revision History:
%   Version 1, Sept 30, 2006-read in raw data
%   Version 2, Oct. 28 2006- correct for time lost from serial port
%       acquisition--for this deployment every five minutes acoustic data
%       not acquired for 4 seconds (plus some time to power up).
%       Thus I use a drift formula of 48 sec/hour..
%
%   Modified by Jit Sarkar
%   2012/12/04
%       Header comments adjusted/formatted for readability
%   2012/12/17
%       Modified parameter checking, adding ability to read whole file
%       Modified all calls to "exist" function to be explicit, e.g.
%       searching for variables, or files
%       Changing directories just to find log file unecessary
%       Applied code-analyzer suggested fixes
%       Adjusted code formatting/spacing for readability




%%  Input parameter checking, and default values
                                                    %   dB re 1uPa/V
if ~exist('sens', 'var') || isempty(sens),      sens    =	-172;	end 
if ~exist('units_voltage', 'var') || isempty(units_voltage),
    units_voltage   =   0;      end
if ~exist('fs', 'var'),                         fs      =   [];     end
if ~exist('ns', 'var'),                         ns      =   [];     end
if ~exist('tstart', 'var') || isempty(tstart),  tstart  =   0;      end
if ~exist('raw', 'var') || isempty(raw),        raw     =   false;  end


%%  Constant values and preassigned variables

VREF    =   2.5;
nc      =   1;

%   An input voltage of 0.05 V p-p produces a 1.2 V p-p at A/D
scale   =   (0.038/1.0)*VREF/(2^16-1);  


%%  Parse input file name and associated log file
[pathstr, name, ~]    =   fileparts(file_name);

%   Check that associated log file exists
log_file    =   fullfile(pathstr, [name '.log']);
if ~exist(log_file, 'file')
    log_file    =   fullfile(pathstr, [name '.LOG']); % check upper case ext
    if ~exist(log_file, 'file')
        error('%s does not exist in same directory as %s',...
            log_file,                               file_name); 
    end
end

fid	=	fopen(log_file, 'r');

line.fs     =   [];
tfs         =   0;    tfe=0;
sec_flag    =   0;

while 1
    tline	=	fgetl(fid);
    
    if ~ischar(tline), break, end
    %   ??  Date parsing can be done with inbuilt matlab functions
    if      strfind(tline, 'Start time')
        line.start  =   tline;
        tfs         =   parse_date(line.start);
        sec_flag    =   1;
    elseif  strfind(tline, 'Stop time')
        line.end    =   tline;
        tfe         =   parse_date(line.end);
        sec_flag    =   2;
    elseif  strfind(tline, 'Sample rate:')>0,
        Idot        =   1 + strfind(tline,':');
        Iend        =   strfind(tline,'Hz')-1;
        if isempty(Iend)
            Iend    =   strfind(tline,'0');
        end
        fs  =   str2double(tline(Idot:Iend(end)));
    elseif  strfind(tline,'Sensitivity:')>0,
        Idot        =   1 + strfind(tline,':');
        Iend        =   strfind(tline,'dB')-1;
        sens        =   str2double(tline(Idot:Iend));
        
    elseif  strfind(tline,'seconds')
        Is          =   strfind(tline,'seconds')-1;
        secs        =   str2double(tline(1:Is));
        if sec_flag == 1
            tfs     =   tfs + datenum(0,0,0,0,0,secs);
        elseif sec_flag == 2
            tfe     =   tfe + datenum(0,0,0,0,0,secs);
        end
        sec_flag = 0;
    end
end
fclose(fid);

if (tfs==0 || tfe==0)
	fprintf('Could not understand %s',log_file);
end


%%  Parse date/time values for reading binary file
if  isempty(fs)
    error('Sampling rate not provided, and not found in log file');
end
if  isempty(ns)     %   assume whole file is requested
    np  =   inf;
else
    np  =   ns * fs;
end


try
    if tstart > tfs,
        offset  =   etime(datevec(tstart), datevec(tfs));
%        offset  =   (offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file
        %disp(sprintf('%i second offset from file start',offset));
    elseif tstart > 0
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset  =   0;
    else
        disp('Tstart less than or equal to zero, will interpret tstart as elasped time in seconds');
        offset  =   abs(tstart);
    end
catch    %#ok<CTCH>
    disp('File name does not contain date, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset  =   tstart;
end

%%  Open binary data file and read in requeste points
fid     =   fopen(file_name,'r','ieee-be');
fskip   =   2*nc*round(fs*offset);
res     =   fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
%disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));

%   Default to raw output
x   =	fread(fid,[nc np],'*uint16');
fclose(fid);

if ~raw
    x   =   scale*double(x);
    
    %Crude acoustic level calibration
    if units_voltage~=1,
        disp('Acoustic calibration executed');
        x   =   x.*(10^(-sens/20));
    end
end



end


function [x,tfs,tfe,fs,x_raw]=read_adi_file(dirname,fname,fs,tstart,ns,units_voltage,sens)
%%%%%%read_adi_file.m%%%%%%%%%%%
% [x,tfs,tfe,fs,x_raw]=read_adi_file(dirname,fname,fs,tstart,nsec,units_voltage)
% Read in a ADIOS file, which is little endian instead of big endian, and different calibration.
%  WARNING!  LOG file must be in same location.
% dirname: directory where DAT and LOG file are located.
% fname, name of *ADI file.  MUST include ADI extension!
% fs, sampling frequency in Hz.
% tstart-datenumber of desired start time.  If zero, data are read from
%   start of file
% nsec-number of seconds to be pulled.
% units_voltage: if 1, then output in terms of voltage, otherwise output in terms of uPa.
%       Default acoustic uPa
% Output:
%   tfs: datenumber of the start of the file
%   tfe: datenumber of the end of the file
% Created Feb. 17, 2012

if ~exist('units_voltage'),
    units_voltage=0;
end
%fs=50000;
x_raw=[];
VREF=3.3;
nc=1;
if ~exist('sens')
    sens=-172; %dB re 1uPa/V
end
%Get start time of file from .log file...
mydir=pwd;
cd(dirname);
fpre=findstr(fname,'.ADI')-1;
log_name=[fname(1:fpre) '.log'];
if ~exist(log_name)
    log_name=[fname(1:fpre) '.LOG'];
    if ~exist(log_name)
        error('read_adi_file: %s does not exist in same directory as %s',log_name,fname);
    end
end
fid=fopen(log_name,'r');

line.fs=[];
tfs=0;tfe=0;
sec_flag=0;
sensitivity_flag=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if findstr(tline,'Start time')
        line.start=tline;
        tfs=parse_date(line.start);
        sec_flag=1;
    elseif findstr(tline,'End time')
        line.end=tline;
        tfe=parse_date(line.end);
        sec_flag=2;
    elseif findstr(tline,'Sample rate:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'Hz')-1;
        if isempty(Iend)
            Iend=findstr(tline,'0');
        end
        fs=str2num(tline(Idot:Iend(end)));
        
    elseif findstr(tline,'Sensitivity:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'dB')-1;
        sens=str2num(tline(Idot:Iend));
        sensitivity_flag=1;
    elseif findstr(tline,'seconds')
        Is=findstr(tline,'seconds')-1;
        secs=str2num(tline(1:Is));
        if sec_flag==1
            tfs=tfs+datenum(0,0,0,0,0,secs);
        elseif sec_flag==2
            tfe=tfe+datenum(0,0,0,0,0,secs);
        end
        sec_flag=0;
    end
end
fclose(fid);
np=ns*fs;


if tfs==0|tfe==0
    disp(sprintf('Could not understand %s',log_name));
end
if sensitivity_flag==0
    disp(sprintf('LOG file did not have sensitivity; using a default value of %6.2f',sens));
end

try
    if tstart>tfs,
        offset=datevec(tstart-tfs);
        offset=(offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file
        %disp(sprintf('%i second offset from file start',offset));
    elseif tstart>0
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset=0;
    else
        disp('Tstart less than or equal to zero, will interpret tstart as elasped time in seconds');
        offset=abs(tstart);
    end
catch
    disp('File name does not contain date, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset=tstart;
end

scale=(1/20)*VREF/(2^16-1);  % Tested with input signal 50 mv pk-pk at 1 kHz (50MV1000.ADI)
fid=fopen(fname,'r','ieee-le');

cd(mydir);
fskip=2*nc*round(fs*offset);
res=fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
%disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));
x=fread(fid,[nc np],'uint16');
fclose(fid);
if nargout>4
    x_raw=x;
end
x=scale*x;


%Crude acoustic level calibration
if units_voltage~=1,
    disp('Acoustic calibration executed');
    x=x.*(10^(-sens/20));
end
end

%Future calibration corrections


function edit_chan_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chan as text
%        str2double(get(hObject,'String')) returns contents of edit_chan as a double

end
% --- Executes during object creation, after setting all properties.
function edit_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function [tfs,tfe,fs,nc,tiltx,tilty,sensitivity]=load_mdat_header(dirname,fname)
cd(dirname);
fpre=findstr(lower(fname),'.mdat')-1;
log_name=[fname(1:fpre) '.LOG'];
fid=fopen(log_name,'r');
sensitivity=-160; %dB re V/uPa default

while 1
    tline=fgetl(fid);
    if ~ischar(tline),break,end
    if findstr(tline,'Channels:')>0,
        Idot=1+findstr(tline,':');
        nc=str2num(tline(Idot:end));
    elseif findstr(tline,'Sample rate:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'Hz')-1;
        fs=str2num(tline(Idot:Iend));
    elseif findstr(tline,'Sensitivity:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'dB')-1;
        sensitivity=str2num(tline(Idot:Iend));
    elseif findstr(tline,'Start time')>0,
        Idot=findstr(tline,'=');
        tfs=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        tfs=tfs+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'Stop time')>0,
        Idot=findstr(tline,'=');
        tfe=parse_date(tline(Idot:end));
        tline=fgetl(fid);
        Iend=findstr(tline,'seconds')-1;
        tfe=tfe+datenum(0,0,0,0,0,str2num(tline(1:Iend)));
    elseif findstr(tline,'Tilt-X')>0,
        Idot=findstr(tline,':')+1;
        Iend=findstr(tline,'Tilt-Y')-1;
        tiltx=str2num(tline(Idot:Iend));
        
        Idot=Iend+8;
        Iend=findstr(tline,'degrees')-1;
        tilty=str2num(tline(Idot:Iend));
        
    end
    
end

fclose(fid);

end

function tabs=parse_date(str)

Istart=findstr(str,'=');
year=str2num(str((end-4):end));
tm=datenum(str((end-13):(end-5)),14)-datenum('00:00:00',14);
day=str2num(str((end-15):(end-14)));
month=(str((end-19):(end-16)));
switch deblank(month),
    case 'Jan'
        mn=1;
    case 'Feb'
        mn=2;
    case 'Mar'
        mn=3;
    case 'Apr'
        mn=4;
    case 'May'
        mn=5;
    case 'Jun'
        mn=6;
    case 'Jul'
        mn=7;
    case 'Aug'
        mn=8;
    case 'Sep'
        mn=9;
    case 'Oct'
        mn=10;
    case 'Nov'
        mn=11;
    case 'Dec'
        mn=12;
end

tabs=tm+datenum(year,mn,day,0,0,0);

end

%%%%%%read_mdat_file.m%%%%%%%%%%%
% [x,fparms]=read_mdat_file(dirname,fname,tstart,nsec,calibrate_acoustic)
% Read in a DAT file.
%  WARNING!  LOG file must be in same location.
% dirname: directory where DAT and LOG file are located.
% fname, name of *DAT file.  MUST include DAT extension!
% fs, sampling frequency in Hz.
% tstart-datenumber of desired start time.  If zero, data are read from
%   start of file
% nsec-number of seconds to be pulled.
% calibrate_acoustic: if 1, then output in terms of acoustics, if 0 output in terms of voltage
%       if 2, output raw count.
%       Default is raw count
% Output:
%   tfs: datenumber of the start of the file

% Version 1, Sept 30, 2006-read in raw data
% Version 2, Oct. 28 2006- correct for time lost from serial port
% acquisition--for this deployment every five minutes acoustic data not
% acquired for 4 seconds (plus some time to power up).  Thus I use
%   a drift formula of 48 sec/hour..


function [x,fparms]=read_mdat_file(dirname,fname,tstart,ns,calibrate_acoustic)

if ~exist('calibrate_acoustic')
    calibrate_acoustic=2;
end
%fs=50000;
VREF=2.5;
Nmax=(2^16)-1;
bias=(Nmax+1)/2;
if calibrate_acoustic<2
    scale=(0.038/1.00)*VREF/Nmax;  %An input voltage of 0.05 V p-p produces a 1.2 V p-p at A/D
elseif calibrate_acoustic==2
    scale=1;
end

%nc=1;
%np=ns*fs;

%Get data from log file
mydir=pwd;

[tfs,tfe,fs,nc,tiltx,tilty,sens]=load_mdat_header(dirname,fname);

fparms=struct('tfs',tfs,'tfe',tfe,'fs',fs,'nc',nc,'tiltx',tiltx,'tilty',tilty);
fid=fopen([fname],'r','ieee-be');
cd(mydir);

try
    if tstart>tfs&&tstart<(tfe-datenum(0,0,0,0,0,ns))
        offset=datevec(tstart-tfs);
        offset=(offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file. DO NOT ROUND
        disp(sprintf('%i second offset from file start',offset));
    else
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset=0;
        
        %offset=4;
    end
catch
    disp('File does not contain wanted time, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset=tstart;
end

%scale=1;
fparms.nstart=round(offset*fs);
%scale=1;
fskip=2*nc*fparms.nstart;
np=round(fs*ns);
fparms.nsamples=np;

res=fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));
x0=fread(fid,[nc np],'uint16');
fclose(fid);
x=scale*(x0-bias);


%Crude acoustic level calibration
if calibrate_acoustic==1
    disp('Acoustic calibration executed');
    x=x.*(10^(-sens/20));
end

%Future calibration corrections

end


%calibrate_GSI_signal.m
% Convert raw A/D value of presure sensor in DASAR to uPa
% using sensitivity of 150 dB re 1uPa/V
% and peak A/D voltage of 2.5 V
function x=calibrate_GSI_signal(xin, keyword,RawFileName)

%keyboard
if strcmp(keyword,'short')||~isempty(findstr(keyword,'DASAR2007'))
    [numd,dend] = DASAR_Shell_2007_equalization(1000,0);
    filt.a=dend;
    filt.b=numd;
    amp_Scale = (2.5/65535)*(10^(149/20));
    Nchan=size(xin,2);
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
        %x(:,I) = amp_Scale*filter(filt.a,filt.b,xt);
    end
elseif strcmp(keyword,'short')||~isempty(findstr(keyword,'DASARC'))
    %     x=xin.*(2.5/(2^16-1));  %x in Volts
    Nchan=size(xin,2);
    %      %%This is a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    amp_Scale = (2.5/65535)*(10^(149/20));
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
    end
elseif strcmp(keyword,'filter')
    error('calibrate_GSI_signal: filter no longer a valid keyword...')
elseif ~isempty(findstr(keyword,'NorthStar08')),
    
    %%This has a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    if ~isempty(findstr(RawFileName,'NS08A0'))|~isempty(findstr(RawFileName,'NA08Cx'))
        amp_Scale = (2.5/65535)*(10^(134/20));
    else
        amp_Scale = (2.5/65535)*(10^(148.8/20));
    end
    %amp_Scale=1;
    x = amp_Scale*filter(filt.b,filt.a,xin);
    %
    %     The Northstar deployments consist of 14 units, 12 of which are DASAR-Cs (built in 2008), and the other two DASAR-As (built in 2003).
    % These DASAR-As are identical to the Liberty DASAR-As with which you are already familiar.
    % By location, the breakdown is as follows (all units are DASAR-C unless marked otherwise):
    %
    % NA08A0 SN45
    % NA08B0 SN51
    % NA08C0 SN36
    % NA08D0 SN37
    % NA08E0 SN48
    % NA08F0 SN47
    % NA08G0 SN65
    % NA08H0 SN52
    % NA08I0 SN49
    % NA08J0 SN50
    % NS08A0 SN2 DASAR-A
    % NS08B0 SN58
    % NS08C0 SN59
    % NA08Cx SN1 DASAR-A
    % Response equalization for both types of DASAR (in the band 10 to 450 Hz) is very similar, with the only difference being a scalar gain value.
    % The hydrophone to ADC gain is
    % -149 dB V/uPa @ 100 Hz for the DASAR-C
    % -134 dB V/uPa @ 100 Hz for the DASAR-A
    %
    % You have the code to equalization both of these. For your convenience I have copied the relevant email message below...
    %
    % ------------------------
    % Aaron -
    %
    % For equalization of the omni channel on the BP DASAR-Cs, at sample rate 1 kHz and good for 10 Hz to 450 Hz analysis, you can use the
    %same IIR filter and coefficient values as were used on the Liberty DASAR-As. The only major difference is a change in gross sensitivity:
    % the DASAR-Cs are about 1/5 as sensitive as the DASAR-As. The same equalization can be used for all DASAR-Cs. Measurements in the Greeneridge
    % in-air loudspeaker test box showed all 65 units to be fairly close in sensitivity, with deviations looking no worse than the uncertainty in the calibration itself.
    %
    % The DASAR-Cs start rolling off below about 4 Hz (with a high-pass characteristic), so I don't recommend doing any analysis below that frequency without different equalization. You should be OK with these coef values since you are staying between 10 and 450 Hz.
    %
    % This equalization as shown includes a second order high-pass filter with break freq 1 Hz. You have the code to change it to 10 Hz if desired.
    %
    % Again, only one line changes going from DASAR-A to DASAR-C...
    %
    % om = x * 2.5/65536; % convert from counts to V
    % oml = 10^(148.8/20)*om; % uPa    convert from V to uPa   (WAS: 134 for DASAR-As, IS: 148.8 for DASAR-Cs)
    %
    % % this cascades the integrator with a 1 Hz high pass to cut low freqs
    % % -0.2 db at 100 Hz (ideally 0 db)
    % % -8.3 db at 250 Hz (ideally 20*log10(100/250)= -7.96
    % b = [0.53503845807280  -0.98982114743468   0.37452692065096   0.08025576871092]; % 11/13/2005
    % a = [1.00000000000000  -2.99111429220165   2.98226788807059  -0.99115359586894]; % 11/13/2005
    % % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % oml =  filter(b,a,oml); % uPa     equalize
elseif ~isempty(findstr(keyword,'Liberty08')>0)
    %om = xin * 2.5/65535; % convert from counts to V
    %om = (10^(134/20))*om; % uPa    convert from V to uPa
    %
    %     % this cascades the integrator with a 1 Hz high pass to cut low freqs
    %     % -0.2 db at 100 Hz (ideally 0 db)
    %     % -8.3 db at 250 Hz (ideally 20*log10(100/250)= -7.96
    % filt.b = [0.53503845807280  -0.98982114743468   0.37452692065096   0.08025576871092]; % 11/13/2005
    %filt.a = [1.00000000000000  -2.99111429220165   2.98226788807059  -0.99115359586894]; % 11/13/2005
    
    
    %This has a 10 Hz high pass filter
    
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    amp_Scale = (2.5/65535)*(10^(134/20));
    %amp_Scale=1;
    x = amp_Scale*filter(filt.b,filt.a,xin);
    
    %     keyboard;
end

if 1==0,
    Fs=1000;
    Nfft=256;
    A=fft(filt.a,Nfft);
    B=fft(filt.b,Nfft);
    F=linspace(0.1,Fs,Nfft);
    F=F(1:(Nfft/2));
    W=20*log10(abs(B./A));
    W=W(1:(Nfft/2));
    freqz(filt.b,filt.a,Nfft,Fs)
    hold on;plot(F,W,'r')
    
    figure;semilogx(F,W);xlim([1 1000]);grid on
    
end
end

function  [a_eqC2, b_eqC2]=get_DASARA_filter(f_hp,plot_data)
% DASAR_A_equalization.m

%%
Fs = 1000; % sample rate, Hz

% equalizer rev C
% used to equalize Sparton DIFAR sonobuoy head omni phone response
% The raw sonobuoy response is like a differentiator, with a +19 to +20dB/dec slope.
% This equalizer is mostly an integrator, with some shaping near Nyquist to get rid
% of aliasing effects and a high-pass filter cascaded to reduce low frequency gain

f1 = 100; % Hz
a1 = 1; % gain @ f1 (ignoring high-pass), V/V

% integrator
a = [1 -1]; % gain of integrator is cos(w/2)/sin(w), w = [0,2*pi] around unit circle

% hf zero
k = (f1/Fs)*2*pi;
b = a1 * sin(k)/cos(k/2) * [1 0.15]/1.15; % 0.15 eyeball adhoc
%b=a1;
% b = a1 * sin(k)/cos(k/2);

% high-pass
%f_hp = 10; % high-pass -3 dB break freq, Hz
[bb,aa] = butter(2,f_hp/(Fs/2),'high');

a_eqC2 = conv(a,aa);
b_eqC2 = conv(b,bb);

if plot_data==1,
    ff = [logspace(-2,log10(500),1000)]';
    [h_eqC2,w] = freqz(b_eqC2,a_eqC2,ff*pi*(2/Fs));
    figure
    semilogx(w/pi *(Fs/2),20*log10(abs(h_eqC2)),'k')
    grid on
    xlabel('Frequency, Hz')
    ylabel('gain, dB V/V')
end

end

function [numd,dend] = DASAR_Shell_2007_equalization(Fs,plot_on)
% DASAR_Shell_2007_equalization
%       numd=b; dend=a;
%   Returns filter coefficients of a digital equalization filter to flatten
%   the (very) low-frequency response of data recorded by the
%   omnidirectional channel in the Shell 2007 DASAR.
%
%   The 2007 Shell DASAR was a 4-channel unit, with a PZT flexural disk
%   hydrophone and three orthogonal geophones. The lower band edge of the
%   hydrophone channel is controlled by two cascaded (and independent)
%   high-pass filters, one formed by the shunt resistor across the PZT
%   ceramic, the other by the preamp (a single non-inverting opamp stage)
%
%   This equalization is only needed if it is desired to reconstruct data
%   below 5 (or 8 Hz): for the "usual" (10 Hz and above) analysis, this
%   equalization to the recorded data is not necessary. The user is advised
%   that attempting to equalize more than a decade below these frequencies
%   (i.e., 0.5 Hz) should be attempted with caution, since ambient noise
%   can be large (especially pressure fluctuations from surface waves in
%   shallow water), and self-noise of the instrument can become visible.
%
%   Unlike the DASAR-A and DASAR-Cs, which employed modified sonobuoy heads
%   having a differentiator-like response (rising ~6 dB/octave), the
%   in-band response of all four channels of the Shell 2007 DASAR were
%   flat, so no equalization is needed in the passband of 10 to 450 Hz

% R Norman Mar 22, 2011

%%
if(Fs ~= 1000)
    error('designed & tested for Fs = 1000 Hz only')
end
%plot_on=1;
%% equalization filter for the high-pass filter formed by PZT ceramic and shunt resistor
R1 = 2e6; % shunt resistance, Ohms
C1 = 15e-9; % ceramic capacitance, F
f1a = 1/(2*pi*R1*C1); % break frequency, Hz

p1a = 2*pi*f1a; % zero location, rad/s

% pole location, rad/s (here, arbitrarily placed one-decade below the zero)
% Don't place this pole at a frequency any lower than necessary)
c = 2.5;
p1b = p1a/c; % pole location, rad/s

nums1 = c * [1/p1a 1]; % s-plane numerator coefs
dens1 = [1/p1b 1]; % s-plane denominator coefs

[numd1,dend1] = bilinear(nums1,dens1,Fs);

if(plot_on)
    hFVT = fvtool(numd1,dend1);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'First Equalization Filter')
    pause
end

%% equalization filter for the high-pass filter formed by preamplifier (non-inverting gain stage)
R2 = 200; % resistance, Ohms
C2 = 100e-6; % capacitance, F
f2a = 1/(2*pi*R2*C2);

p2a = 2*pi*f2a; % zero location, rad/s

% pole location, rad/s (here, arbitrarily placed one-decade below the zero)
% Don't place this pole at a frequency any lower than necessary)
p2b = p2a/c; % pole location, rad/s

nums2 = c * [1/p2a 1];
dens2 = [1/p2b 1];

[numd2,dend2] = bilinear(nums2,dens2,Fs);

if(plot_on)
    hFVT = fvtool(numd2,dend2);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'Second Equalization Filter')
    pause
end

%% cascade the two filters
numd = conv(numd1,numd2);
dend = conv(dend1,dend2);

% future fancier
% H1 = dfilt.df2(numd1,dend1); % df2t for transposed
% H2 = dfilt.df2(numd2,dend2); % df2t for transposed
% Hcascade = dfilt.cascade(H1,H2)
% [n1,d1] = tf(Hcascade)
% fvtool(Hcascade)

%%
if(plot_on)
    % freqress(nums,dens,logspace(-3,2,1000))
    hFVT = fvtool(numd,dend);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'Final (composite) Equalization Filter')
    orient landscape
    print -djpeg DASAR2007_equalization.jpg
end

end

% --- Executes on button press in pushbutton_GSIbearing.
function pushbutton_GSIbearing_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_GSIbearing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get_GSI_bearing(hObject,eventdata,handles);

end

function [thet,kappa,tsec]=get_GSI_bearing(hObject,eventdata,handles)
%tsec: seconds into time series that data are selected...
thet=-1;  %Start with failed result
kappa=-1;
tsec=-1;
tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
%Ichan=str2num(get(handles.edit_chan,'String'));
[x,t,Fs,tstart,tend,head]=load_data(handles.filetype,handles.min_time,tdate_start,tlen,'all',handles);

disp('Click on two extreme corners, click below axis twice to reject:');
tmp=ginput(2);
if (isempty(tmp)||any(tmp(:,2)<0))
    return
end
tsec=min(tmp(:,1));
n=round(Fs*sort(tmp(:,1)));

freq=1000*sort(tmp(:,2));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

[thet0,kappa,sd]=extract_bearings(x(n(1):n(2),:),0.25,Nfft,Fs,freq(1),freq(2),50);

if ~isempty(findstr('T2007',handles.myfile))
    cal07flag=1;
    handles.calibration_DASAR2007_dir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/RawData';
    
    brefa_table=calibrate_bearing_Shell2007(handles.calibration_DASAR2007_dir,handles.myfile);
else
    cal07flag=0;
end

if cal07flag==0
    thet=bnorm(thet0+head.brefa);
else
    [junk,Icol]=calibrate_bearing_Shell2007(handles.calibration_DASAR2007_dir,handles.myfile,1);
    thet= interp1(0:360,brefa_table(:,Icol),bnorm(thet0));
end

%handles.bearing=thet;
%guidata(hObject, handles);
%titlestr=get(handles.text_filename,'String');
%Iend=findstr(titlestr,'.gsi')-1;
set(handles.text_filename,'String',sprintf('%s/%s %6.2f degrees... ',handles.mydir,handles.myfile,thet));
%keyboard;
end

% --- Executes on button press in pushbutton_GSI_localization.
function pushbutton_GSI_localization_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_GSI_localization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NDasar=7;
theta=-2*ones(1,NDasar);
kappa=theta;
strr=upper('abcdefghijklm');
for Idasar=1:NDasar
    %%Assume final directoryname containing files is of the form 'S510G0/'
    handles.mydir(end-2)=strr(Idasar);
    %%S510G0T20100831T000000.gsi form
    handles.myfile(5)=strr(Idasar);
    disp(sprintf('Directory %s contains %s',handles));
    try
        load_and_display_spectrogram(hObject,eventdata,handles);
        [theta(Idasar),kappa(Idasar),tsec]=get_GSI_bearing(hObject,eventdata,handles);
    catch
        disp(sprintf('Directory does not exist'));
    end
    
    
end

prompt1={'File of locations','Site','UTM location to compute range'};
dlgTitle1='Parameters for GSI localization...';
def1={'/Volumes/ThodePortable2/2010_Beaufort_Shell_DASAR/DASAR_locations_2010.mat', '5','[4.174860699660919e+05 7.817274204098196e+06]'};
answer=inputdlg(prompt1,dlgTitle1,1,def1);
locs=load(answer{1});
Isite=str2num(answer{2});
Ikeep=find(~isnan(theta)&theta>0);  %Remove absent or unselected DASARs
Igood=find(~isnan(theta)&theta>-2);  %Remove absent DASARS (used to compute center of array)
%theta=theta(Ikeep);
%kappa=kappa(Ikeep);
%Note that the Site contains all DASAR locations,
%  whether they collected data or not
DASAR_coords=[locs.Site{Isite}.easting locs.Site{Isite}.northing];
[VM,Qhat,w0,outcome] = vmmle_r(theta(Ikeep)',DASAR_coords(Ikeep,:),'h',kappa(Ikeep)');
mean_coords=mean(DASAR_coords(Igood,:));
CRITVAL=4.60517; %chi2inv(0.90,2);

[AREA,A,B,ANG,Baxis] = ellipsparms(Qhat,CRITVAL,mean_coords,VM);

figure
Ikeep=find(~isnan(theta(Igood))&theta(Igood)>0);
[Dn,xg,yg]=plot_location(DASAR_coords(Igood,:),theta(Igood),Ikeep,VM,A,B,ANG);
hold on
VA_cords=str2num(answer{3})/1000;
tmp=(VA_cords-[xg yg]);
plot(tmp(1),tmp(2),'o');
range=sqrt(sum((VA_cords-VM/1000).^2));
title(sprintf('Range of source from chosen location: %6.2f km',range));
yes=menu('Print?','Yes','No');
if yes==1
    orient landscape
    tstart=datestr(datenum(0,0,0,0,0,tsec)+datenum(get(handles.edit_datestr,'String')),30);
    figure(1);
    print(1,'-djpeg',sprintf('Localization_G_%s.jpg',tstart));
    print('-djpeg',sprintf('Spectrogram_G_%s.jpg',tstart));
end
close(1)
end


%function [thet,kappa,sd]=extract_bearings(y,bufferTime,Nfft,Fs,fmin,fmax,Nsamples)
% Input:
%    y: time series, with channels arranged as columns
%    bufferTime: how much time exists before and after signal proper.
%           Needed because time-domain filtering requires a signal buffer.
%    Nfft: FFT size to use if CSDM is to be estimated.
%    Fs: Sampling frequency, Hz
%    fmin,fmax: minimum and maximum frequency of signal, Hz.
%    algchc:  String containing one of three possibilities:
%           'spectrogram: make spectrogram, estimate bearing from each
%                   time-frequency cell of contour, provide standard
%                   deviation
%           'sel_ratio':  time-domain method that takes rms value of
%           bandpassed signals.
%           'sel_ratio_FFT': takes the FFT of all signal components and
%               computes active intensity in frequency domain.
%           'CSDM':  With a small Nfft value, construct multiple shapshots
%               of signal structure, yielding a CSDM that in turn gives
%               and active intensity.  May also be modified to permit
%               robust beamforming for weak signals.
%    Nsamples:  number of bootstrap samples for estimating kappa and sd
%  Output: thet: angle in degrees, increasing clockwise from true north...,
%  sd standard deviation in degrees
function [thet,kappa,sd]=extract_bearings(y,bufferTime,Nfft,Fs,fmin,fmax,Nsamples)
kappa=[];sd=[];thet=[];

filter_chc='butter';

if strcmp(filter_chc,'FIR')
    transband=0.1*(fmax-fmin); %Hz
    filter_min=max([0.5*transband 0.8*fmin]);
    filter_max=min([500-0.5*transband 1.2*fmax]);
    [n,fo,mo,w] = firpmord( [filter_min+0.5*transband*[-1 1] filter_max+0.5*transband*[-1 1]], [0 1 0], [0.01 0.1 0.01], Fs );
    b = firpm(n,fo,mo,w);
    
    % design and apply filter to pass only between e1 and e2
    
    for I=1:3,
        y(:,I)=y(:,I)-mean(y(:,I));
        x(:,I)=filtfilt(b,1,y(:,I));
    end
    
else
    w1=max([0.01 fmin*2/Fs]);
    w2=min( [fmax*2/Fs 0.99]);
    
    [B,A]=butter(2,[w1 w2]);
    for I=1:3,
        y(:,I)=y(:,I)-mean(y(:,I));
        x(:,I)=filter(B,A,y(:,I));
    end
end
if bufferTime>0,
    x=x(ceil(1+bufferTime*Fs):(end-floor(bufferTime*Fs)),:);
end

vx=(x(:,1).*x(:,2));
vy=(x(:,1).*x(:,3));
%thet=atan2(sum(vx),sum(vy))*180/pi;


[thet,kappa,sd]=get_vmests([vx vy],Nsamples);


end

function [mu,kappa,sd] = get_vmests(x,B)
%GET_VMESTS Calculates bootstrap estimate of mean and standard error of bearings.
%  [MU,KAPPA,SD] = GET_VMESTS(X,B) assumes that X is an n-by-2 matrix of x-y
%  coordinate pairs and B is a scalar denoting the number of bootstrap iterations.
%  MU is the mean bearing, KAPPA is an estimate of the standard error of the mean
%  expressed as the von Mises concentration parameter, and SD is the standard error
%  estimate in degrees expressed as for a linear (not a circular variable).
%
%  Note: B may be reduced for greater speed.

n = size(x,1);
flag = 0;           % Set to 1 to calculate estimates of kappa and Rbar for each
%   bootstrap sample, 0 for mean vector only.
mux_hat = zeros(B,2);
options = optimset('Display','off');
mux = vm_ests_uneq(x,options,flag);             % Get mean angle from the data.
mu = 180/pi*atan2(mux(1),mux(2));
%tic
%for i = 1:B,
%  U = ceil(n .* rand(n,1));          % The guts of UNIDRND.M without error checking.
%  xb = x(U,:);
%  mux_hat(i,:) = vm_ests_uneq(xb,options,flag); % Estimation accounting for lengths.
%end

U = ceil(n .* rand(n,B));          % The guts of UNIDRND.M without error checking.
for i = 1:B,
    %U = ceil(n .* rand(n,1));          % The guts of UNIDRND.M without error checking.
    %xb = x(U,:);
    mux_hat(i,:) = vm_ests_uneq(x(U(:,i),:),options,flag); % Estimation accounting for lengths.
end

[junk,kappa,sd] = vm_ests(mux_hat,options);     % Estimation ignoring lengths.
%toc
end

function [mux,Rbar,kappa] = vm_ests_uneq(x,options,flag)
%VM_ESTS_UNEQ Maximum likelihood estimates of Von Mises parameters from data.
%  [MUX,RBAR,KAPPA] = VM_ESTS_UNEQ(X,OPTIONS,FLAG) estimates the mean vector (MUX),
%  the length of the mean vector (RBAR), and the concentration parameter (KAPPA)
%  of the Von Mises distribution.  X is assumed to be an n-by-2 matrix, where
%  each row represents an x,y coordinate pair, which defines a bearing from
%  origin (0,0).  OPTIONS is an argument for FMINBND, for instance as set by
%  OPTIMSET.  FLAG is a boolean: 0 indicates that only MUX will be calculated;
%  1 indicates that RBAR and KAPPA will also be calculated
%
%  Calculates mean vector length, and thus KAPPA, taking account of individual
%  vector lengths.

%n1 = size(x,1);

lx = sqrt(sum(x.^2,2));                % Get length of each vector
idx = find(lx);                        % Get rid of 0-length vectors
%n = length(idx);
%if n<n1,
x = x(idx,:);
%end
r = sum(x);
mux = r/sum(lx);                       % Mean x,y coordinate


end


function [mu,kappa,sd] = vm_ests(x,options)
%VM_ESTS Maximum likelihood estimates of Von Mises parameters from data.
%  [MU,KAPPA,RBAR,SD] = VM_ESTS(X,OPTIONS) estimates the angular mean (MU)
%  and concentration parameter (KAPPA) of the Von Mises distribution, the
%  mean vector length (RBAR), and the approximate standard deviation (SD)
%  of the normal distribution from data X.  X is assumed to be an n-by-2
%  matrix, where each row represents an x,y coordinate pair, which defines
%  a bearing from origin (0,0).  OPTIONS is an argument for FMINBND, for
%  instance as set by OPTIMSET.

n1 = size(x,1);
lx = sqrt(sum(x.^2,2));                % Get length of each vector
idx = find(lx);                        % Get rid of 0-length vectors
n = length(idx);
if n<n1,
    x = x(idx,:);
end
x = x./repmat(lx,1,2);                 % Make unit vectors
r = sum(x);
mux = r/n;                             % Mean x,y coordinate
mu = 180/pi*atan2(mux(1),mux(2));      % Mean angle (use y/x for math convention)
Rbar = norm(mux);                      % Length of mean vector
sd = 180/pi*sqrt(2*(1-Rbar));          % Angular standard deviation

kappa=0;
% if Rbar<=(1/sqrt(n)),
%     kappa = 0;
% else
%     kappa = fminbnd('diffkr',eps,5e5,options,Rbar,n);
% end
%     function d = diffkr(k,r,n)
%         %DIFFKR Difference of a function of Von Mises kappa and R, mean vector length.
%         %  D = DIFFKR(K,R) is used by Matlab function FMINBND to find the maximum
%         %  likelihood estimate of kappa, the concentration parameter of the Von Mises
%         %  distribution. The MLE of kappa is the value K such that |R - A(K)| is
%         %  minimum, where R is the length of the mean vector, A(K) =  I1(K)/I0(K), and
%         %  I1 and I0 are modified bessel functions of order 1 and 0, respectively.
%         %
%         %  A bias correction for small sample size (see Batschelet, 1981, p. 47) is
%         %  included so that K is sought for min{|R - A(K)/A(KRN)|} , where N is the
%         %  sample size.
%         %
%         %  Rather than call the Matlab function, BESSELI, faster polynomial
%         %  approximations from Abramowitz and Stegun (1965, p. 378) are used.  Also,
%         %  the approximations allow arguments > 700 (which cause overflow in BESSELI).
%         %  In the code below, I0 and I1 represent functions of I0(X) and I1(X),
%         %  respectively, depending on the value of the argument X.  For X>3.75, the
%         %  leading factor,  cancels in the ratio allowing calculation of A(X) without
%         %  numeric overflow.
%
%         %  Abramowitz, M. and I.A. Stegun. 1965.  Handbook of Mathematical Functions.
%         %  Dover Publications, New York.
%
%         %  Batschelet, E. 1981. Circular Statistics in Biology. Academic Press, London.
%
%         krn = k*r*n;
%         t = krn/3.75;
%         if krn<=3.75,
%             I0 = 1 + 3.5156229*t^2 + 3.0899424*t^4 + 1.2067492*t^6 + 0.2659732*t^8 + ...
%                 0.0360768*t^10 + 0.0045813*t^12;
%             I1 = krn * (0.5 + 0.87890594*t^2 + 0.51498869*t^4 + 0.15084934*t^6 + ...
%                 0.02658733*t^8 + 0.00301532*t^10 + 0.00032411*t^12);
%             Akrn = I1/I0;
%         else
%             I0 = 0.39894228 + 0.01328592/t + 0.00225391/t^2 - 0.00157565/t^3 + ...
%                 0.00916281/t^4 - 0.02057706/t^5 + 0.02635537/t^6 - 0.01647633/t^7 + ...
%                 0.00392377/t^8;
%             I1 = 0.39894228 - 0.03988024/t - 0.00362018/t^2 + 0.00163801/t^3 - ...
%                 0.01031555/t^4 + 0.02282967/t^5 - 0.02895312/t^6 + 0.01787654/t^7 - ...
%                 0.00420059/t^8;
%             Akrn = I1/I0;
%         end
%         t = k/3.75;
%         if k<=3.75,
%             I0 = 1 + 3.5156229*t^2 + 3.0899424*t^4 + 1.2067492*t^6 + 0.2659732*t^8 + ...
%                 0.0360768*t^10 + 0.0045813*t^12;
%             I1 = k * (0.5 + 0.87890594*t^2 + 0.51498869*t^4 + 0.15084934*t^6 + ...
%                 0.02658733*t^8 + 0.00301532*t^10 + 0.00032411*t^12);
%             A = I1/I0;
%         else
%             I0 = 0.39894228 + 0.01328592/t + 0.00225391/t^2 - 0.00157565/t^3 + ...
%                 0.00916281/t^4 - 0.02057706/t^5 + 0.02635537/t^6 - 0.01647633/t^7 + ...
%                 0.00392377/t^8;
%             I1 = 0.39894228 - 0.03988024/t - 0.00362018/t^2 + 0.00163801/t^3 - ...
%                 0.01031555/t^4 + 0.02282967/t^5 - 0.02895312/t^6 + 0.01787654/t^7 - ...
%                 0.00420059/t^8;
%             A = I1/I0;
%         end
%         d = abs(r - A./Akrn);
%
%     end


end

function bearings=bnorm(bearings)
Ibig=find(bearings>=360);
bearings(Ibig)=bearings(Ibig)-360;

Ilow=find(bearings<0);
bearings(Ilow)=bearings(Ilow)+360;

end
% %tic
% U=ceil(n.*rand(n,B));
% for I=1:B
%    mux_hat(I,:)=vm_ests_uneq(x(U(:,I)),options,flag);
% end
% [junk,kappa,sd] = vm_ests(mux_hat,options);     % Estimation ignoring lengths.
% %toc


% --- Executes on button press in pushbutton_CSDM.
function pushbutton_CSDM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_CSDM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
mydir=pwd;
Ichan='all';  %Hardwire first channel
[x,t,Fs,tstart,junk,head]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,Ichan,handles);
x=x.';

%%For MDAT file, channel 1 is already shallowest channel, so remove data below...
%x=fliplr(x);
%head.rd=fliplr(head.rd);

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;
fmin=str2num(get(handles.edit_fmin,'String'));
fmax=str2num(get(handles.edit_fmax,'String'));
chann=input('Enter channel indicies (remember that channels have been flipped, 1 is now shallowest) (1:Nchan):');
if isempty(chann)
    chann=1:head.Nchan; %For 2008 short-term vertical array
end

Ichc=menu('Enter type of frequency processing:','Contour','All');

if Ichc==1
    Nf=input('Enter number of harmonics:');
    if isempty(Nf)
        Nf=1;
    end
    
    ftmp=ginput(Nf);
    frange=ftmp(:,2)*1000;
    
    Nbad=input('Enter number of bad frequencies: ');
    if ~isempty(Nbad)
        fbad=ginput(Nbad);
        fbad=fbad(:,2)*1000;
    else
        fbad=[];
    end
    fr=input('Enter frequency band to search (10):');
    if isempty(fr)
        fr=10;
    end
    
    [Ksout,Ns]=extractKsbest_contour(x,ovlap,Nfft,chann, frange,fr,fbad,Fs,Nfft,'keep_zeros');
    %    Ksout: structure array containing
    %     Kstot CSDM size(Nel,Nel,Nfreq)
    %     freq: frequencies corresponding to Ks in third dimension
    %     fcontour: frequencies detected, (Is,If)--can use to make contours..
    %     fcount: Number of times each frequency has been averaged..
    %     fpower: Power in each bin
    
    %[Kstot,Ks_eig,Nsnap,freq_all,pwr_est,SNRest]=extractKsexact_MDAT(x,ovlap,Nfft,chann,1000*[fmin fmax],Fs,-1,Nfft);
    
    axes(handles.axes1)
    hold on
    plot(Ksout.tcontour,Ksout.fcontour/1000,'o');
    hold off
elseif Ichc==2  %extractKsexact
    %%process all frequencies over a given range
    disp('Select two frequencies from image:')
    ftmp=ginput(2);
    frange=sort(ftmp(:,2)*1000);
    disp(sprintf('frange: %6.2f to %6.2f Hz',frange));
    threshold=input('Enter threshold in dB:');
    if isempty(threshold)
        threshold=-Inf;
    end
    [Ksout.Kstot,Ksout.freq,Ksout.VV,Ksout.EE]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,-1,Nfft,0,threshold);
    if ~isempty(Ksout.EE)
        Ksout.SNR=Ksout.EE(1,:)./Ksout.EE(2,:);
    else
        Ksout.SNR=Inf*ones(size(Ksout.freq));
    end
end

figure
plot(Ksout.freq,10*log10(Ksout.SNR),'-x');ylabel('dB SNR');xlabel('Freq (Hz)')
grid on;title('estimated SNR of CSDM frequency components')


yes=menu('Beamform?','No','Conventional','MV','Both','Reflection Coefficient Estimation','MFP');

if yes>1
    
    if yes<6
        angles=input('Enter vector of angles (-90:90):');
        if isempty(angles)
            angles=-90:90;
        end
    end
    switch yes
        case 2
            B=conventional_beamforming(Ksout.Kstot,angles,Ksout.freq,head.rd,1495);
        case 3
            B=MV_beamforming(Ksout.Kstot,angles,Ksout.freq,head.rd,1495);
        case 4
            B=conventional_beamforming(Ksout.Kstot,angles,Ksout.freq,head.rd,1495);
            
            B2=MV_beamforming(Ksout.Kstot,angles,Ksout.freq,head.rd,1495);
        case 5
            R=derive_reflection_coefficient2(Ksout.Kstot,angles,Ksout.freq,head.rd,1495);
        case 6
            %Simple Pekeris waveguide determined from DASAR cutoff frequencies
            
            prompt1={'Model file..','tilt offsets (deg)','ranges (m)', 'depths (m):','plot intermediate images?', ...
                'Nfft for plotting:','SNR_cutoff or list of frequencies..','Primary Eigenvector only?','receiver depths:'};
            
            dlgTitle1='Parameters for matched field processing...';
            
            [MFP_replica_file,range_str,default_tilt,default_SNR]=load_MFP_scenario;
            
            if length(MFP_replica_file)>1
                error('Cannot select multiple environmental models during MFP');
            end
            head.rd=get_VA_2010_depths_and_offset(2010,handles.myfile,MFP_replica_file{1});
            
            def1={MFP_replica_file{1}, default_tilt{1},range_str{1},'1:1:55','0','2048',['[' num2str(default_SNR{1}) ']'],'yes',num2str(head.rd)};
            answer=inputdlg(prompt1,dlgTitle1,1,def1);
            model_name=answer{1};
            tilt_offset=str2num(answer{2});
            ranges=str2num(answer{3});
            depths=str2num(answer{4});
            plot_chc=str2num(answer{5});
            Nfft2=str2num(answer{6});
            SNRmin=str2num(answer{7});  %May also be frequency
            head.rd=str2num(answer{9});
            Ksout.Nfft=Nfft;
            Ksout.ovlap=ovlap;
            
            if ~isempty(Ksout.VV)
                V1=squeeze(Ksout.VV(:,1,:));
                for If=1:length(Ksout.freq)
                    Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';
                    
                end
            else
                Ksout.SNR=Inf*ones(size(Ksout.freq));
                Ksout.Kstot_eig=[];
            end
            
            srcname=sprintf('CSDM_Nfft%i_%s',Nfft,datestr(tdate_start,30));
            
            if strcmp(answer{8},'yes')&&~isempty(Ksout.VV)
                save_result=matched_field_processor(model_name,tilt_offset,ranges,depths,Ksout.Kstot_eig,Ksout.freq,10*log10(Ksout.SNR),head.rd,SNRmin,srcname);
            else
                save_result=matched_field_processor(model_name,tilt_offset,ranges,depths,Ksout.Kstot,Ksout.freq,10*log10(Ksout.SNR),head.rd,SNRmin,srcname);
            end
            
            
            if ~isempty(save_result)
                save(srcname,'Ksout','head','tdate_start','tlen');
                write_covmat(Ksout.freq,Ksout.Kstot,head.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
                if ~isempty(Ksout.VV)
                    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
                    write_covmat(Ksout.freq,Ksout.Kstot_eig,head.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
                end
            end
            return
    end
    figure
    if yes==4
        subplot(2,1,1)
    end
    imagesc(Ksout.freq,angles,10*log10(B'));
    colorbar
    cmap=colormap;
    caxis([60 100])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Frequency (Hz)');ylabel('Angle from horizontal (deg)');grid on;
    title(sprintf('%s, %i FFT, %i elements',datestr(tdate_start),Nfft,length(head.rd)));
    set(gcf,'colormap',cmap(1:4:64,:));
    
    if yes==4
        subplot(2,1,2);
        imagesc(Ksout.freq,angles,10*log10(B2'));
        colorbar
        %cmap=colormap;
        caxis([60 100])
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Frequency (Hz)');ylabel('Angle from horizontal (deg)');grid on;
        title('MV processor');
        %set(gcf,'colormap',cmap(1:4:64,:));
        
    end
    orient tall
    print(gcf,'-djpeg',sprintf('BeamformingExample_%s',datestr(tdate_start,30)));
    keyboard
    
end

% yes=input('Type ''y'' to write a file:','s');
% if strcmp(yes,'y')
%     disp(['Printing ' srcname]);
%
%     save(srcname,'Ksout','head');
% end
%Save only non-zero CSDM to *.in file
if Ichc==2
    
    Ksout.Nfft=Nfft;
    Ksout.ovlap=ovlap;
    V1=squeeze(Ksout.VV(:,1,:));
    for If=1:length(Ksout.freq)
        Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';
        
    end
    
    srcname=sprintf('CSDM_Nfft%i_%s',Nfft,datestr(tdate_start,30));
    save(srcname,'Ksout','head','tdate_start','tlen');
    write_covmat(Ksout.freq,Ksout.Kstot,head.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    
    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
    write_covmat(Ksout.freq,Ksout.Kstot_eig,head.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    
    %     srcname=sprintf('CSDM_Nfft%i_%s_flipped',Nfft,datestr(tdate_start,30));
    %     write_covmat(Ksout.freq,Ksout.Kstot_eig,head.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    %
    figure
    imagesc(Ksout.freq,[],10*log10(abs(Ksout.EE)));
    xlabel('frequency (Hz)')
    ylabel('eigenvalue:')
    title(sprintf('eigenvectors of %s',srcname))
    orient landscape
    for I=1:3
        figure(I+1)
        imagesc(abs(squeeze(Ksout.VV(:,I,:))));
        
    end
    
end

keyboard;

% figure;
% axx=plotyy(freq_all,SNRest,freq_all,pwr_est);grid on;
% %set(axx(1),'xlim',xlimm)
% %set(axx(2),'xlim',xlimm)
% set(axx(1),'ylim',[0 40]);
% set(axx(1),'ytick',[0 10 20 30]);
% xlabel('Frequency (Hz)');
% ylabel('dB SNR');

end






%%%%%%%%%%%%%%%%%%%%%%%extractKsbest_contour.m%%%%%%%%%%%
% [Ksout,Ns]=extractKsbest_contour(x,ovlap,Nfft,goodel,frange,fr,fbad,Fs,M,nowin);
%  Aaron Thode
%  April 2, 2004
%  Generates averaged cross-spectral outputs
%    from an FM contour
% INPUT:
%x=array of data, rows are time, columns are channels
%  ovlap=overlap of time samples
%  Nfft-number of points used in FFT
%  chann-vector containing element indicies: referenced to *bottom* element
%  frange-Vector of frequencies: first is initial start of contour,
%   second is the start of second harmonic, etc.
%  fr-Search space in terms of +-freq
%  fbad-frequencies of constant interference, etc.
% M time window
% nowin-if exists, don't window the data before using fft.used for source signatureestimates
% keep_zeros: if exists, keep frequency bins that have no samples.  Useful for plotting...
% OUTPUT:
%    Ksout: structure array containing
%     Kstot CSDM size(Nel,Nel,Nfreq)
%     freq: frequencies corresponding to Ks in third dimension
%     fcontour: frequencies detected, (Is,If)--can use to make contours..
%     fcount: Number of times each frequency has been averaged..
%     fpower: Power in each bin
%The replica must also be conjugated to work properly with Kraken.
%	Don't know why.
% Feb 9-remove bad elements before summing power
% April 1, 2004-normalize by Fs*Nfft to put units as power spectral density

function [Ksout,Ns,EE_sort,VV]=extractKsbest_contour(x,ovlap,Nfft,chann, frange,fr,fbad,Fs,M,keep_zeros,nowin)
EE_sort=[];VV=[];
figure;
if ~exist('nowin'),
    nowin=0;
elseif  nowin==1
    nowin=1;
else
    nowin=0;
end

Nel=length(chann);
if ~exist('M')|M<0, M=Nfft;end

%Compute number of time snapshots
Ns=floor(((size(x,1)/M)-ovlap)/(1-ovlap));
disp(['Ns: ' int2str(Ns)]);
%if Ns==0,Ns=1;end

f=linspace(0,Fs,Nfft+1);
df=diff(f);df=df(1);
nbins=ceil(fr/df);
bins=-nbins:nbins;

fcount=zeros(size(f));
fpower=fcount;
Kstot=zeros(Nel,Nel,length(f));
%Kstot=zeros(Nel,Nel,length(findex));

%Select appropriate frequency bin
for I=1:length(frange),
    [junk,findex(I)]=min(abs(f-frange(I)));
    frange(I)=f(findex(I));
end

for I=1:length(fbad)
    [tmp,findexjunk(I)]=min(abs(f-fbad(I)));
    fbad(I)=f(findexjunk(I));
end

if Ns<0,
    disp('Signal too short for one shapnot, will center pad for FFT:');
    %pause;
    Ns=1;Nx=size(x,1);
    x0=zeros(Nfft,size(x,2));
    index=floor(Nfft/2-(Nx/2));
    index=(index:(index+Nx-1));
    x0(index,:)=x;
    x=x0;
    clear x0
    M=size(x,1);
end
if nowin==0,
    win=kaiser(M,2.5);
else
    win=ones(M,1);
end

%Determine the frequency with greatest average power near your bin!
Pt=[];
t=[];
for I=0:(Ns-1),
    index=round(I*M*(1-ovlap)+1);
    t(I+1)=index(1)/Fs;
    xindex=(index:(index+M-1));
    xh=x(xindex,chann);
    for Ic=1:size(xh,2),
        xh(:,Ic)=xh(:,Ic)-mean(xh(:,Ic));
        xh(:,Ic)=xh(:,Ic).*win;
    end
    Xh=fft(xh,Nfft);
    Pwr=abs(Xh).^2;
    Pt=sum(Pwr,2);
    if length(fbad)>0
        Pt(findexjunk)=0; %Remove bad freqencies.
        Pt(findexjunk+1)=0;
        Pt(findexjunk-1)=0;
    end
    %Pt=cat(2,Pt,Pwr);
    %end
    %Pt=sum(Pt,2);
    for If=1:length(findex),
        subplot(length(findex),1,If);
        
        %plot(f(findex(If)+bins),Pt(findex(If)+bins));
        [junk,fi]=max(Pt(findex(If)+bins));
        findex(If)=findex(If)+bins(fi);
        fcount(findex(If))=fcount(findex(If))+1;
        fcontour(I+1,If)=f(findex(If));
        fpower(findex(If))=fpower(findex(If))+Pt(findex(If));
        disp(f(findex(If)));
        pgoal=Xh(findex(If),:);
        %Make pgoal a vertical array
        %pgoal=conj(pgoal);   %Test to see if conjugation is the problem
        %   PREVIOUS STATEMENT COMMENTED OUT BY A THODE MARCH 15, 2004
        %     HE ALSO CHANGED write_covmat.m to remove conjugation as well
        
        pgoal=pgoal.';        %Rotate so vector is vertical, like KRAKEN
        Kstot(:,:,findex(If))=Kstot(:,:,findex(If))+pgoal*pgoal';
        %Ksout.pgoal
    end
    %title(int2str(Ns));
    %pause(0.25);
    disp('');
end

%%Collapse Kstot to non-zero components if desired
Igood=find(fcount>0);

if exist('keep_zeros') %%Keep all frequency bins, even if no power...
    %Igood=[min(Igood):max(Igood)];
    Igood=1:length(fcount);
end
fcount=fcount(Igood);
fpower=fpower(Igood);
freq=f(Igood);
Kstot=Kstot(:,:,Igood);

%Normalize by sample size, if necessary
Iavg=find(fcount>1);
for I=1:length(Iavg),
    Kstot(:,:,Iavg(I))=Kstot(:,:,Iavg(I))/fcount(Iavg(I));
    fpower(Iavg(I))=fpower(Iavg(I))/fcount(Iavg(I));
end

%Normalize the |X(f)|^2 term to make units power spectral density (power
%per Hz)
Kstot=Kstot/(Fs*Nfft);
%keyboard;

if Ns>1
    for J=1:length(freq)
        Ks=squeeze(Kstot(:,:,J));
        [V,D,FLAG]=eigs(Ks,4,'LM');
        %keyboard;
        D=real(diag(D));
        %%sqrt(D(1)*V(:,1)) gives scaled eigenvector
        %D_sort=sort(D);
        %D=flipud(D_sort)
        pgoal(:,J)=V(:,1)*sqrt(abs(D(1)));
        SN(1,J)=10*log10(abs(D(1)/D(2)));
        if fcount(J)>0
            disp(['Est. S/N for ' num2str(freq(J)) ' is ' num2str(SN(1,J)) 'dB']);
        end
    end
end

Ksout.Kstot=Kstot;
Ksout.freq=freq;
Ksout.fcount=fcount;
Ksout.fpower=fpower;
Ksout.fcontour=fcontour;
Ksout.tcontour=t;
Ksout.pgoal=pgoal;
Ksout.SN=SN;

close;
end


%%%%%%%%%%%%%%%%%%%%%%%extractKsexact.m%%%%%%%%%%%
%function [Kstot,Ks_eig,Ns,f,power]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin);
%[Kstot,pgoal,Ns,f,pwr_est,SNRest,pphase]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin);
%  Aaron Thode
%  July 3, 1996
%  Generates averaged cross-spectral outputs
%    from data input
%  CRUCIAL:  There is a line that flips pgoal so that the first
%		element represents the top phone.
%x=array of data, rows are time, columns are channels
%ovlap=overlap of time samples
%Nfft-number of point desired taken
%chann-vector containing element number to be used:
%frange-pairs of frequencies that define desired ranges If odd, last
%	element is bin spacing (i.e.) evaluate every other bin
%Isnap %Select the Isnap window.  If negative, average
%M amount of dataused for the fft snapshot
% nowin-if exists, don't window the data before using fft.used for source signature estimates
% threshold-dB threshold of power (sum of energy across all frequencies and channels) to reject a snapshot.
%    Set to Inf to ensure all snapshots used.
% tiltdata: tilt: estimated vertical tilt of array in degrees.
%           rd: element depths in m
% Output: Kstot-CSDM
%          %V: eigenvectors of CSDM
%          %ev: eigenvalues of CSDM
%  April 1, 2004: normalize CSDM to have units of power spectral density:
%  |X(f)|^2/(Fs*Nfft), Power/Hz

function [Kstot,f,VV,EE_sort]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin,threshold,tiltdata)

VV=[];
EE_sort=[];
MAXF=Inf;
if ~exist('threshold','var')
    threshold = -Inf;
end

if ~exist('tilt','var')
    tiltdata.tilt=0; %
    tiltdata.rd=ones(length(chann),1);
end
if size(tiltdata.rd,2)>1
    tiltdata.rd=tiltdata.rd';
end
sintilt=sin(tiltdata.tilt*pi/180);
if ~exist('nowin','var')
    nowin=0;
    disp('The signal will be windowed');
elseif nowin==1,
    nowin=1;
    disp('The signal will NOT be windowed');
else
    nowin=0;
    disp('The signal will be windowed');
end

Nel=length(chann);
Ns=floor(((size(x,1)/M)-ovlap)/(1-ovlap));
disp(Ns);
%pause;
%Select appropriate frequency bin
%[f,findex]=makefaxis(Nfft,Fs,frange);
frange=round((Nfft/Fs)*frange);
findex=max([1 frange(1)]):min([Nfft/2 frange(2)]);
f=findex*(Fs/Nfft);

if length(findex)<MAXF
    Kstot=zeros(Nel,Nel,length(findex));
else
    disp('frange too long to make Kstot\n just making pgoal');
end
power=zeros(Ns,1);
Nf=length(findex);
if Isnap<0|~exist('Isnap'),
    Isnap=0:Ns-1; %average all
end

if Ns<=0,
    disp('Signal too short for one shapnot, will center pad for FFT:');
    %pause;
    Ns=1;Isnap=0;Nx=size(x,1);
    x0=zeros(Nfft,size(x,2));
    index=floor(Nfft/2-(Nx/2));
    index=(index:(index+Nx-1));
    x0(index,:)=x;
    x=x0;
    clear x0
    M=size(x,1);
end
if nowin==0,
    win=kaiser(M,2.5);
else
    win=ones(M,1);
end

for I=Isnap
    index=round(I*M*(1-ovlap)+1);
    xindex=(index:(index+M-1));
    xh=x(xindex,chann);
    for Ic=1:length(chann)
        xh(:,Ic)=xh(:,Ic)-mean(xh(:,Ic));
        xh(:,Ic)=xh(:,Ic).*win;
    end
    Xh=fft(xh,Nfft);
    pgoal=Xh(findex,:);
    %Make pgoal a vertical array
    %pgoal=pgoal(:,chann);   %Reject elements ;
    %pgoal=conj(pgoal);   %Test to see if conjugation is the problem
    pgoal=pgoal.';        %Columns are now single-frequency array snapshots
    %pgoal=flipud(pgoal);  %Puts topmost element first, according to lewis
    
    
    power(I+1)=(Fs/Nfft)*sum(sum(abs(pgoal).^2))/(Fs*Nfft);
    
    if length(findex)<MAXF&&threshold<=10*log10(abs(power(I+1)))
        for J=1:Nf
            %disp(f(J));
            
            tiltvec=exp(1i*2*pi*f(J)*sintilt*tiltdata.rd/1500);
            ptemp=pgoal(:,J).*tiltvec;
            Kstemp=ptemp*ptemp'; %Top LH cornertop element autocor
            Kstot(:,:,J)=Kstot(:,:,J)+Kstemp;
            
        end
    else
        power(I+1)=NaN;
    end
end  %I=Isnap

figure
tt=Isnap*(1-ovlap)*Nfft/Fs;
plot(tt,10*log10(abs(power)));
grid on
if ~isinf(threshold)
    hold on
    line([min(tt) max(tt)],threshold*[1 1]);
end
xlabel('Time (s)');ylabel('dB power');

if length(findex)>MAXF,
    Kstot=[];
end

%Normalize to have units of power spectral density...pwr per Hz
Kstot=Kstot/(Fs*Nfft);

%Compute eigenvalues of CSDM, and the ratio of the high-power to the
%low-power sound be like SNR


%SNR=zeros(1,Nf);
if Ns>1
    EE_sort=zeros(Nel,Nf);
    VV=zeros(Nel,Nel,Nf);
    for If=1:Nf
        [VV0,EE]=eig(Kstot(:,:,If));
        EE=diag(EE);
        [EE_sort(:,If),Isort]=sort(EE,1,'descend');
        VV(:,:,If)=VV0(:,Isort);
        
    end
end

% if Ns<=1
%     SNRest=-SNRest;
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_Mode.
function pushbutton_Mode_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%Only called if MDAT exists...

tmp=ginput(1);
twant=tmp(1);
%fwant=tmp(2)*1000; %in Hz


tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%for Ichan=1:8
[x,t,Fs,tstart,junk,head]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,'all',handles);
if size(x,2)>1
    x=x';
end
%xall{Ichan}=x;
%end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');
for Ichan=1:head.Nchan
    
    %%Test 1
    [S,FF,TT,B] = spectrogram(x(:,Ichan)-mean(x(:,Ichan)),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
    [junk,Iwant]=min(abs(twant-TT));
    Bslice(:,Ichan)=B(:,Iwant);
    
    
    %Sign test
    index=round(twant*Fs)+(1:Nfft);
    Bslice2(:,Ichan)=fft((x(index,Ichan)-mean(x(index,Ichan))).*(hanning(Nfft)));
    
end
Bslice2=Bslice2((1:(Nfft/2+1)),:);
figure(1)
for I=1:2
    subplot(2,1,I)
    if I==1
        contourf(FF,head.rd,10*log10(abs(Bslice2')), 'linestyle','none');axis('ij')
        %caxis([80 120])
        title(sprintf('log magnitude Raw FFT %s: %s plus %6.2f sec',get(handles.text_filename,'String'),datestr(tdate_start),twant));
    else
        imagesc(FF,head.rd,abs(Bslice2'));  title('magnitude of FFT');
        %caxis([1e9 1e11])
    end
    
    xlim([0 1000])
    grid
    colorbar
    xlabel('Hz')
    ylabel('Channel');
end


figure(2)
for I=1:2
    subplot(2,1,I)
    if I==1
        contourf(FF,head.rd,10*log10(Bslice'), 'linestyle','none');axis('ij')
        caxis([80 120])
        title(sprintf('Log spectrogram %s: %s plus %6.2f sec',get(handles.text_filename,'String'),datestr(tdate_start),twant));
    else
        imagesc(FF,head.rd,(Bslice'));title('Spectrogram slice');
        %contourf(FF,head.rd,(Bslice'),1e10:1e10:10e10,'linestyle','none');
        %caxis([1e9 1e11])
    end
    
    xlim([0 1000])
    grid
    colorbar
    xlabel('Hz')
    ylabel('Channel');
end

keyboard
figure(1)
disp('Select a frequency:');
tmp=ginput(1);
freq=tmp(1);
[junk,Ibest]=min(abs(freq-FF));
mode=Bslice(Ibest,:);
figure(3)
subplot(3,1,1);
plot(head.rd,mode/max(abs(mode)),'-x');grid on;title(sprintf('Depth: %6.2f Frequency: %6.2f', head.D,freq))

subplot(3,1,2)
angg=180*(unwrap(angle(Bslice2(Ibest,:))))/pi;

plot(head.rd,angg,'-x');grid on;
title('unwrapped angle in degrees');

subplot(3,1,3)
angg=diff(unwrap(angle(Bslice2(Ibest,:))));
angg=[0 angg];
%%since angg=kL =2*pi*f*L/c in case of tilt
tilt_offset=angg*1490./(2*pi*freq);
plot(head.rd,cumsum(tilt_offset),'-x');grid on;
title('horizontal meter offset from top element');

yes=input('Save a mode example?','s');
if ~isempty(yes)
    modeno=input('Enter mode number:');
    save(deblank(sprintf('Mode%iExample%3.0fHz_%s',modeno,FF(Ibest),datestr(tdate_start+datenum(0,0,0,0,0,twant),30))),'mode','head','freq','angg','tilt_offset');
end

close(1:3)


end

function [B,wout]=conventional_beamforming(Ks,angles,freq,Lz,c,yesnorm)
B=zeros(length(freq),length(angles));
if size(Lz,2)>1
    Lz=Lz.';
end

if nargout==2
    wout=zeros(length(Lz),length(freq),length(angles));
end
winn=hanning(length(Lz));
for If=1:length(freq)
    for Iang=1:length(angles)
        % lambda=1500/freq(If);
        w=exp((-1i*2*pi*Lz*freq(If)/c)*sin(angles(Iang)*pi/180));
        w=w.*winn;
        
        w=w/norm(w);
        K=squeeze(Ks(:,:,If));
        if exist('yesnorm')
            K=K/norm(K);
        end
        B(If,Iang)=real(w'*K*w);
        if nargout==2
            wout(:,If,Iang)=w;
        end
    end
    
    
end

%plot(angles,10*log10(B(If,:)));
end


function B=MV_beamforming(Ks,angles,freq,Lz,c)
B=zeros(length(freq),length(angles));
if size(Lz,2)>1,
    Lz=Lz.';
end
for If=1:length(freq),
    for Iang=1:length(angles)
        % lambda=1500/freq(If);
        w=exp((i*2*pi*Lz*freq(If)/c)*sin(angles(Iang)*pi/180));
        
        %B(If,Iang)=real(w'*Ks{If}*w)/(norm(Ks{If})*norm(w).^2);
        w=w/norm(w);
        B(If,Iang)=1./real(w'*inv(squeeze(Ks(:,:,If)))*w);
        
    end
    
    
end

%plot(angles,10*log10(B(If,:)));
end


% --- Executes on button press in pushbutton_tilt.
function pushbutton_tilt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%for Ichan=1:8
[x,t,Fs,tstart,junk,head]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,'all',handles);

%x(1) is shallowest element...
if size(x,2)>1
    x=x';
end
%xall{Ichan}=x;
%end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

disp('Select min and max frequency:')
tmp=ginput(2);
frange=sort(tmp(:,2)*1000);
minfreq=frange(1);maxfreq=frange(2);
frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
[N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.1],Fs);
Bfilt = firpm(N,Fo,Ao,W);

for I=1:head.Nchan
    y(:,I)=filtfilt(Bfilt,1,x(:,I)-mean(x(:,I)));
end

offsets=zeros(head.Nchan,1);
offsets_cum=offsets;
maxval=zeros(head.Nchan,2);
Ioffset=1;
for I=1:(head.Nchan-Ioffset)
    [xc,lags]=xcov(y(:,I),y(:,I+Ioffset));  %lags are positive if array tilting towards signal
    [maxval(I,1),Ibest]=max(xc);
    offsets(I+1)=1.45*1000*lags(Ibest)/Fs;
    [maxval(I,2),junk]=max(-xc);
    
    figure(1);
    subplot(3,1,1)
    plot(lags/Fs,xc);xlim(4*[-0.05 0.05]);
    grid on
    title(sprintf('Element %i, leads previous by %10.6f m with %0.5g, max negative value is %i percent of max',I+1,offsets(I+1), maxval(I,1),round(100*maxval(I,2)/maxval(I,1))));
    
    subplot(3,1,2)
    [S,FF,TT,B] = spectrogram(y(:,I),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
    imagesc(TT,FF/1000,10*log10(B));%
    axis('xy');ylim([minfreq maxfreq]/1000);caxis([40 120]);
    
    subplot(3,1,3)
    [S,FF,TT,B] = spectrogram(y(:,I+Ioffset),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
    imagesc(TT,FF/1000,10*log10(B));%
    axis('xy');ylim([minfreq maxfreq]/1000);caxis([40 120]);
    
    [xc,lags]=xcov(y(:,1),y(:,I+Ioffset));  %lags are positive if array tilting towards signal
    [junk,Ibest]=max(xc);
    offsets_cum(I+1)=1.45*1000*lags(Ibest)/Fs;
    
    pause
end

% figure(1);
% subplot(3,1,1)
% plot(lags/Fs,xc);xlim(4*[-0.05 0.05]);
% grid on
% title(sprintf('Element %i, leads element 1 by %10.6f m ',head.Nchan,offset_max));

close;
figure(1)
plot(cumsum(offsets),head.rd,'-x',offsets_cum,head.rd,'-ro');
set(gca,'fontweight','bold','fontsize',14);
grid on;ylabel('receiver depth (m)');xlabel('tilt offset (m)');
legend('offset accumulated from adjacent elements','direct comparison to element 1','Location','best');
xlim([-5 5]);
hold on;axis('ij')
title(sprintf('Depth: %6.2f m Times: %s, freq range: %i to %i Hz', head.D,datestr(tdate_start,30),round(minfreq),round(maxfreq)));

yes=menu('Save a tilt example?','Yes','No');
if yes==1
    figure(1)
    orient landscape
    print(1,'-djpeg',sprintf('TiltExample_%s_%ito%iHz',datestr(tdate_start,30),round(minfreq),round(maxfreq)));
end
%close(1);
%%%Test to understand xcorr lags
%x1=[0 0 0 1 0 0 0 0 0 0];
%x2=[0 0 0 0 0 0 1 0 0 0];  %x2 lags x1 (signal arrives later on x2)
%[cc,lags]=xcov(x1,x2); %lags max at -4

%Thus if first vector leads second, lags are negative
%  If first vector lags second, lags are positive

end

function R=derive_reflection_coefficient2(Kstot,angles,freq,rd,c)
%function [B,wout]=conventional_beamforming(Ks,angles,freq,Lz,c,yesnorm)

tilt=input('Enter tilt estimate in deg (0):');
if isempty(tilt)
    tilt=0;
end

% angles=angles0+tilt;
%
% Igood=find(angles<=90&angles>=-90);
% angles=angles(Igood);
%
% Nup=Iang-1;
% Ndown=length(Iang:length(angles));
% if Ndown>=Nup
%    angles=angles(1:(Iang+Nup));
% else
%     angles=angles((Nup-Ndown+2):end);
%end
Lz=rd-rd(1);

B=conventional_beamforming(Kstot,angles,freq,Lz,c,'yesnorm');
%B=MV_beamforming(Kstot,angles,freq,Lz,c,'yesnorm');
B=real(B).';  % B now [angles freq]
figure
imagesc(freq,angles,10*log10(abs(B)))
caxis([-20 0]);
axis('xy')
set(gca,'fontweight','bold','fontsize',14);
xlabel('Frequency (Hz)');
ylabel('Angle (positive is towards surface');
grid on
%title(sprintf('Theoretical conventional beamforming response, even spacing, source angle %6.2f',angle))
colorbar('fontweight','bold','fontsize',14)


for It=1:length(tilt)
    [junk,Iang]=min(abs(angles-tilt(It)));
    
    %R is bottom/surface, or negative over positive
    %upp=flipud(B(1:(Iang-1),:));  %Negative angles point towards bottom
    %downn=B((Iang+1):end,:);
    %Nmin=min([size(downn,1) size(upp,1)]);
    %R=downn(1:Nmin,:)./upp(1:Nmin,:);
    
    downn=flipud(B(1:(Iang-1),:));  %Negative angles point towards bottom
    upp=B((Iang+1):end,:); %Positive angles point towards surface
    Nmin=min([size(downn,1) size(upp,1)]);
    R=downn(1:Nmin,:)./upp(1:Nmin,:);
    
    %%Note that R is the power reflection coefficient (Harrison and Simmons, Eq. (2)), so 10*logR please..
    R=-10*log10(abs(R))';
    angle_graze=angles(Iang+(1:Nmin))-tilt(It);
    figure
    
    subplot(3,1,1)
    imagesc(10*log10(abs(upp')))
    subplot(3,1,2)
    imagesc(10*log10(abs(downn')))
    subplot(3,1,3)
    imagesc(angle_graze,freq,R);
    caxis([0 15]);
    set(gca,'fontweight','bold','fontsize',14);
    grid on;colorbar;axis('xy')
    xlabel('Grazing angle (deg)');
    ylabel('Frequency (Hz)')
    title(sprintf('reflection loss (dB), tilt %6.2f',tilt(It)))
    
end

yes=input('Save?');
if ~isempty(yes)
    save(sprintf('Refl_%s.mat',num2str(mean(tilt))),'angle_graze','freq','R','tilt','rd','c','Kstot','angles','Lz');
end
% disp('Select a frequency:');
% tmp=ginput(1);
% [junk,Iwant]=min(abs(tmp(2)-freq));
% figure
% plot(angle_graze,R(Iwant,:));
% keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function pushbutton_modalfiltering_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_modalfiltering.
function pushbutton_modalfiltering_Callback(hObject, eventdata, handles)
%% Estimate distance and range by isolating modes on a vertical arrray...
% handles    structure with handles and user data (see GUIDATA)

figss=get(0,'child');
Iclose=find(figss-floor(figss)==0);
if ~isempty(figss(Iclose))
    close(figss(Iclose));
end
tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%for Ichan=1:8
[data.x,t,Fs,tstart,junk,head]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,'all',handles);

%xall{Ichan}=x;
%end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

disp('Select min and max frequency:')
tmp=ginput(2);
frange=sort(tmp(:,2)*1000);
minfreq=frange(1);maxfreq=frange(2);

frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
[N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.05],Fs);
Bfilt = firpm(N,Fo,Ao,W);

Nchan=size(data.x,1);
y=zeros(size(data.x));
for I=1:Nchan
    try
        y(I,:)=filtfilt(Bfilt,1,data.x(I,:)-mean(data.x(I,:)));
        %[B(:,:,I),FF,TT] = specgram(y(I,:),Nfft,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
    catch
        keyboard
    end
end

%%Convert time series into one long FFT
Np=size(y,2);
Nfft=2^ceil(log10(Np)/log10(2));
%B=permute(B,[2 3 1]);  %time, element, freq
B=fft([zeros((Nfft-Np)/2,Nchan); y'.*(hanning(Np)*ones(1,Nchan)); zeros((Nfft-Np)/2,Nchan)],Nfft);
B=B.'; %(Nel, freq);



[MFP_scenario,range_str,default_tilt]=load_MFP_scenario;

if length(MFP_scenario)>1  %Multiple runs..
    prompt1={'range guesses (m)', 'maxdelay for mode 1:2 range estimation (s)[-1 to automatically select]:', ...
        'maxdelay for mode 1:N range estimation (s): [-1 to automatically select]','plot intermediate images?', ...
        'Nfft for plotting:','depths','tilt','max modes'};
    for J=1:length(default_tilt)
        Islash=max(findstr(MFP_scenario{J},'/'))+1;
        disp(sprintf('tilt: %s File: %s ',default_tilt{J},MFP_scenario{J}(Islash:end)));
    end
    dlgTitle1='Parameters for modal beamforming...';
    def1={range_str{1},'0.1','0.025','0','2048','1:1:55',default_tilt{1},'3'};
    answer=inputdlg(prompt1,dlgTitle1,1,def1);
    range_guess=str2num(answer{1});
    maxdelay(1)=str2num(answer{2});
    maxdelay(2)=maxdelay(1);
    maxdelay(3)=str2num(answer{3});
    plot_chc=str2num(answer{4});
    Nfft2=str2num(answer{5});
    depths=str2num(answer{6});
    tilt_offset=str2num(answer{7});
    Maxmodes=str2num(answer{8});
    
    
end


for Imodel=1:length(MFP_scenario)
    
    if length(MFP_scenario)==1
        
        prompt1={'Model file..','tilt offsets (deg)','range guesses (m)', ...
            'maxdelay for mode 1:2 range estimation (s):', ...
            'maxdelay for mode 1:N range estimation (s):', ...
            'plot intermediate images?', ...
            'Nfft for plotting:','depths','max modes'};
        
        dlgTitle1='Parameters for modal beamforming...';
        
        def1={MFP_scenario{Imodel}, default_tilt{Imodel},range_str{Imodel},'0.1','0.025','0','2048','1:0.5:55','3'};
        
        answer=inputdlg(prompt1,dlgTitle1,1,def1);
        model=load(answer{1});
        disp(sprintf('Loaded %s',answer{1}));
        tilt_offset=str2num(answer{2});
        range_guess=str2num(answer{3});
        maxdelay(1)=str2num(answer{4});
        maxdelay(2)=maxdelay(1);
        maxdelay(3)=str2num(answer{5});
        plot_chc=str2num(answer{6});
        Nfft2=str2num(answer{7});
        depths=str2num(answer{8});
        Maxmodes=str2num(answer{9});
    else
        
        model=load(MFP_scenario{Imodel});
        disp(sprintf('Loaded %s',MFP_scenario{Imodel}));
        %tilt_offset=str2num(default_tilt{Imodel});
        
    end
    maxdelay(4:Maxmodes)=maxdelay(3);
    
    disp('Save signal sample if you wish...');
    keyboard
    %following operation is now done in load_data...
    %data.x=flipud(data.x); %Shallowest data now first...
    
    %%Remove bad channels, if any
    %[junk,Igood]=intersect(round(model.rd),round(head.rd));
    %model.rd=model.rd(Igood);
    %model.U=model.U(Igood,:,:);
    
    Igood=zeros(1,length(head.rd));
    for I=1:length(head.rd)
        [junk,Igood(I)]=min(abs(model.rd-head.rd(I)));
    end
    
    %%depths are desired modeled source depths.
    Igood_source=zeros(1,length(depths));
    for I=1:length(depths)
        [junk,Igood_source(I)]=min(abs(model.rd-depths(I)));
    end
    
    model.sd=model.rd(Igood_source);
    model.Urd=model.U(Igood,:,:);
    
    model.U=model.U(Igood_source,:,:);
    model.rd=model.rd(Igood);
    
    
    %%%Match frequency bins with frequencies computed in model
    If_min=floor(minfreq*Nfft/Fs);
    If_max=ceil(maxfreq*Nfft/Fs);
    index=If_min:If_max;
    FF=(Fs/Nfft)*index;
    Iwant=zeros(length(index),1);  %Iwant will be indicies in model frequencies that match FFT bins in index
    for If=1:length(index)
        [junk,Iwant(If)]=min(abs(model.freq-FF(If)));
    end
    kk=2*pi*model.freq(Iwant)./1495;
    
    %%Compute expected group delays between mode pairs
    if any(maxdelay<0)
        for Im=1:Maxmodes
            Igg1=find(model.cgroup(Iwant,Im))>0;
            Igg1=Iwant(Igg1);
            for Im2=(Im+1):Maxmodes
                Igg2=find(model.cgroup(Iwant,Im2))>0;
                Igg2=Iwant(Igg2);
                Sij(Im,Im2)=median(abs(1./model.cgroup(Igg1,Im)-1./model.cgroup(Igg2,Im2)));
            end
            
        end
    end
    
    
    windowpwr=sum(hanning(Nfft).^2)/Nfft;
    normm=1./sqrt(Nfft*Fs*windowpwr);
    %normm=1;
    
    %Make initial normalization of modes such that numerical integration approximations analytical result
    depth_weights=sqrt([model.rd(1) diff(model.rd) ]');  %Set dz
    
    
    %%Gosh, this has been the incorrect formula for a linear tilt for a long time...
    %rd_cum=model.rd./cumsum(model.rd);
    
    rd_cum=-(model.rd-min(model.rd))./(max(model.rd)-min(model.rd));
    
    %%Depth estimator setups...
    magg{Imodel}=zeros(length(tilt_offset),length(range_guess),Maxmodes);
    talign{Imodel}=magg{Imodel};
    talign_allModeOne{Imodel}=zeros(length(tilt_offset),length(range_guess));
    talign_all{Imodel}=zeros(length(tilt_offset),length(range_guess));
    
    %%Cross-correlation has different FFT size than orignal processing, so adjust relevent variables...
    %Nxc=2*round(Fs*maxdelay(1))+1;
    Nxc=2*round(Fs*max(maxdelay))+1;
    Nfft_xc=2^ceil(log10(Nxc)/log10(2));
    If_min1=floor(minfreq*Nfft_xc/Fs);
    If_max1=ceil(maxfreq*Nfft_xc/Fs);
    index2=If_min1:If_max1;
    Iwant2=zeros(length(index2),1);  %Iwant will be indicies in model frequencies that match FFT bins in index
    for If=1:length(index2)
        [junk,Iwant2(If)]=min(abs(model.freq-Fs*index2(If)/Nfft_xc));
    end
    
    
    
    for Itilt=1:length(tilt_offset)  %Each tilt offset will get a separate figure...
        disp(sprintf('Tilt: %6.2f',tilt_offset(Itilt)));
        
        dvec=zeros(length(range_guess),Maxmodes,length(index2));
        for Ir=1:length(range_guess)
            if rem(Ir,10)==0, disp(sprintf('%6.2f percent ranges done',100*Ir/length(range_guess)));end
            Bout=zeros(Nfft,Maxmodes);
            Bout_corrected=Bout;
            tilt_matrix=exp(-1i*tilt_offset(Itilt)*rd_cum'*kk);
            for Im=1:Maxmodes
                kg=real(model.kr{1}(Iwant,Im))-(2*pi*model.freq(Iwant)'./max(model.cgroup(Iwant,1)));
                %Last term is maximum group velocity of group 1 across all frequencies.
                % As mode 1 is fastest mode, this is fastest velocity across all
                % frequencies.
                
                %kg=ones(size(kg));
                
                %%From Buck 1998 paper on modal extraction, weight SMS filter by 1./U(z).^2 to ensure U'*U
                %%diagonal is unbiased (all ones)
                
                %U=squeeze(model.Urd(:,Iwant,Im)).*(depth_weights*ones(1,length(index)));  %(rd,freq)
                U=squeeze(model.Urd(:,Iwant,Im));
                U=U./(ones(Nchan,1)*(sum(U.^2)));
                Bout(index,Im)=trapz(B(:,index).*U.*tilt_matrix);
                
                %%%Alternate PI inversion....
                
                kg_dispersed=exp(1i*range_guess(Ir)*ones(length(model.rd),1)*kg.');
                Bout_corrected(index,Im)=trapz(B(:,index).*U.*tilt_matrix.*kg_dispersed);
                
                
            end  %Im
            
            yout=2*real(ifft(Bout,Nfft));
            yout_corrected=2*real(ifft(Bout_corrected,Nfft));
            
            %tilt=asin(tilt_offset(Itilt)./(max(model.rd)-min(model.rd)))*180/pi;
            
            if plot_chc==1
                tilt_check(Itilt)=plot_mode_filtering;
                if length(tilt_offset)>1
                    figure(1)
                    pause(0.1);
                    disp('capturing frame');
                    MM(Itilt)=getframe(gcf);
                    
                end
            end
            
            for Im=2:Maxmodes
                %maxdelay(Im)=2*range_guess(Ir).*(1./model.cgroup
                
                if any(maxdelay<0)
                    maxdelayt(1,Im)=Sij(1,Im)*range_guess(Ir);
                else
                    maxdelayt(1,Im)=maxdelay(Im);
                end
                [xc,lags]=xcov(yout_corrected(:,1),yout_corrected(:,Im),round(Fs*maxdelayt(1,Im)));  %lags are positive if array tilting towards signal
                [magg{Imodel}(Itilt,Ir,Im),Ibest]=max(abs(xc));
                magg{Imodel}(Itilt,Ir,Im)=sign(xc(Ibest))*magg{Imodel}(Itilt,Ir,Im);
                talign{Imodel}(Itilt,Ir,Im)=lags(Ibest)/Fs;  %talign estimates time difference between modes 1 and Im
                %%%Add square of time offsets between phase-shifted modal arrivals, weighted by peak of
                %%%correlation (which should be proportional to modal excitation).
                %%% To avoid situations where a low peak reduces error, normalize weighting factor by peak
                %%% power (so relative weighting between modes preserved)
                
                if Im==2 %first round
                    magg_norm=abs(magg{Imodel}(Itilt,Ir,Im));  %Mode 12 interference weight
                    
                end
                
                wght{Imodel}(Itilt,Ir,Im)=abs(magg{Imodel}(Itilt,Ir,Im))./magg_norm;
                talign_allModeOne{Imodel}(Itilt,Ir)=talign_allModeOne{Imodel}(Itilt,Ir)+wght{Imodel}(Itilt,Ir,Im).*talign{Imodel}(Itilt,Ir,Im).^2;
                talign_all{Imodel}(Itilt,Ir)=talign_all{Imodel}(Itilt,Ir)+wght{Imodel}(Itilt,Ir,Im).*talign{Imodel}(Itilt,Ir,Im).^2;
                
                if Im>1
                    for Im2=(Im+1):Maxmodes
                        %fprintf('Computing xcorr of %i and %i\n',Im,Im2);
                        if any(maxdelay<0)
                            maxdelayt(Im,Im2)=Sij(Im,Im2)*range_guess(Ir);
                            
                        else
                            maxdelayt(Im,Im2)=maxdelay(Im);
                        end
                        
                        [xc2,lags2]=xcov(yout_corrected(:,Im2),yout_corrected(:,Im),round(Fs*maxdelayt(Im,Im2)));  %lags are positive if array tilting towards signal
                        [magg2,Ibest2]=max(abs(xc2));
                        wght2{Imodel}(Itilt,Ir,Im,Im2)=abs(sign(xc2(Ibest2))*magg2/magg_norm);
                        talign2=lags2(Ibest2)/Fs;  %talign estimates time difference between modes 1 and Im
                        talign_all{Imodel}(Itilt,Ir)=talign_all{Imodel}(Itilt,Ir)+wght2{Imodel}(Itilt,Ir,Im,Im2).*talign2.^2;
                        
                    end
                end
                
                %if min(maxdelay)~=maxdelay(Im)  %When estimating depth use full xc
                [xc,lags]=xcov(yout_corrected(:,1),yout_corrected(:,Im),round(Fs*max(maxdelay)));  %lags are positive if array tilting towards signal
                % end
                %%Try to keep the function as even as possible (as small a phase as possible)
                %%Try to multiply Bout directly?  But need to make sure no sidelobes are included...
                zpad=floor(0.5*(Nfft_xc-length(xc)));
                xc_pad=fftshift([zeros(zpad,1); hanning(length(xc)).*xc ;zeros(zpad,1)]);
                xc_FFT=fft(xc_pad,Nfft_xc);
                dvec(Ir,Im,:)=xc_FFT(index2);
                
                if plot_chc==1
                    plot_xcorr;
                end
            end  %Im
            
            
        end %Ir
        
        
        
        
        %%%Estimate best range from each tilt computation...
        compute_depth_estimate;
        
    end %Itilt
    %%Plot summary of range estimation time lags...
    plot_range_estimates;
    
    if plot_chc==1
        if exist('MM');
            movie2avi(MM,'tiltTest');
        end
        keyboard
        for I=1:12
            if any(get(0,'child')==I)
                figure(I)
                orient tall
                print(I,'-djpeg',sprintf('ModalBeamform%i',I));
            end
        end
    end
end %Imodel

%%%Inner function for plotting range estimates
    function plot_range_estimates
        %%%Plot time alignment vs hypothesized source range and look for zero crossing
        persistent best_range fitt best_range2 fitt2
        
        if Imodel==1
            best_range=[]; fitt=[]; best_range2=[]; fitt2=[];
        end
        plotstr='oxs+*';
        plotstr_tilt='kbgrckbgrckbgrckbgrc';
        
        if length(range_guess)>1  %Plot only once...
            for Itilt2=1:length(tilt_offset)
                figure(10+Itilt2)
                
                %%Plot individual mode pair combinations
                for Imm=2:Maxmodes
                    subplot(Maxmodes+1,1,Imm-1);
                    %amplitudes=1000*abs(magg{Imodel}(Itilt2,:,Imm));
                    %[ax,h1,h2]=plotyy(range_guess/1000,1000*(squeeze(talign{Imodel}(Itilt2,:,Imm))),range_guess/1000,amplitudes);
                    
                    amplitudes=squeeze(wght{Imodel}(Itilt2,:,Imm));
                    [ax,h1,h2]=plotyy(range_guess/1000,1000*(squeeze(talign{Imodel}(Itilt2,:,Imm))),range_guess/1000,amplitudes);
                    
                    
                    set(h1,'markeredgecolor',[plotstr_tilt(Itilt2) ])
                    set(h1,'marker', plotstr(Imodel));set(ax(2),'xtick',[])
                    set(gca,'fontsize',14,'fontweight','bold');
                    set(ax(2),'fontsize',14,'fontweight','bold');
                    ylabel(ax(2),'Amplitude');
                    grid on
                    xlabel('Range guessed (km)');ylabel('lag (ms)');
                    
                    title(sprintf('legend: %s, Cross-correlation between modes 1 and %i, Tilt: %6.2f Freq range: %6.1f-%6.1f Hz', ...
                        plotstr,Imm, tilt_offset(Itilt2),minfreq,maxfreq));
                    hold on
                    line(range_guess([1 end])/1000,maxdelayt(1,Imm)*[1 1]*1000);
                    line(range_guess([1 end])/1000,-maxdelayt(1,Imm)*[1 1]*1000);
                    %if Imm==2, ylim([-25 25]);end
                    ylim(maxdelayt(1,Imm)*[-1 1]*1000);
                    ylimm=ylim;
                    set(ax(1),'ytick',linspace(ylimm(1),ylimm(2),5));
                end
                
                
                %%Plot MSE error vs. range estimate, mode 1 combos only%%%%%%%%
                [fitt(Itilt2,Imodel),Ir_chc]=min(talign_allModeOne{Imodel}(Itilt2,:));
                best_range(Itilt2,Imodel)=range_guess(Ir_chc)/1000;
                
                subplot(Maxmodes+1,1,Maxmodes);
                plot(range_guess/1000,1000*sqrt(squeeze(talign_allModeOne{Imodel}(Itilt2,:))),[plotstr_tilt(Itilt2) '-' plotstr(Imodel)]);
                set(gca,'fontsize',14,'fontweight','bold');
                grid on;xlabel('Range guessed (km)');ylabel('MSE (ms)');
                title(sprintf('legend: %s, MSE, all mode 1 combos: %6.3g ms at %6.2f km, Tilt: %6.2f Freq range: %6.1f-%6.1f Hz',plotstr, ...
                    1000*sqrt(fitt(Itilt2,Imodel)), range_guess(Ir_chc)/1000, tilt_offset(Itilt2),minfreq,maxfreq));
                hold on
                %line(range_guess([1 end])/1000,maxdelay(Imm)*[1 1]*1000);
                %line(range_guess([1 end])/1000,-maxdelay(Imm)*[1 1]*1000);
                plot(best_range(Itilt2,Imodel),1000*sqrt(fitt(Itilt2,Imodel)),[plotstr_tilt(Itilt2) 's' ],'markerfacecolor',plotstr_tilt(Itilt2));
                %ylim([0 50]);
                
                %%%Plot MSE error vs. range estimate for all mode combinations...
                subplot(Maxmodes+1,1,Maxmodes+1);
                plot(range_guess/1000,1000*sqrt(squeeze(talign_all{Imodel}(Itilt2,:))),[plotstr_tilt(Itilt2) '-' plotstr(Imodel)]);
                set(gca,'fontsize',14,'fontweight','bold');
                grid on
                xlabel('Range guessed (km)');ylabel('MSE (ms)');
                title(sprintf('legend: %s, MSE, all modes, Tilt: %6.2f Freq range: %6.1f-%6.1f Hz',plotstr, tilt_offset(Itilt2),minfreq,maxfreq));
                hold on
                [fitt2(Itilt2,Imodel),Ir_chc2]=min(talign_all{Imodel}(Itilt2,:));
                best_range2(Itilt2,Imodel)=range_guess(Ir_chc2)/1000;
                plot(best_range2(Itilt2,Imodel),1000*sqrt(fitt2(Itilt2,Imodel)),[plotstr_tilt(Itilt2) 's' ],'markerfacecolor',plotstr_tilt(Itilt2));
                
                %line(range_guess([1 end])/1000,maxdelay(Imm)*[1 1]*1000);
                %line(range_guess([1 end])/1000,-maxdelay(Imm)*[1 1]*1000);
                
            end  %Itilt2
            
        end  %if range_guess length is 1
        
        
        if Imodel==length(MFP_scenario)
            yes=menu('Save and Print result?','Yes','No');
            if yes==1
                for Itilt2=1:length(tilt_offset)
                    savename_base=sprintf('%s_%ikm_tilt%4.2f_%ito%iHz_%imodels',datestr(tdate_start,30),round(range_guess(1)/1000), ...
                        tilt_offset(Itilt2),round(minfreq),round(maxfreq),Imodel);
                    
                    figure(10+Itilt2)
                    orient landscape
                    print(gcf,'-djpeg',sprintf('RangeEstimate_%s.jpg',savename_base));
                    
                    %%%%%Plot best range vs. fit%%%%%%
                    figure(1)
                    
                    subplot(2,1,1)
                    for Imodel2=1:length(MFP_scenario)
                        %subplot(length(MFP_scenario),1,Imodel2);
                        plot(best_range(Itilt2,Imodel2),1000*sqrt(fitt(Itilt2,Imodel2)),[plotstr_tilt(Itilt2) plotstr(Imodel2)]);
                        hold on
                    end
                    xlabel('Best range (km)');
                    ylabel('MSE (msec)');
                    title(sprintf('Mode 1 contribution only: Color %s is tilt, %s is inversion model',plotstr_tilt,plotstr));grid on
                    
                    subplot(2,1,2)
                    for Imodel2=1:length(MFP_scenario)
                        %subplot(length(MFP_scenario),1,Imodel2);
                        plot(best_range2(Itilt2,Imodel2),1000*sqrt(fitt2(Itilt2,Imodel2)),[plotstr_tilt(Itilt2) plotstr(Imodel2)]);
                        hold on
                    end
                    xlabel('Best range (km)');
                    ylabel('MSE (msec)');
                    title(sprintf('All mode contributions: Color %s is tilt, %s is inversion model',plotstr_tilt,plotstr));grid on
                    
                end
                
                %%%%%Plot min error vs array tilt...
                figure(2)
                
                subplot(2,1,1)
                for Imodel2=1:length(MFP_scenario)
                    plot(tilt_offset,1000*sqrt(fitt(:,Imodel2)),['k-' plotstr(Imodel2)]);hold on;grid
                    
                end
                xlabel('tilt offset (m)');
                ylabel('mse error (msec)');
                title(sprintf('Mode 1 only: %s is inversion model',plotstr));grid on
                
                subplot(2,1,2)
                for Imodel2=1:length(MFP_scenario)
                    plot(tilt_offset,1000*sqrt(fitt2(:,Imodel2)),['k-' plotstr(Imodel2)]);hold on;grid
                    
                end
                xlabel('tilt offset (m)');
                ylabel('mse error (msec)');
                title(sprintf('All modes: %s is inversion model',plotstr));grid on
                
                
                figure(1)
                orient landscape
                legend('full model','Pekeris model','Averaged Pekeris')
                print(gcf,'-djpeg',sprintf('RangeError1_%s.jpg',savename_base));
                
                figure(2)
                orient landscape
                legend('full model','Pekeris model','Averaged Pekeris')
                print(gcf,'-djpeg',sprintf('RangeError2_%s.jpg',savename_base));
                
                
                save(sprintf('ModalEstimateRange_%s.mat',savename_base), ...
                    'range_guess','tilt_offset','talign','talign_all','talign_allModeOne','MFP_scenario',...
                    'minfreq','maxfreq','maxdelay','plotstr','plotstr_tilt','Maxmodes', ...
                    'fitt','best_range','fitt2','best_range2','magg','wght');
                
                
            end
        end
    end

%%%%%Inner function for depth estimation
    function compute_depth_estimate
        persistent depth_chc
        depth_chc=1;
        
        if Itilt==1&&Imodel==1
            %depth_chc=menu('Select depth estimation choice:','No depth estimation','Depth estimate for all ranges','Depth estimate at best mode one only range estimate');
        end
        %depth_chc=1;
        
        if depth_chc==2
            Ir_chc=1:length(range_guess);
        else
            Ir_chc=1;
            fitt=Inf;
            if length(range_guess)>1
                %dt=diff(sign(squeeze(talign{Imodel}(Itilt,:,2))));dt=[0 dt];
                [fitt,Ir_chc]=min(talign_allModeOne{Imodel}(Itilt,:));
                
            end
            if isempty(Ir_chc)
                Ir_chc=-1;
                fprintf('No best range for Model %i..\n',Imodel);
            else
                fprintf('Best range for Model %i: %6.2f km, error %6.2f msec\n',Imodel,range_guess(Ir_chc)/1000,sqrt(fitt));
            end
            if depth_chc==1
                return
            end
        end
        
        %%%Estimate source depth by relative amplitudes of cross correlations with mode 1 extraction...
        mvec=model.U(:,Iwant2,1:Maxmodes);  %All real values
        sqkr=sqrt(real(model.kr{1}(Iwant2,1:Maxmodes)));
        for Is=1:length(model.sd)
            mvec(Is,:,:)=(squeeze(model.U(Is,Iwant2,1))'*ones(1,Maxmodes)).*squeeze(mvec(Is,:,:))./sqkr;  %U1(zs)*Um(zs)/sqrt(km)
        end
        mvec(isnan(mvec))=0;
        amb{Imodel}=zeros(length(model.sd),length(Iwant2));
        amb_abs{Imodel}=amb{Imodel};
        amb_nowght{Imodel}=amb{Imodel};
        
        for Ir_best=Ir_chc
            
            dd=squeeze(dvec(Ir_best,:,:));
            for I=1:Maxmodes
                freq_wght(I,:)=abs(dd(I,:))./max(abs(dd(I,:)));
            end
            freq_wght=mean(freq_wght);  %Frequency weight for all modes
            for J=1:length(Iwant2)
                for I=1:length(model.sd)
                    mm=real(squeeze(mvec(I,J,:)));
                    Iexist=find(mm~=0);
                    mm(Iexist)=mm(Iexist)/norm(mm(Iexist));
                    %Don't penalize absent modes...
                    dd(Iexist,J)=dd(Iexist,J)/norm(dd(Iexist,J));
                    if Ir_best>0
                        amb{Imodel}(I,J)=freq_wght(J).*abs(dd(Iexist,J)'*mm(Iexist)).^2;  %FT of xcorr data, unnormalized...,
                        amb_abs{Imodel}(I,J)=freq_wght(J).*abs(dd(Iexist,J))'*abs(mm(Iexist));  %FT of xcorr data, unnormalized...,
                        
                        amb_nowght{Imodel}(I,J)=abs(dd(Iexist,J)'*mm(Iexist)).^2; %No frequency weighting
                        %amb2{Imodel}(I,J)=(mm'*dvec1).^2;   %Real/real comparison...
                    else
                        amb{Imodel}(I,J)=NaN;
                        amb_abs{Imodel}=NaN;
                        amb_nowght{Imodel}(I,J)=NaN;
                        %amb2{Imodel}(I,J)=NaN;
                    end
                end
            end
            
            amb{Imodel}=amb{Imodel}/max(max(amb{Imodel}));
            amb_abs{Imodel}=amb_abs{Imodel}/max(max(amb_abs{Imodel}));
            amb_nowght{Imodel}=amb_nowght{Imodel}/max(max(amb_nowght{Imodel}));
            
            %%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            savename_base=plot_depth_estimates;
        end  %Ir_best
        
        
        if Itilt==1&&Imodel==1
            
            yes=menu('Save result?','Yes','No');
            if yes==1
                
                
                save(sprintf('ModalEstimateDepthTilt%i_%s.mat',tilt_offset(Itilt),savename_base), ...
                    'tdate_start','range_guess','Ir_best','Iwant','tilt_offset','Itilt','amb','amb_abs','amb_nowght', ...
                    'model','Imodel','MFP_scenario','minfreq','maxfreq','maxdelay','plotstr','plotstr_tilt','Maxmodes','fitt');
                
            end
        end
        
        %%%%%%%%%Inner function plot_range_estimates%%%%%%%%%%%%
        
        %%%Inner function plot_depth_estimates%%%%%%
        function savename_base=plot_depth_estimates
            %%Plot graphs of depth estimators
            
            savename_base=sprintf('%s_%ikm_tilt%4.2f_%ito%iHz_%imodels',datestr(tdate_start,30),round(range_guess(min(Ir_best))/1000), ...
                tilt_offset(Itilt),round(minfreq),round(maxfreq),Imodel);
            
            plotstr='oxs+*';
            plotstr_tilt='kbgrc';
            
            %%%Plot 2-D depth estimator
            figure(30+10*Imodel+Itilt+Ir_best);
            
            subplot(2,1,1)
            imagesc(model.freq(Iwant),model.sd,10*log10(amb{Imodel}));
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');
            ylabel('Source depth (m)');
            title(sprintf('Depth estimate of call starting at %s, using FT of xc function',datestr(tdate_start)));
            plot_letter_label('a)');
            
            subplot(2,1,2)
            imagesc(model.freq(Iwant),model.sd,10*log10(amb_nowght{Imodel}));
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');
            ylabel('Source depth (m)');
            title('No frequency weighting');
            orient tall
            print(gcf,'-djpeg',sprintf('DepthEstimate2D_%s.jpg',savename_base));
            
            
            %Plot frequency-averaged depth estimator%
            
            figure(100+Itilt+Ir_best);
            plot1=10*log10(sum(amb{Imodel}'))/length(Iwant);
            plot_norm=10*log10(sum(amb_nowght{Imodel}'))/length(Iwant);
            plot_abs=10*log10(sum(amb_abs{Imodel}'))/length(Iwant);
            plot(plot1/max(abs(plot1)),model.sd,[plotstr_tilt(Itilt) '-' plotstr(Imodel)]);hold on
            %plot(plot_norm/max(abs(plot_norm)),model.sd,[plotstr_tilt(Itilt) '-' plotstr(Imodel+1)]);
            plot(plot_abs/max(abs(plot_abs)),model.sd,[plotstr_tilt(Itilt) '-' plotstr(Imodel+1)]);
            
            set(gca,'fontweight','bold','fontweight','bold');grid on
            ylabel('Source depth (m)');xlabel('Normalized depth response');
            try
                title(sprintf('legend: %s, Depth estimate, Range: %6.2f km, Tilt: %6.2f m, Freq range: %6.1f-%6.1f Hz', ...
                    plotstr_tilt,range_guess(Ir_best)/1000,tilt_offset(Itilt),minfreq,maxfreq));
            end
            axis('ij')
            hold on
            legend('frequency weighted','abs value');
            print(gcf,'-djpeg',sprintf('DepthEstimate_%s.jpg',savename_base));
            
            
            
            
            
        end  %Inner function plot_depth_estimates
        
    end  %inner function compute_depth_estimates
%%%%%%%%%%%%%%%%%


%%Inner function for plotting
    function tltchk=plot_mode_filtering
        persistent timmes_signal timmes_noise pwrr yes
        
        if Itilt==1
            %yes=menu('What intermediate data do you want?','Filtered spectrograms','SNR estimate');
            yes=1;
            pwrr=[];
        end
        if yes==1
            
            y1=yout(:,1);
            y2=yout(:,2);
            tltchk=y1'*y2./(norm(y1).*norm(y2));
            fprintf('Xcor for bottom sensor, tilt %6.2f: %6.4f, \n',  ...
                tilt_offset(Itilt),tltchk);
            
            figure(1)
            
            for Imm=1:3
                [Bbeam,FF1,TT] = specgram(yout(:,Imm),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
                figure(1)
                set(gcf,'units','norm','pos',[ 0.1589    0.3345    0.2216    0.6102]);
                
                subplot(3,1,Imm);
                imagesc(TT,FF1,20*log10(abs(Bbeam).*normm));%
                orient landscape
                grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
                axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);
                caxis([-150 -70])
                if Imm==1
                    %title(sprintf('Mode %i, modal beamforming only: vertical tilt: %6.2f deg',Imm, tilt));
                end
                ylabel('Hz');
                xlabel('Time (s)');
                if Imm<3
                    set(gca,'xticklabel',[]);
                    xlabel('');
                end
                set(gca,'units','norm');
                poss=get(gca,'pos');
                poss(4)=.27;
                %set(gca,'pos',poss,'xtick',0:0.25:3);
                
                text(0.1,0.8,sprintf('Mode %i',Imm),'units','norm','color','w','fontsize',14,'fontweight','bold');
                title(sprintf('Modal Filtering, Tilt: %6.2f deg, correlation: %6.4f',tilt_offset(Itilt),tltchk));
            end
            hh=colorbar('south','fontweight','bold','fontsize',14);
            set(hh,'pos',[0.1369    0.0239    0.7462    0.0234]);
            
            
            
            for Imm=1:3
                [Bbeam1,FF1,TT] = specgram(yout_corrected(:,Imm),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
                figure(2)
                subplot(3,1,Imm);
                imagesc(TT,FF1,20*log10(abs(Bbeam1).*normm));%
                orient landscape
                grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
                ylabel('Hz'); xlabel('Time (s)');
                
                axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);grid on
                caxis([-150 -70])
                if Imm==1
                    %title(sprintf('Mode %i, modal beamforming and dispersion correction: vertical tilt: %6.2f deg',Imm, tilt));
                end
                if Imm<3
                    set(gca,'xticklabel',[]);
                    xlabel('');
                    
                end
                set(gca,'units','norm');
                poss=get(gca,'pos');
                poss(4)=.27;
                
                %set(gca,'pos',poss,'xtick',0:0.25:3);
                text(0.25,110,sprintf('Mode %i',Imm),'color','w','fontsize',14,'fontweight','bold');
                title(sprintf('Phase Correction, Range: %6.2f km, Tilt: %6.2f deg',range_guess(Ir)/1000,tilt_offset(Itilt)));
            end
            colorbar('south','fontweight','bold','fontsize',14);
            
            
            
        else
            %%Estimate SNR of signal...
            
            [tmp,lags]=xcov(yout(:,1),y(1,:)');
            [junk,Imax]=max(abs(tmp));
            Imax=lags(Imax);
            
            figure(10);
            
            [Bbeam0,FF1,TT] = specgram(y(1,:),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
            subplot(2,1,1);
            
            imagesc(TT,FF1,20*log10(abs(Bbeam0).*normm));%
            orient landscape
            grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
            axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);
            xlimm=xlim;
            title('Single element');
            
            subplot(2,1,2);
            set(gcf,'units','norm','pos',[ 0.1589    0.3345    0.2216    0.6102]);
            yout=yout(Imax:end,:);
            
            [Bbeam,FF1,TT] = specgram(yout(:,1),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
            
            imagesc(TT,FF1,20*log10(abs(Bbeam).*normm));%
            orient landscape
            grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
            axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);
            xlim(xlimm)
            title(sprintf('Mode 1 filter, tilt %6.2f',tilt_offset(Itilt)));
            
            if Itilt==1
                disp('Select two times for signal estimate...');
                
                timmes_signal=ginput(2);
            end
            In=round(Fs*timmes_signal(:,1));
            dn=diff(In);
            pwrr.single.signal(Itilt)=sum(y(1,In(1):In(2)).^2)/dn;  %This has been prefiltered...
            pwrr.filtered.signal(Itilt)=sum(yout(In(1):In(2),1).^2)/dn;
            
            y1=yout(In(1):In(2),1);
            y2=yout(In(1):In(2),2);
            pwrr.filtered.xcor(Itilt)=y1'*y2./(norm(y1).*norm(y2));
            
            %             figure(11)
            %             subplot(2,1,1)
            %             hold on
            %             plot(y(1,In(1):In(2)),'r');
            %             subplot(2,1,2);
            %             hold on;
            %             plot(yout(In(1):In(2),1),'r');
            
            %         soundsc(yout(In(1):In(2),1),Fs)
            %         soundsc(y(1,In(1):In(2)),Fs)
            %
            %             figure(10);
            %             if Itilt==1
            %                 disp('Select two times for background noise estimate...');
            %                 timmes_noise=ginput(2);
            %             end
            %             In=round(Fs*timmes_noise(:,1));
            %             dn=diff(In);
            %             pwrr.single.noise=sum(y(1,In(1):In(2)).^2)/dn;  %This has been prefiltered...
            %             pwrr.filtered.noise=sum(yout(In(1):In(2),1).^2)/dn;
            %
            %             figure(11)
            %             subplot(2,1,1)
            %             plot(y(1,In(1):In(2)),'k');
            %             subplot(2,1,2);
            %             plot(yout(In(1):In(2),1),'k');
            
            %             fprintf('SNR for bottom sensor, tilt %6.2f: %6.2f, power: %6.2f \n',  ...
            %                 tilt_offset(Itilt),10*log10(abs(pwrr.single.signal(Itilt)/pwrr.single.noise)),10*log10(pwrr.single.signal(Itilt)));
            %             fprintf('SNR for mode 1 filter, tilt %6.2f: %6.2f, power: %6.2f \n', ...
            %                 tilt_offset(Itilt),10*log10(abs(pwrr.filtered.signal(Itilt)/pwrr.filtered.noise)),10*log10(pwrr.filtered.signal(Itilt)));
            %
            
            fprintf('Xcor for bottom sensor, tilt %6.2f: %6.4f, power: %6.2f \n',  ...
                tilt_offset(Itilt),abs(pwrr.filtered.xcor(Itilt)),10*log10(pwrr.single.signal(Itilt)));
            
            if Itilt==length(tilt_offset)
                figure(45)
                plot(tilt_offset,10*log10(pwrr.filtered.signal));grid on
                xlabel('Tilt Range Estimate');
                ylabel('dB re 1uPa');
                title('Measured power of filtered mode 1 signal');
                keyboard
            end
            
        end
        
    end  %plot_mode_filtering..

%Inner function for plotting xcorr
    function plot_xcorr
        pause(0.5);
        MMM(Itilt)=getframe(gcf);
        figure(3)
        % for Im=1:3
        subplot(3,1,Im)
        plot(lags/Fs,xc);xlim(10*[-0.05 0.05]);
        grid on
        title(sprintf('Mode %i, leads mode 1 by %10.6f ms',Im,1000*talign{Imodel}(Itilt,Ir,Im)));
        xlabel('Relative delay (s)');
        ylabel('Cross-correlation');
        %end
    end
end  %function modal_filtering

%%%%%%%%%
%%%%%matched_field_processor(model_name,tilt_offset,ranges,depths,Ksout.Kstot,Ksout.freq,head.rd);
%%%%%%%%%

function save_result=matched_field_processor(model_name,tilt_offset,ranges,depths,Kstot,freq,SNR,rd,SNRmin,data_name)
conj_flag=1;  %If one, conjugate data...

%%%%%Prepare data
%SNRmin=15;keyboard
%%%If SNRmin is a sclar, it is a dB SNR value threshold to select frequencies
%       If SNRmin is a vector, it is a list of frequencies.
if length(SNRmin)==1
    Igood=find(SNR>=SNRmin);
    freq=freq(Igood);
    disp(sprintf('Frequencies to process at %i dB SNR: %s',SNRmin,mat2str(freq,4)));
    
else
    
    Igood=zeros(length(SNRmin),1);  %Iwant will be indicies in model frequencies that match FFT bins in index
    
    for If=1:length(SNRmin)
        [junk,Igood(If)]=min(abs(SNRmin(If)-freq));
    end
    freq=freq(Igood);
    disp(sprintf('Frequencies to process, manually selected: %s',mat2str(freq,4)));
    if length(unique(Igood))~=length(Igood),disp('WARNING! Frequencies requested via SNRmin not available...');end
end

Kstot=Kstot(:,:,Igood);
disp(sprintf('Loading model: %s',model_name));
model=load(model_name);


%%Extract model frequencies...

Iwant=zeros(length(freq),1);  %Iwant will be indicies in model frequencies that match FFT bins in index
for If=1:length(freq)
    [junk,Iwant(If)]=min(abs(model.freq-freq(If)));
end


%%%%%Prepare model%%%%
%keyboard
%%Remove bad channels, if any
Igood=zeros(1,length(rd));
for I=1:length(rd)
    [junk,Igood(I)]=min(abs(model.rd-rd(I)));
end

Igood_source=zeros(1,length(depths));
for I=1:length(depths)
    [junk,Igood_source(I)]=min(abs(model.rd-depths(I)));
end

model.sd=model.rd(Igood_source);
model.Urd=model.U(Igood,:,:);

model.U=model.U(Igood_source,:,:);
model.rd=model.rd(Igood);

kk=2*pi*model.freq(Iwant)./1495;

%Mistake in linear tilt calculation
%rd_cum=rd./cumsum(rd);

rd_cum=-(rd-min(rd))./(max(rd)-min(rd));
%%generate replicas
% weight modes by mode excitation

% Irtest=25;
% Iztest=5;
for Itilt=1:length(tilt_offset)
    fprintf('tilt: %6.2f\n',tilt_offset(Itilt));
    tilt_matrix=exp(-1i*tilt_offset(Itilt)*rd_cum'*kk); %[depths frequency]
    
    for If=1:length(freq)
        ck=(model.kr{1}(Iwant(If),:).');
        tiltt=tilt_matrix(:,If);
        Ks=squeeze(Kstot(:,:,If));
        Ks=Ks/trace(Ks);
        if conj_flag==1
            Ks=conj(Ks);  %Need to do this since MATLAB and SAGA definitions of Fourier Transform different than Comp. Ocean. Acoust.
        end
        for Iz=1:length(depths)
            
            [junk,isd]=min(abs(model.sd-depths(Iz)));
            a = squeeze(model.U( isd, Iwant(If), : ));
            phi = squeeze(model.Urd(:,Iwant(If),:)) * diag( a, 0 );	% scale modes by a
            
            % ******************************************************
            % form pressure field
            % ******************************************************
            
            phase = diag( 1.0 ./ sqrt( ck ) ) * exp( 1i * ck * ranges ) * diag( sqrt( 2 * pi ./ ranges ) );
            
            Ibad = find(isnan(real(phase)));
            phase(Ibad)=0;
            
            p0 = phi * phase;  % [rd Nm] x [Nm ranges] = [rd ranges]
            
            p=(tiltt*ones(1,length(ranges))).*p0;
            %Add tilt correction here..
            
            %temp=squeeze(Kstot(:,:,If))*p;  % [rd rd] x [rd ranges] =[rd ranges]
            for Ir=1:length(ranges)
                replica=p(:,Ir)/norm(p(:,Ir));
                amb(Iz,Ir,If)=real(replica'*Ks*replica);
                
                %             if Iztest==Iz&&Irtest==Ir
                %                 testdata(:,:,If)=p(:,Ir)*p(:,Ir)';
                %             end
            end
            
        end
    end
    
    max_corr=squeeze(max(max(amb)))';
    amb_tot=sum(amb,3)/length(freq);
    %Iamb=find(max_corr>0.7);amb_tot=sum(amb(:,:,Iamb),3)/length(Iamb);  %Trim poor results
    figure
    %     for I=1:length(freq)
    %         imagesc(ranges,depths,squeeze(amb(:,:,I)));
    %         colorbar;caxis([0 1]);
    %         title(sprintf('Frequency: %6.2f',freq(I)));
    %         pause
    %     end
    %imagesc(ranges,depths,amb_tot);caxis([0 1]);
    contourf(ranges/1000,depths,amb_tot,0:0.1:1);
    caxis([0.2 1]);
    cmap=colormap;colormap(cmap(1:4:end,:));
    set(gca,'fontweight','bold','fontsize',14);xlabel('range (km)');ylabel('depth (m)');grid on
    colorbar;title(sprintf('Frequencies: %s',mat2str(round(freq),4)));
    axis('ij')
    if length(SNRmin)>1
        SNRmin=NaN;
    end
    if length(freq)<6
        titlestr=sprintf('Tilt: %6.2f, Frequencies to process at %i dB SNR: %s',tilt_offset(Itilt),SNRmin,mat2str(freq,4));
    else
        titlestr=sprintf('Tilt: %6.2f, %i frequencies to process at %i dB SNR',tilt_offset(Itilt),length(freq),SNRmin);
        
    end
    %     if conj_flag==1
    %         titlestr=[titlestr ' Conjugated Kstot'];
    %
    %     end
    title(titlestr);
    tmp=[freq' squeeze(max(max(amb)))];
    fprintf('Frequency: %6.2f, Correlations:  %6.2f\n',tmp');
    fprintf('Total corr: %6.2f \n',max(max(amb_tot)));
    
    save_result=input('Enter any character to save images and write CSDM to file: ');
    if ~isempty(save_result)
        savename=sprintf('MFP_result_%s_tilt%4.2f_Nfreqs%i_SNRmin%i',data_name,tilt_offset(Itilt),length(freq),SNRmin);
        save([savename '.mat'], 'ranges','depths','amb','amb_tot','SNRmin','freq','data_name','model_name');
        
        figure(gcf)
        orient tall
        print('-djpeg',[savename '.jpg']);
    end
end %tilt
end

function [MFP_replica_file_out,range_str_out,default_tilt_out,default_SNR_out]=load_MFP_scenario

%%MFP replicas can be divided into Pekeris models, and more complex models...

if strcmp(getenv('USER'),'thode')
    model_dir='/Users/thode/Projects/Arctic_2010/PropagationModeling.dir/';
elseif strcmp(getenv('USER'),'shabadi');
    model_dir='/Users/shabadi/Documents/Shima-SIO/PropagationModeling.dir/';
end
%%%%%%%%Short-Range Estimate%%%%%%%%%%%%%%%%
MFP_replica_file{2}{1}=[model_dir 'Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_10to500Hz.mat'];
default_tilt{2}{1}='1.94';
range_str{2}{1}='100:5:5000';
default_SNR{2}{1}=[45.78 76.29 88.5 170.9 210.6 253.3];

MFP_replica_file{2}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_ShortRange_arctic_pekeris1707_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{2}{2}='2.41';
range_str{2}{2}=range_str{2}{1};
default_SNR{2}{2}=default_SNR{2}{1};

%%%%%%%%%%%%%Medium range 7 km est%%%%%%%%%%%%%%%%
MFP_replica_file{3}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_10to500Hz.mat'];
default_tilt{3}{1}='1.35';
range_str{3}{1}='5000:100:10000';
%default_SNR{3}{1}=[112.9 116 119 122.1 164.8 167.8 170.9 174 177 180.1 235 238];  %Original 12 frequencies
%used for inversion
default_SNR{3}{1}=[116 119 122.1 164.8 167.8 170.9 174 177 180.1 ]; %Plotting results

MFP_replica_file{3}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{3}{2}='1.585';
range_str{3}{2}=range_str{3}{1};
default_SNR{3}{2}=default_SNR{3}{1};

%%%%%%%%%%%%%%%%%Long range estimate (17 km)%%%%%%%%%%%%%%%%%%%%%%%%
MFP_replica_file{4}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_10to500Hz.mat'];
default_tilt{4}{1}='1.35';
range_str{4}{1}='5000:100:20000';
%default_SNR=12;
default_SNR{4}{1}=[74.01 78.58 85.45 88.5 90.03 92.32 93.08 95.37 99.95 103 107.6 109.9];

MFP_replica_file{4}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_FarRange_arctic_pekeris1692_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{4}{2}='1.585';
range_str{4}{2}=range_str{4}{1};
default_SNR{4}{2}=default_SNR{4}{1};

%%%%%%%%%%%%%%%%%35 km range estimate%%%%%%%%%%%%%%%%
MFP_replica_file{5}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105748_10to500Hz.mat'];
default_tilt{5}{1}='3.93';
range_str{5}{1}='30000:100:40000';
default_SNR{5}{1}=[99.1800  102.2000  105.3000  108.3000  111.4000  114.4000  117.5000  120.5000 ...
    123.6000  126.6000  129.7000  132.8000  135.8000  138.9000  141.9000  145.0000 ...
    148.0000  151.1000  154.1000  157.2000  160.2000 ];
%default_SNR=default_SNR(1:(end-1));
% default_SNR=default_SNR((end-3):end)

MFP_replica_file{5}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_UltraRange_arctic_pekeris1553_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{5}{2}='3.34';
range_str{5}{2}=range_str{5}{1};
default_SNR{5}{2}=default_SNR{5}{1};


%Averaged Pekeris model...
for I=2:5
    %%Averaged Pekeris model
    MFP_replica_file{I}{3}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortPekeris_arctic_pekeris1673_KRAKEN_flat55m_10to500Hz.mat'];
    default_tilt{I}{3}=default_tilt{I}{2};
    range_str{I}{3}=range_str{I}{2};
    default_SNR{I}{3}=default_SNR{I}{2};
    
end

yes=menu('Which inversion?','Blank...','Whale Close Range inversion', ...
    'Whale 7km Range inversion','Whale 15km Range inversion','Whale 35km Range inversion');

MFP_replica_file_out=MFP_replica_file{yes};
default_tilt_out=default_tilt{yes};
range_str_out=range_str{yes};
default_SNR_out=default_SNR{yes};

yes=menu('Which model?','Detailed model','Inverted Pekeris','Average Pekeris','All');

if yes<4
    MFP_replica_file_out=MFP_replica_file_out(yes);
    default_tilt_out=default_tilt_out(yes);
    range_str_out=range_str_out(yes);
    default_SNR_out=default_SNR_out(yes);
    
end



end

function [mean_corr,tindex,TT_plot,pwr,pwr_tot,yscale]= create_incoherent_correlogram(TT,FF,B,param,flo,fhi)
%function [mean_corr,tindex,TT_plot,pwr,pwr_tot,yscale]= create_incoherent_correlogram(TT,FF,B,param,flo,fhi)

ici_range=param.ici_range;  %Autocorrelation lag to examine
time_sample=param.time_sample;  %sec,  amount of spectrogram time columns to process for autocorrelation..
ovlap=param.ovlap;  %How much to shift the autocorrelation window between samples
teager=param.teager;  %Apply teager-kaiser operation to spectrogram before processing...

% A single event may thus persist for up to time_sample/dX times, 1/(1-ovlap) times



%freq=flo:bandwidth:fhigh;
B=10*log10(B+eps);
if teager
    B=B(:,2:end-1).^2-(B(:,1:end-2).*B(:,3:end));
end
dT= TT(2)- TT(1);
dF= FF(2)- FF(1);
Ncol=length(find( TT<=time_sample));  %Number of columns to process at a time:
Iindex=1:floor((1-ovlap)*Ncol):length( TT);
tindex=dT*Iindex;  %Time axis of absolute time
dX=tindex(2)-tindex(1);  %X axis of new correlation image
Ntime=length(Iindex)-1;
yscale=dT*(0:(Ncol-1));

Igood=find( FF>flo& FF<fhi);
maxlag = Ncol - 1;  %maximum autocorrelation lag
laggs=-maxlag:maxlag;
I_range=find(laggs*dT>ici_range(1)&laggs*dT<ici_range(2));
TT_plot=laggs(I_range)*dT;
mean_corr=zeros(length(I_range),Ntime);
pwr=zeros(Ntime,Ncol);
pwr_tot=zeros(Ntime);
xcovv=zeros(length(Igood),length(I_range));

for I=1:Ntime
    if rem(I,100)==0,disp(sprintf('%6.2f percent done',100*I/Ntime));end
    
    try
        A= B(Igood,Iindex(I)+(0:(Ncol-1)));
        
        [N,M] = size(A);  %Autocorrelate columns
        X=fft(A'-ones(M,1)*mean(A'), 2^nextpow2(2*M-1));
        tmp=ifft(abs(X).^2);
        
        tmp1=[tmp(end-maxlag+1:end,:);tmp(1:maxlag+1,:)];
        tmp=tmp1./(ones(size(tmp1,1),1)*tmp1(maxlag+1,:));
        
        xcovv=tmp(I_range,:);
        %end
        mean_corr(:,I)=median(xcovv')';
        pwr(I,:)=sum(A)/length(Igood);
        pwr_tot(I)=max(pwr(I,:));
    catch
        xcovv=[];
        %mean_corr(I,:)=[];
    end
    %     figure(2);
    %     subplot(2,1,1)
    %     imagesc([], FF(Igood)/1000,A);
    %     subplot(2,1,2)
    %     imagesc(TT_plot, FF(Igood)/1000,xcovv{I});
    %     caxis([0 1]);
    %     colorbar
    %     pause(0.25);
    
    
end  %Ntime
end


function [XC_eq,XC,Trel,tt,pwrr]= create_coherent_correlogram(x,fs,param,flo,fhi)
%function [XC_eq,XC,Trel,tt,pwr]= create_coherent_correlogram(x,fs,param,flo,fhi)
% Trel: relative time
% tt: lag times...

ovlap=round(param.ovlap*param.Nfft);
Nfft=param.Nfft;
teager=param.teager;  %Apply teager-kaiser operation to spectrogram before processing...
ici_range=param.ici_range;

%noiseT=param.noiseT;
%alpha=param.alpha;
noiseT=0.01;
alpha=0.01;

frange=[0.8*flo flo fhi fhi+0.2*flo];
[N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.05],fs);
Bfilt = firpm(N,Fo,Ao,W);
y=filtfilt(Bfilt,1,x-mean(x));
%Adaptive filter option

if teager==1
    y=y(2:end-1).^2-(y(1:end-2).*y(3:end));
    
end

index1=1:round((Nfft-ovlap)):length(y);
index1=index1(index1<=length(y));
Ncoll=length(index1)-1;
XC=zeros(Nfft+1,Ncoll);
pwrr=zeros(1,Ncoll);
for I=1:Ncoll
    try
        Indexx=index1(I)+(0:(Nfft-1));
        tmp=xcov((hilbert(hanning(length(Indexx)).*y(Indexx))),'coeff');  %Delphine, note that this should be windowed...
        
        XC(:,I)=tmp((Nfft-1):end);
        pwrr(I)=tmp(Nfft-1);
        %%Mulitpath can generate strong correlations.  Surface-reflected paths should have negative
        %%correlations.  Therefore, only keep levels greater than zero.
        %Iplus=find(XC(:,I)<=0);
        %XC(Iplus,I)=0;
    catch
        disp('Problem');
    end
end
XC=10*log10(abs(XC));

tt=(0:Nfft)/fs;
dT=ovlap/fs;
%Tabs=twant(Itimes)+datenum(0,0,0,0,0,index1/fs);
Trel=index1/fs;

%%Adaptively equalize along y-axis
Inoise=round(fs*noiseT);  %y-axis bins to estimate autocorrelation level
XC_eq=zeros(size(XC));
eq=mean(XC(2:Inoise,:));  %Original equalization estimate

for I=1:size(XC,1)
    eq=(1-alpha).*eq+alpha.*XC(I,:);
    XC_eq(I,:)=XC(I,:)-eq;
end

%%Trim away portions of autocorrelation that will not have creaks (ici_range)
Iwant=find(tt>ici_range(1)&tt<ici_range(2));

%Remove times not of interest, normalize zero lag to one.
XC=XC(Iwant,:);
XC_eq=XC_eq(Iwant,:);
%XC0=XC0(Iwant,:);
%signn=signn(Iwant,:);
tt=tt(Iwant);
%Quick click removal...





end  %create_coherent_correlogram




function  [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,linel)

%LL=3;


if nargin==3,
    VM=[];
    A=[];
    B=[];
    ANG=[];
elseif nargin==4,
    A=[];
    B=[];
    ANG=[];
end
%Convert to km
VM=VM/1000;
DASAR_coords=DASAR_coords/1000;
A=A/1000;
B=B/1000;

if ~exist('linel')
    linel=35; %length of bearing lines in km
end
%subplot(3,1,LL);
xg=mean(DASAR_coords(:,1));
yg=mean(DASAR_coords(:,2));
DASAR_coordsn(:,1)=DASAR_coords(:,1)-xg;
DASAR_coordsn(:,2)=DASAR_coords(:,2)-yg;

plot(DASAR_coords(:,1)-xg,DASAR_coords(:,2)-yg,'r^','markersize',5,'markerfacecolor',[1 0 0]);hold on
set(gca,'fontweight','bold','fontsize',14);
xlabel('Easting (km)');
ylabel('Northing (km)');
grid on;

%Convert bearings from nautical to mathematical frame.
if ~isempty(bearings)
    bearings=(90-bearings(Igood))*pi/180;
    for I=1:length(Igood)
        XX=DASAR_coords(Igood(I),1)+[0 linel*cos(bearings(I))]-xg;
        YY=DASAR_coords(Igood(I),2)+[0 linel*sin(bearings(I))]-yg;
        line(XX,YY);
    end
end
if ~isempty(VM)
    VMn=[VM(:,1)-xg VM(:,2)-yg];
    plot(VMn(1),VMn(2),'ks','markerfacecolor',[0 0 0],'markersize',5);
end


%Plot error elipps
if ~isempty(A)
    %ELLIPSE(ra,rb,ang,x0,y0)
    h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
    set(h,'linewidth',0.5);
end

hold off;

end

function h=ellipse(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same.
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value.
% If no color is specified, it makes automatic use of the colors specified by
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% as a sample of how ellipse works, the following produces a red ellipse
% tipped up at a 45 deg axis from the x axis
% ellipse(1,2,pi/8,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%

% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original
% written by Peter Blattner, Institute of Microtechnology, University of
% Neuchatel, Switzerland, blattner@imt.unine.ch


% Check the number of input arguments

if nargin<1,
    ra=[];
end;
if nargin<2,
    rb=[];
end;
if nargin<3,
    ang=[];
end;

%if nargin==1,
%  error('Not enough arguments');
%end;

if nargin<5,
    x0=[];
    y0=[];
end;

if nargin<6,
    C=[];
end

if nargin<7,
    Nb=[];
end

% set up the default values

if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;

% work on the variable sizes

x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);

if isstr(C),C=C(:);end;

if length(ra)~=length(rb),
    error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
    error('length(x0)~=length(y0)');
end;

% how many inscribed elllipses are plotted

if length(ra)~=length(x0)
    maxk=length(ra)*length(x0);
else
    maxk=length(ra);
end;

% drawing loop

for k=1:maxk
    
    if length(x0)==1
        xpos=x0;
        ypos=y0;
        radm=ra(k);
        radn=rb(k);
        if length(ang)==1
            an=ang;
        else
            an=ang(k);
        end;
    elseif length(ra)==1
        xpos=x0(k);
        ypos=y0(k);
        radm=ra;
        radn=rb;
        an=ang;
    elseif length(x0)==length(ra)
        xpos=x0(k);
        ypos=y0(k);
        radm=ra(k);
        radn=rb(k);
        an=ang(k)
    else
        rada=ra(fix((k-1)/size(x0,1))+1);
        radb=rb(fix((k-1)/size(x0,1))+1);
        an=ang(fix((k-1)/size(x0,1))+1);
        xpos=x0(rem(k-1,size(x0,1))+1);
        ypos=y0(rem(k-1,size(y0,1))+1);
    end;
    
    co=cos(an);
    si=sin(an);
    the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
    %  x=radm*cos(the)*co-si*radn*sin(the)+xpos;
    %  y=radm*cos(the)*si+co*radn*sin(the)+ypos;
    h(k)=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
    set(h(k),'color',C(rem(k-1,size(C,1))+1,:));
    
end;

end

function c = cond2(A)
%COND2   Condition number with respect to inversion.
%   COND2(X) returns the 2-norm condition number (the ratio of the
%   largest singular value of X to the smallest).  Large condition
%   numbers indicate a nearly singular matrix.
%
%   Modified to handle only the 2-norm, eliminate some error checking,
%   and other condition testing.
%      Chris Nations  8/20/02

s = svd(A);
if any(s == 0)   % Handle singular matrix
    c = Inf;
else
    c = max(s)./min(s);
end
end

function [area,a,b,ang,angax] = ellipsparms(S,critval,MD,VM)
%ELLIPSPARMS  Calculates parameters of bivariate normal ellipse from covariance.
%   [AREA,A,B,ANG] = ELLIPSPARMS(S,CRITVAL,MD,VM) where
% S is the covariance matrix,
% CRITVAL is the critical value from the chi-squared distribution,
% MD is the mean position of the DASAR array,
% VM is the estimated location.
% It returns
%   AREA of the ellipse,
%   A and B are lengths of major and minor full axes respectively, and
% ANG is angle in radians ccw of  major semi-axis.
%   The angle is referenced to an imaginary vector emanating from the center of the
%   DASAR array and directed toward the estimated location.
%   Note: index matrices using vec representation.
% ANGAX is angle in radians of major axis measured ccw from +X axis
%
% The adjustments to ANG (last set of if/then statements) reflect the
%   angle of the major axis 180 degrees about the ray directed from MD to
%   VM (if necessary) so that the major axis is always directed away from
%   MD.  The absolute orientation of the ellipse is unchanged; only its
%   "direction" relative to the center of the DASAR array may be adjusted.

dp = VM-MD;
theta = atan2(dp(2),dp(1));     % 4-quadrant inverse tangent: angle from MD to VM
detS = det(S);
if isnan(S(1)) | (detS<=0) | min(diag(S))<0,
    area = nan;  a = nan;  b = nan;  ang = nan;  angax = nan;
else
    area = pi*sqrt(detS)*critval;
    [V,L] = eig(S);
    d = 2*sqrt(critval*diag(L));  % Lengths of both axes.
    [d,i] = sort(d);
    a = d(2);  b = d(1);          % a & b: major & minor axis lengths, respectively.
    v = V(:,i(2));                % dominant eigenvector
    angax = atan(v(2)/v(1));      % Absolute orientation of ellipse major axis
    if angax<0,                   % Returns an angle between 0 and pi following
        angax = angax+pi;           %    math convention (0 radians = 90 degrees, East;
    end                           %    pi radians = 270 degrees, West).
    ang = angax;
    if ((theta>-pi/2) & (theta<pi/2)  & (ang>(pi/2+theta))), % Adjust angle for
        ang = ang-pi;                                          %   estimated position
    elseif ((theta>=pi/2) & (ang<(theta-pi/2))),             %   relative to center
        ang = ang+pi;                                          %   of DASAR array.
    elseif ((theta<=-pi/2) & (ang<(3*pi/2+theta))),
        ang = ang+pi;
    end
end
end

function [VM,Qhat,w,outcome] = vmmle_r(angle,dasar,r,k)
%[VM,QHAT,W,OUTCOME] = VMMLE(ANGLE,DASAR,R,K)
%   Von Mises MLE estimate of location of animal target, based on:
%
%   Lenth,R.V. 1981.  On finding the source of a signal.
%     Technometrics 23:149-154.
%
%   [VM,QHAT,W,QHAT] = VMMLE(ANGLE,DASAR,R,K) processes a single observation
%   comprised of ANGLE (an n-by-1 column vector of bearings), DASAR
%   (an n-by-2 matrix of coordinates for the DASAR source of each
%   bearing), and K (an n-by-1 vector of bearing standard error
%   estimates, expressed as the von Mises concentration parameter kappa,
%   associated with each DASAR).  R is a character code for method -
%   either 'm' for the standard (non-robust) MLE, 'a' for the Andrews
%   robust procedure, or 'h' for the Huber robust procedure.
%   VM is the von Mises estimate of location, QHAT is the associated
%   covariance matrix estimate and W is a 1-by-n row vector of weights:
%   for the MLE, either all 1's if the procedure converged or 0's if not;
%   for the Andrews, generally a mix of 1's and 0's; and for the Huber,
%   generally a mix of values between 0 and 1.  OUTCOME is a character
%   array identifying whether a good solution was found, and if not,
%   the reason for failure.

if nargin<3,
    r = 'm';
end
if nargin<4||all(k==0)          % Estimate kappa from the data
    kest = logical(1);
else
    kest = logical(0);
end
r = lower(r);
robust = strfind('mah',r);
if isempty(robust),
    error('Input r must be either "m", "a", or "h"');
end
failed = strvcat('less than 2 bearings','negative variance estimates',...
    'solution behind DASARs','failed to converge');
n = length(angle);

if n<=1  %If only one set of bearings present...
    outcome = failed(1,:);
else
    cond_num = 1e15;       % For test of singularity.
    tc = 1.5;              % Tuning constant for robust versions.
    dist1 = 1;             % Was 0.1
    maxiter = 50;
    x = dasar(:,1);
    y = dasar(:,2);
    iter = 0;
    theta = (90-angle)*pi/180;  %Mathematical angle defined
    theta = (theta<-pi)*2*pi + theta;
    s = sin(theta);
    c = cos(theta);
    z = s.*x - c.*y;
    sstar = s';
    cstar = c';
    w = ones(1,n);
    M1 = [sstar; -cstar]*[s -c];
    converge = 0;
    if cond2(M1)<cond_num,
        M2 = [sstar; -cstar]* z;
        xyhat = M1\M2;
        while (~converge)&(iter<maxiter),
            iter = iter+1;
            xyold = xyhat;
            % d = dist(xyhat', [x y]');
            
            for JJ=1:length(x)
                d(JJ)=sqrt((xyhat(1)-x(JJ)).^2+(xyhat(2)-y(JJ)).^2);
            end
            if (robust>1) & (n>2), % Need 3 or more bearings to calculate weights for
                dxy = repmat(xyhat',n,1)-[x y];      % robust methods
                muhat = cart2pol(dxy(:,1),dxy(:,2));
                Cd = cos([theta-muhat]');
                if kest,
                    Cbar = abs(w*Cd'/sum(w))^(n/(n-2));  % Abs is ad hoc but may avoid temporary numeric problems in iteration
                    k = inv(2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar);
                end
                t = sqrt(2*k'.*(1-Cd));
                if robust==2,
                    phit = tc*sin(t/tc).*(abs(t)<tc*pi);
                else
                    phit = sign(t).*min(abs(t),tc);
                end
                same = (Cd==1);     % Take care when estimated & observed bearings
                t(same) = 1;        %   are identical; avoid division by zero.
                phit(same) = 1;     %   Usually occurs with just 2 intersecting bearings.
                w = phit./t;
            end    %  if robust>1
            sstar = w.*(xyhat(2)-y')./(d.^3);
            cstar = w.*(xyhat(1)-x')./(d.^3);
            M1 = [sstar; (-cstar)]*[s -c];
            if ((n-sum(~w))>1) & (cond2(M1)<cond_num),
                M2 = [sstar; -cstar]*z;
                xyhat = M1\M2;
                converge = sum(abs(xyhat-xyold)<dist1)==2;
            else
                break    % If either condition above occurs, convergence will very
            end        %    likely fail, so break out of while loop.
        end  % while
    end    % if cond2
    if converge,
        dxy = repmat(xyhat',n,1)-[x y];
        muhat = cart2pol(dxy(:,1),dxy(:,2));
        Cd = cos([theta-muhat]');
        if kest & (n>2),
            Cbar = (w*Cd'/sum(w))^(n/(n-2));   % Exponent is small sample size correction
            k = inv(2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar);
        end
        if Cd*w'>0,             % Weighted average of cosine differences (check
            VM = xyhat';          %   on bad solutions behind DASARs)
            if kest & (n==2),     % Cannot estimate Qhat with only 2 bearings
                Qhat = nan*ones(2);
                outcome = 'successful; 2 bearings; no Kappa';
            else
                k = k';
                cv = -(k.*sstar*c + k.*cstar*s)/2;
                M3 = [k.*sstar*s cv; cv k.*cstar*c];
                Qhat = inv(M3);
            end
            if ~kest | (n>2),
                if all(diag(Qhat)>0),
                    outcome = 'successful';       % Successful solution
                else
                    outcome = failed(2,:);        % Implausible variance estimate(s)
                end
            end
        else
            outcome = failed(3,:);          % Bad solution behind DASARs
        end % if all(Cd>0)
    else
        outcome = failed(4,:);            % No convergence
    end   % if converge
end     % if n<=1
if ~isempty(strmatch(outcome,failed)),
    VM = [nan nan];
    Qhat = nan*ones(2);
    w = zeros(1,n);
end
end

%function [brefa_table,Icol]=calibrate_bearing_Shell2007(cal_dir,fname,no_load_table)
% Return nonlinear calibration data for 2007 Dasars
% Input:
%          cal_dir: base directory of calibration information
%          fname: string similar in structure to 'S108H0T20080921T000000'
%          no_load_table:  if exists, don't load table.
% Note:  the calibration is assigned based on name of file, not name of enclosing directory..

function [brefa_table,Icol]=calibrate_bearing_Shell2007(cal_dir,fname,no_load_table)
brefa_table=[];
Islash=1+max(findstr(fname,'/'));
if ~isempty(Islash)
    fname=fname(Islash:end);
end
site=str2num(fname(2));
dasar=lower(fname(5));
Icol=double(dasar)-96;


%%Date checking...
%% In 2007 if a DASAR was redeployed it was reassigned a number...
% if site==3&&strcmp(dasar,'g')
%     datte=datenum(fname(8:end),'yyyymmddTHHMMSS');
%     if datte>=datenum(2007,9,18,3,15,0)
%         Icol=8;
%     end
% elseif site==5&&strcmp(dasar,'d')
%     datte=datenum(fname(8:end),'yyyymmddTHHMMSS');
%     if datte>=datenum(2007,9,18,19,0,0)
%         Icol=8;
%     end
% elseif site==5&&strcmp(dasar,'e')
%     datte=datenum(fname(8:end),'yyyymmddTHHMMSS');
%     if datte>=datenum(2007,9,19,0,0,0)
%         Icol=9;
%     end
% end


if nargin==2
    cal_dir=sprintf('%s/S%i07gsif/Caloutput07s%i/Bgrid07s%i.mat',cal_dir,site,site,site);
    load(cal_dir);
    brefa_table=bgridt;
    
    %     if site==1  %DASAR 1B failed, not reflected in bgrid table...
    %        brefa_table=[brefa_table(:,1) zeros(size(bgridt,1),1) brefa_table(:,2:end)];
    %     end
end

end

function status=write_covmat(freqvec,Cov,depth,srcname,titlestr)
% write_covmat(farray,Ks,depth,srcname,titlestr)
%
% Write Covariance Matrices to file cov_dpss.in
%
% farry = vector of frequencies in Hertz
% Ks     = matrix of concatenated covariance matrices
%           this must be of size [Nel Nel nfreq]
% depth  = hydrophone depths.  First element is shallowest depth
% srcname =file name
% titlestr = title to write to top of input file.
%[n]= size(Cov,1);
%nf=size(Cov,3);
%if(nf ~= length(freqvec))
%   fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
%           size(Cov),length(freqvec));
%  error('sizes of Cov and freqvec are not compatible');
%end

[n,m,F]= size(Cov);
if(m ~= n)
    fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
        size(Cov),length(freqvec));
    error('sizes of Cov and freqvec are not compatible');
end
if(F ~=length(freqvec))
    fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
        size(Cov),length(freqvec));
    error('sizes of 3Cov and freqvec are not compatible');
end

if(n ~= length(depth))
    fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
        size(Cov),length(freqvec));
    error('sizes of Cov and depth are not compatible');
end


%depth=rd;
%disp('WARNING! CONJUGATING COV to be used by SAGA!');
%Cov=conj(Cov);

fidout =  fopen(sprintf('%s.in',srcname),'w');

for ifreq=1:length(freqvec)
    fprintf(fidout,' %s\n',titlestr);
    fprintf(fidout,' %f     0.000 dB\n',freqvec(ifreq));
    fprintf(fidout,' %d\n',n);
    fprintf(fidout,' %8.2f\n',depth);
    %   cols = [ (ifreq-1)*n+1:ifreq*n ];
    %   Cx = Cov(:,cols); % cut one cov-matrix out of the bunch
    Cx=squeeze(Cov(:,:,ifreq));
    for row=1:n
        for col=1:n
            fprintf(fidout,'%10d%10d (%E,%E) \n',row,col, ...
                real(Cx(row,col)),imag(Cx(row,col)));
        end
    end
end
fprintf(fidout,'!  \n');
status=fclose(fidout);
unix('ls -l *.in');
end


% --- Executes on button press in checkbox_teager.
function checkbox_teager_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_teager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_teager
end



% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	Only bother, if specgram data exists
%	! need to add check for whether it is current
if isfield(handles, 'Specgram')

	if ~isfield(handles, 'hf_shipdet') || isempty(handles.hf_shipdet)
		handles.hf_shipdet	=	[];
	end

	%	External function to perform calculations and plots
	%	Can then be separated later
    %   Use default parameters for peak picking
	[metrics,handles.hf_shipdet]	=	...
	shipdet_metrics(handles.Specgram, [], handles.axes1, handles.hf_shipdet);
	
	guidata(hObject, handles);
end

end


% --- Executes during object deletion, before destroying properties.
function pushbutton17_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'hf_shipdet') && ~isempty(handles.hf_shipdet)
	hf	=	handles.hf_shipdet;
	close(hf(ishandle(hf)));
end
end
