
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

% Last Modified by GUIDE v2.5 03-Oct-2009 11:30:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
%Set GSIbearing button invisible
set(handles.pushbutton_GSIbearing,'Vis','off');

set(handles.edit_maxfreq,'String',num2str(10000/2-200));
set(handles.edit_minfreq,'String',num2str(200));
cd(mydir);
%%%Create list of function handles to be used in other specialized
%%%applications

cd(startupinfo.base_directory);
fnames=dir('*mat');
no_handles_stored=1;
for II=1:length(fnames),
    if strcmp(fnames(II).name,startupinfo.function_handles_filename),
        no_handles_stored=0;
    end
end
if no_handles_stored,
    III=1;function_handles{III}=@pushbutton_update_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@OpenMenuItem_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@PrintMenuItem_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@CloseMenuItem_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_datestr_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_datestr_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_windowlength_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_windowlength_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@popupmenu_Nfft_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@popupmenu_Nfft_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@popupmenu_ovlap_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@popupmenu_ovlap_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@slider_datestr_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@slider_datestr_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@pushbutton_playsound_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_mindB_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_mindB_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_dBspread_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_dBspread_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_minfreq_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_minfreq_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_maxfreq_CreateFcn;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@edit_maxfreq_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@pushbutton_playufsound_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@pushbutton_save_Callback;function_names{III}=func2str(function_handles{III});
    % III=III+1;function_handles{III}=@togglebutton_getwav_Callback;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@set_slider_controls;function_names{III}=func2str(function_handles{III});
    III=III+1;function_handles{III}=@pushbutton_print_Callback;function_names{III}=func2str(function_handles{III});
    save(startupinfo.function_handles_filename,'function_handles','function_names');
end
cd(mydir);

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
cla;

tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
mydir=pwd;
Ichan=str2num(get(handles.edit_chan,'String'));  %Hardwire first channel
[x,t,Fs,tstart]=load_data(handles.filetype,handles.min_time , tdate_start,tlen,Ichan,handles.mydir,handles.myfile);
if size(x,2)>1,
    x=x';
end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2num(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2num(contents{get(handles.popupmenu_ovlap,'Value')})/100;

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if strcmp(handles.display_view,'Spectrogram')
    %[B,FF,TT]=specgram(x(:,1),Nfft,Fs,hanning(Nfft),round(ovlap*Nfft));
    [S,FF,TT,B] = spectrogram(x(:,1),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
    %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
    axes(handles.axes1);
    imagesc(TT,FF/1000,10*log10(B));%
    axis('xy')
    fmax=str2num(get(handles.edit_fmax,'String'));
    if fmax==0,
        ylim([0 Fs/2000]);
        set(handles.edit_fmax,'String',num2str(Fs/2000));
    else
        ylim([0 fmax]);
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
    colorbar;
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
else %%Time series
    t=(1:length(x(:,1)))/Fs;
    plot(t,x(:,1));grid on;
    xlabel('Time (sec)');ylabel('Amplitude');
    
end

handles.x=x;
handles.Fs=Fs;
%tmp=ginput(2)

if strcmp(lower(handles.filetype),'gsi')
    set(handles.pushbutton_GSIbearing,'vis','on');
else
    set(handles.pushbutton_GSIbearing,'vis','off');
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
menustr={['*.' handles.filetype];['*.' lower(handles.filetype)];'*.*'};
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

tmin=datenum(get(handles.text_mintime,'String'));
tmax=datenum(get(handles.text_maxtime,'string'));

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
if tlen<120,
    maxx=get(handles.slider_datestr,'max');
    minn=get(handles.slider_datestr,'min');
    slider_step=get(handles.slider_datestr,'sliderstep');
    slider_step(1)=datenum(0,0,0,0,0,tlen)/(maxx-minn);
    set(handles.slider_datestr,'sliderstep',slider_step)
else
    errordlg('window length greater than 120 sec');
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
tmin=datenum(get(handles.text_mintime,'String'));
tmax=datenum(get(handles.text_maxtime,'string'));

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


end
function edit_minfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_minfreq as a double
minfreq=str2num(get(hObject,'String'));
maxfreq=str2num(get(handles.edit_maxfreq,'String'));
frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
[n,fo,mo,w] = remezord(frange, [0 1 0], [0.1 0.01 0.1], handles.Fs );
handles.b = remez(n,fo,mo,w);
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
maxfreq=str2num(get(hObject,'String'));
minfreq=str2num(get(handles.edit_minfreq,'String'));
frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
[n,fo,mo,w] = remezord(frange, [0 1 0], [0.1 0.01 0.1], handles.Fs );
handles.b = remez(n,fo,mo,w);
guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_playufsound.
function pushbutton_playufsound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playufsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.Fs>48000
    n=ceil(handles.Fs/48000);
    soundsc(downsample(handles.x,n),handles.Fs/n);
else
    disp('Resampling sound');
    soundsc(handles.x,handles.Fs);
end
end
% --- Executes on button press in pushbutton_playsound.
function pushbutton_playsound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.Fs>48000
    n=ceil(handles.Fs/48000);
    x=filter(handles.b,1,handles.x);
    soundsc(downsample(x,n),handles.Fs/n);
else
    disp('Resampling sound');
    soundsc(filter(handles.b,1,handles.x),handles.Fs);
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
[x,t,Fs]=load_data(handles.filetype,handles.min_time,tdate_start,tlen,Ichan,handles.mydir,handles.myfile);

Islash=findstr(handles.mydir,'/');
chan=get(handles.edit_chan,'String');

save_name=sprintf('soundsamp%s_%s',handles.mydir((Islash(end)+1):end),datestr(tdate_start,30));
disp(sprintf('Saving %s ...',save_name));
%wavwrite(x/(1.1*max(x)),Fs,save_name);
save(save_name,'x','Fs','tdate_start');
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
    [x,t,Fs,minn,maxx]=load_data(filetype,-1,-1,1,1,handles.mydir,handles.myfile);
catch %no file selected
    errordlg(sprintf('No %s file selected',filetype));
    errorflag=1;cd(mydir);
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

axes(handles.axes1);
tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
chan=get(handles.edit_chan,'String');
Islash=findstr(handles.mydir,'/');
save_name=sprintf('soundsamp%s_%s_%s',handles.mydir((Islash(end)+1):end),datestr(tdate_start,30),chan);
disp(sprintf('Printing %s ...',save_name));

orient landscape
print('-djpeg',save_name);
print('-dtiff',save_name);
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

fmax=str2num(get(hObject,'String'));
ylim([0 fmax]);

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
    
    startup_info.base_directory='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/RawData/S208gsif/S208G0';
    startup_info.default_directory=startup_info.base_directory;
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    startup_info.annotation_file='annotated.txt';
    
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    startup_info.calibration_DASAR2007_dir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/RawData';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
else
    startup_info.base_directory='/Volumes/';
    startup_info.default_directory='/Volumes/';
    startup_info.default_inputfiledir='/Volumes/';
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
[x,t,Fs,tstart]=load_data(handles.filetype,handles.min_time,tdate_start,tlen,Ichan,handles.mydir,handles.myfile);

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
fname_bounds=[fname_bounds Idot(1)];

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

function [x,t,Fs,tmin,tmax,head]=load_data(filetype,tstart_min,tdate_start,tlen,Ichan,mydir,myfile)
persistent fs keyword
x=[];
t=[];
tmin=[];
tmax=[];
head=[];
switch filetype
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
            keyword=input('Enter a keyword for GSI calibration [Shell08]:','s');
            if isempty(keyword)
                keyword='Shell08';
            end
        end
        x=calibrate_GSI_signal(x, keyword);
        
        Fs=head.Fs;
        
        tmin=datenum(1970,1,1,0,0,head.ctbc);
        tmax=tmin+datenum(0,0,1,0,0,0);
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
            [x,t]=load_mt([mydir '/' myfile],nsec,tlen,1);
        end
        x=1000*x;  %Converts millipascal to micropascal
    case 'DAT'
        [x,tmin,tmax,fs]=read_dat_file(mydir,myfile,[],-1,tlen,1);
        
        if isempty(fs)
            fs=input('Enter sampling rate in Hz:');
        end
        Fs=fs;
        [x,tmin,tmax]=read_dat_file(mydir,myfile,Fs,tdate_start,tlen,1);
        
    case 'MDAT'
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
        
    case 'WAV'
        sizz=wavread([mydir '/' myfile],'size');
        [Y,Fs]=wavread([mydir '/' myfile],1);
        handles.Fs=Fs;
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
            tdate_vec=datevec(tdate_start-tmin);
            nsec=tdate_vec(6)+60*tdate_vec(5)+3600*tdate_vec(4);
            N1=1+round(nsec*handles.Fs);
            N2=N1+round(tlen*handles.Fs);
            [x,Fs]=wavread([mydir '/' myfile],[N1 N2]);
            if ~strcmp(Ichan,'all')
                x=x(:,Ichan);
            end
            
            t=(1:length(x))/Fs;
        end
        
end

end


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
%%%%%%read_dat_file.m%%%%%%%%%%%
% [x,tfs,tfe,fs]=read_dat_file(dirname,fname,fs,tstart,nsec,units_voltage)
% Read in a DAT file.
%  WARNING!  LOG file must be in same location.
% dirname: directory where DAT and LOG file are located.
% fname, name of *DAT file.  MUST include DAT extension!
% fs, sampling frequency in Hz.
% tstart-datenumber of desired start time.  If zero, data are read from
%   start of file
% nsec-number of seconds to be pulled.
% units_voltage: if 1, then output in terms of voltage, otherwise output in terms of uPa.
%       Default acoustic uPa
% Output:
%   tfs: datenumber of the start of the file
%   tfe: datenumber of the end of the file

% Version 1, Sept 30, 2006-read in raw data
% Version 2, Oct. 28 2006- correct for time lost from serial port
% acquisition--for this deployment every five minutes acoustic data not
% acquired for 4 seconds (plus some time to power up).  Thus I use
%   a drift formula of 48 sec/hour..


function [x,tfs,tfe,fs]=read_dat_file(dirname,fname,fs,tstart,ns,units_voltage)

if ~exist('units_voltage'),
    units_voltage=0;
end
%fs=50000;
VREF=2.5;
nc=1;
np=ns*fs;
sens=-182; %dB re 1uPa/V

%Get start time of file from .log file...
mydir=pwd;
cd(dirname);
fpre=findstr(fname,'.DAT')-1;
log_name=[fname(1:fpre) '.log'];
fid=fopen(log_name,'r');

line.fs=[];
tfs=0;tfe=0;
sec_flag=0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if findstr(tline,'Start time')
        line.start=tline;
        tfs=parse_date(line.start);
        sec_flag=1;
    elseif findstr(tline,'Stop time')
        line.end=tline;
        tfe=parse_date(line.end);
        sec_flag=2;
    elseif findstr(tline,'Sample rate:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'Hz')-1;
        fs=str2num(tline(Idot:Iend));
    elseif findstr(tline,'Sensitivity:')>0,
        Idot=1+findstr(tline,':');
        Iend=findstr(tline,'dB')-1;
        sens=str2num(tline(Idot:Iend));
        
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


if tfs==0|tfe==0
    disp(sprintf('Could not understand %s',log_name));
end

try
    if tstart>tfs,
        offset=datevec(tstart-tfs);
        offset=round(offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file
        %disp(sprintf('%i second offset from file start',offset));
    elseif tstart>0
        %error('cannot read this time from this file');
        disp('desired start is before file start, setting offset to 0');
        offset=0;
    else
        disp('Tstart less than zero, exiting gracefully...');
        x=[];
        return
    end
catch
    disp('File name does not contain date, will interpret tstart as elasped time in seconds');
    %keyboard;
    offset=tstart;
end

scale=(0.038/1.0)*VREF/(2^16-1);  %An input voltage of 0.05 V p-p produces a 1.2 V p-p at A/D
fid=fopen(fname,'r','ieee-be');
cd(mydir);
fskip=2*nc*fs*offset;
res=fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
%disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));
x=fread(fid,[nc np],'uint16');
fclose(fid);
x=scale*x;


%Crude acoustic level calibration
if units_voltage~=1,
    disp('Acoustic calibration executed');
    x=x.*(10^(-sens/20));
end
end
%Future calibration corrections

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
sensitivity=-182; %dB re 1uPa/V

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
        offset=round(offset(:,6)+60*offset(:,5)+3600*offset(:,4)+24*3600*offset(:,3));  %offset in seconds from file
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
fskip=2*nc*round(fs*offset);
np=round(fs*ns);
res=fseek(fid,fskip,-1);
if res<0,
    keyboard;
end
disp(sprintf('Offset %10.4f s into %s',ftell(fid)/(2*nc*fs),fname));
x0=fread(fid,[nc np],'uint16');
fclose(fid);
x=scale*(x0-bias);


%Crude acoustic level calibration
if calibrate_acoustic==1,
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
if strcmp(keyword,'short')|~isempty(findstr(keyword,'Shell08')),
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

% --- Executes on button press in pushbutton_GSIbearing.
function pushbutton_GSIbearing_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_GSIbearing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tdate_start=datenum(get(handles.edit_datestr,'String'));
tlen=str2num(get(handles.edit_windowlength,'String'));
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
%Ichan=str2num(get(handles.edit_chan,'String'));
[x,t,Fs,tstart,tend,head]=load_data(handles.filetype,handles.min_time,tdate_start,tlen,'all',handles.mydir,handles.myfile);

disp('Click on two extreme corners:');
tmp=ginput(2);
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

%titlestr=get(handles.text_filename,'String');
%Iend=findstr(titlestr,'.gsi')-1;
set(handles.text_filename,'String',sprintf('%s/%s %6.2f degrees... ',handles.mydir,handles.myfile,thet));
%keyboard;
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

