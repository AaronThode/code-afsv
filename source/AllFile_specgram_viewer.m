
function varargout = AllFile_specgram_viewer(varargin)
%Edited considerably by Jit Sarkar
%	Version history is in associated GIT repository
%
%
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

% Last Modified by GUIDE v2.5 01-Oct-2014 16:39:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AllFile_specgram_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @AllFile_specgram_viewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
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
function AllFile_specgram_viewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<VANUS>
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

if ~ isdeployed
    path('~/Desktop/Ulysses',path);
    path('~/Desktop/Ulysses/deps',path);
end

startupinfo		=	gui_startup_information;

%Make a splash screen
try
    s = SplashScreen( 'Splashscreen', '~/Desktop/Ulysses/104703.JPG', ...
        'ProgressBar', 'on', ...
        'ProgressPosition', 5, ...
        'ProgressRatio', 0.4 );
    s.addText( 30, 50, 'Ulysses/Ribbit', 'FontSize', 30, 'Color', [0 0 0.6] )
    s.addText( 30, 80, 'v1.1, August 26, 2014', 'FontSize', 20, 'Color', [0.2 0.2 0.5] )
    s.addText( 30, 110, 'Aaron Thode and Jit Sarkar', 'FontSize', 20, 'Color', [0.2 0.2 0.5] )
    s.addText( 30, 150, 'Adventure through Analysis', 'FontSize', 20, 'Color', [0.2 0.2 0.5] ,'Fontangle','Italic')
    
    s.addText( 300, 270, 'Loading...', 'FontSize', 20, 'Color', 'white' )
    
    pause(2)
    delete( s )
catch
    disp('Splash Screen failure...');
end

%%%%Set up times based on available times in input directory...%%%%
handles.mydir			=	startupinfo.default_directory;
handles.inputdir		=	startupinfo.default_inputfiledir;
handles.outputdir=uigetdir(handles.mydir,'Select a default output directory for analyses and printouts:');
if ~isempty(handles.outputdir)
    cd(handles.outputdir);
end
handles.annotation_file	=	startupinfo.annotation_file;
try
    handles.calibration_DASAR2007_dir	=	startup_info.calibration_DASAR2007_dir;
catch
    disp('2007 calibration directory not defined');
end
set(handles.text_filename,'String',[handles.mydir ]);

mydir	=	pwd;
try
    cd(handles.mydir);
end
handles.filetype	=	'mt';
[handles,errorflag] =	set_slider_controls(handles,handles.filetype); %#ok<*NASGU>

%Group controls into multielement (arrays) and linked (distributed array)
%categories, for easier times in hiding or disabling buttons

handles.buttongroup.linked(1)=handles.pushbutton_next_linked_annotation;
handles.buttongroup.linked(2)=handles.pushbutton_previous_linked_annotation;
handles.buttongroup.linked(3)=handles.pushbutton_next_station;
handles.buttongroup.linked(4)=handles.pushbutton_prev_station;

handles.buttongroup.array(1)=handles.pushbutton_CSDM;
handles.buttongroup.array(2)=handles.pushbutton_Mode;
handles.buttongroup.array(3)=handles.pushbutton_tilt;
handles.buttongroup.array(4)=handles.pushbutton_modalfiltering;

handles.buttongroup.multichan(1)=handles.togglebutton_ChannelBeam;
handles.buttongroup.multichan(2)=handles.edit_chan;

handles.buttongroup.GSI(1)=handles.pushbutton_GSIbearing;
handles.buttongroup.GSI(2)=handles.pushbutton_GSI_localization;

handles.buttongroup.accelerometer(1)=handles.edit_normal_rotation;
handles.buttongroup.accelerometer(2)=handles.text_normal_rotation;

%Make GSIbearing button, CSDM, and accelerometer buttons invisible
for I=1:length(handles.buttongroup.linked)
    set(handles.buttongroup.linked(I),'Vis','off');
end
%set(handles.pushbutton_next_linked_annotation,'Vis','off');
%set(handles.pushbutton_previous_linked_annotation,'Vis','off');

for I=1:length(handles.buttongroup.GSI)
    set(handles.buttongroup.GSI(I),'Vis','off');
end

for I=1:length(handles.buttongroup.multichan)
    set(handles.buttongroup.multichan(I),'Vis','off');
end
%set(handles.pushbutton_GSIbearing,'Vis','off');
%set(handles.pushbutton_GSI_localization,'Vis','off');

for I=1:length(handles.buttongroup.array)
    set(handles.buttongroup.array(I),'Vis','off');
end

for I=1:length(handles.buttongroup.accelerometer)
    set(handles.buttongroup.accelerometer(I),'Vis','off');
end


%	Set prev/next buttons inactive initially
set(handles.pushbutton_next,'Enable','off');
set(handles.pushbutton_prev,'Enable','off');

set(handles.edit_maxfreq,'String','0');
set(handles.edit_minfreq,'String','0');
cd(mydir);
%%%Create list of function handles to be used in other specialized
%%%applications

%	Check for required dependencies, i.e. subfunctions
status	=	dependency_check();
if ~status
    error('You do not seem to have all the required subfunctions on your path');
end



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

% UIWAIT makes AllFile_specgram_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = AllFile_specgram_viewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


end

%%	GUI menu item callbacks

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

%	List of supported file types + descriptions/names
File_types	=	{'MT';'SIO';'WAV';'GSI';'ADI';'DAT';'MDAT';'PSD';'MAT';'M4V'};
File_descs	=	...
    {	'MT files';...
     'SIO files';...
    'WAV files';...
    'GSI files';...
    'ADI files';...
    'DAT files';...
    'MDAT files';...
    'PSD binary files'; ...
    'Simulations'; ...
    'M4A,M4V (GoPro) files'};

%	 Setup menu file extensions
menustr		=	cell(length(File_types)+1,1);
all_types	=	[];
for	ff = 1:length(File_types)
    ftype	=	File_types{ff};
    fdesc	=	File_descs{ff};
    menustr{ff+1,1}	=	['*.' lower(ftype), ';*.' upper(ftype)];
    menustr{ff+1,2}	=	[fdesc ' (' menustr{ff+1,1} ')'];
    
    all_types	=	[all_types ';' menustr{ff+1,1}];
end
%	Prepend all supported files
all_types(1)	=	[];			%	Get rid of initial ;
menustr{1,1}	=	all_types;
menustr{1,2}	=	['All supported files' ' (...)'];


%	Specify default file
default_file	=	[];
if ~isempty(handles.mydir)
    disp(handles.mydir);
    default_file	=	handles.mydir;
    try
        default_file	=	fullfile(default_file, handles.myfile);
    end
end

%	Open dialog box
[filename, pathname, ~]	=	...
    uigetfile(menustr, 'Select a file:', default_file);

%	Parse selected file name
if	isnumeric(filename) || isnumeric(pathname)
    return;		%	window canceled
end

[~,fname,myext]		=	fileparts(filename);

if	~isempty(myext)
    myext	=	myext(2:end);
    file_ind=	find(strcmpi(myext,File_types));
end

if	isempty(myext) || isempty(file_ind)
    disp([myext ' File type not supported']);
    return;
end

%	Assume valid data type at this point
handles.mydir		=	pathname;
handles.myfile		=	filename;
handles.myext		=	myext;
handles.filetype	=	File_types{file_ind};
handles.filedesc	=	File_descs{file_ind};

%yes_wav=get(handles.togglebutton_getwav,'value');

[handles,~]		=	set_slider_controls(handles, handles.filetype);
set(handles.text_filename,'String',fullfile(handles.mydir, handles.myfile));
set(handles.text_filetype,'String',handles.filetype);

%	Disable buttons that require an initial update
toggle_initial_buttons(handles, 'off');
set(handles.pushbutton_update, 'Enable', 'on');
handles.fig_updated		=	false;

%	Enable notes folder selection, in case it's not already enabled
set(handles.pushbutton_notes_select, 'Enable', 'on');
%	set notes file name and dir accordingly
if	isempty(handles.notes.folder_name)
    notes_folder	=	pathname;
else
    notes_folder	=	[];
end

%	Jit: Loading should always occur if valid files are present in given
%	directory
%	Otherwise same folder has to be reselected everytime new data file in
%	same series is loaded
handles		=	load_notes_file(handles, notes_folder);  %OpenMenu Callback



% Update handles structure
guidata(hObject, handles);
end

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
set(handles.pushbutton_fileread, 'Enable', 'On');
cd(mydir);
guidata(hObject, handles);

end

% --------------------------------------------------------------------
function Menu_processors_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_processors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function MenuItem_boatdet_Callback(hObject, eventdata, handles)
% hObject    handle to MenuItem_boatdet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%[metrics,hf]	=	shipdet_metrics(PSD_struct, param, ha_ref, hf,want_hough)
%PSD_struct	:	{struct}	   Power spectral density structure,
%								must contain the following	fields:
%		PSD		:	{dbl}(t,f)	The 2d power spectral density matrix, linear (not dB) units
%		T		:	{dbl}(t)	Time vector
%		F		:	{dbl}(f)	Frequency vector
%       Idebug:                 scalar that determines if debug output given...
%
persistent Batch_vars Batch_desc resett

Batch_type	=	get(hObject, 'Label');
Batch_mode	=	input_batchmode(Batch_type);

if	isempty(Batch_mode)
    errordlg('Processing cancelled!');
    return;
end

if isempty(resett)
    resett=0;  %Allows context-dependent Batch_vars, based on Batch_mode selection
end

if isempty(Batch_vars)||resett  %Default
    Batch_vars.threshold	=	'5';	Batch_desc{1}	=	'dB threshold between peak and surrounding values at +/-df_search';
    Batch_vars.df_search	=	'0.2';	Batch_desc{2}	=	'maximum half-bandwidth (kHz) to examine around each local maximum';
    Batch_vars.f_min	=	get(handles.edit_fmin,'string');	Batch_desc{3}	=	'minimum frequency to examine (kHz)';
    Batch_vars.f_max	=	get(handles.edit_fmax,'string');	Batch_desc{4}	=	'maximum frequency to examine (kHz)';
    Batch_vars.sec_avg	=	'2';	Batch_desc{5}	=	'Averaging time of spectrogram (sec)--ignored when bulk loading';
    resett=0;
end

if ~isempty(strfind(Batch_mode,'Start Bulk Processing'))
    clear Batch_vars Batch_desc
    Batch_vars.sec_avg	=	'2';	Batch_desc{1}	=	'Averaging time of spectrogram (sec)';
    resett=1;
    
end
Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, Batch_type);


if	isempty(Batch_vars)
    uiwait(errordlg('Vessel Processing cancelled!'));
    return;
end

try
    sec_avg=str2double(Batch_vars.sec_avg);
    
    if isfield(Batch_vars,'threshold')
        threshold=str2double(Batch_vars.threshold);
        df_search=1000*str2double(Batch_vars.df_search);
        f_min=1000*str2double(Batch_vars.f_min);
        f_max=1000*str2double(Batch_vars.f_max);
        if f_min>f_max
            return
        end
    end
catch
    uiwait(errordlg('Vessel Processing cancelled! Bad batch parameters'));
    return
end

close_all_figures;
switch	Batch_mode
    case 'Process Visible Window'
        uiwait(msgbox(['Processing all data within window for ' Batch_type]));
        yes=1;
        while yes
            dT=handles.sgram.T(2)-handles.sgram.T(1);
            Ncol=max([1 floor(sec_avg/dT)]);
            Itime=1:Ncol:length(handles.sgram.T);
            Tnew=handles.sgram.T(Itime);
            Tnew=Tnew(1:(end-1));
            F=handles.sgram.F;
            PSD=zeros(length(F),length(Itime)-1);
            for I=1:(length(Itime)-1)
                PSD(:,I)=mean(handles.sgram.B(:,Itime(I):Itime(I+1)),2);
            end
            PSD_dB=10*log10(abs(PSD));
            
            
            Nt	=	length(Tnew);
            Nf	=	length(F);
            dF	=	mean(diff(F));
            
            
            PSD_avg		=	sum(PSD,2)/Nt;
            PSD_med     =   median(PSD')';
            P_tot		=	sum(PSD_avg)*dF;
            
            %PSD_dB		=	abs(10*log10(abs(PSD)));
            PfdB_avg	=	sum(PSD_dB,2)/Nt;
            PfdB_CV      =   std(PSD_dB')'./PfdB_avg;
            PfdB_skew    =   skewness(PSD_dB');
            PfdB_kurt   =   kurtosis(PSD_dB');
            
            %Look for mean second derivative.
            %PfdB_delta2=   diff(PSD_dB',2)./((F(2)-F(1)).^2);
            
            %Plot metrics
            figure(1);subplot(1,4,1)
            plot(PfdB_avg,F,'k',10*log10(PSD_avg)-10,F,'r');grid on;hold on
            xlabel('mean dB PSD re 1 uPa^2/Hz','fontsize',14,'fontweight','bold');
            ylabel('Frequency (Hz)','fontsize',14,'fontweight','bold');
            [peakss]=peak_picker_Thode(PfdB_avg,F,df_search,[f_min f_max],threshold,[]);
            if ~isempty(peakss{1}.adp.F)
                plot( peakss{1}.adp.PdB,peakss{1}.adp.F,'bx','markersize',10);
            end
            if ~isempty(peakss{1}.isi.F)
                plot( peakss{1}.isi.PdB,peakss{1}.isi.F,'go','markersize',10);
            end
            legend('averaged dB power','dB of average power-10','adaptive peaks','isi peaks')
            ylim([f_min f_max]);grid on
            xlim(str2num(get(handles.edit_mindB,'String'))+[0 str2num(get(handles.edit_dBspread,'string'))]);
            
            %Plot coefficient of variation of band over time period in question
            subplot(1,4,2);
            plot(PfdB_CV,F,'k');xlabel('C.V. of dB power','fontsize',14,'fontweight','bold');
            ylabel('Frequency (Hz)','fontsize',14,'fontweight','bold');
            ylim([f_min f_max]);
            hold on
            plot_peak_values(peakss{1}.adp.F,F,PfdB_CV,'ro');
            plot_peak_values(peakss{1}.isi.F,F,PfdB_CV,'bo');
            grid on
            
            subplot(1,4,3);
            plot(PfdB_skew,F,'k');xlabel('skewness ','fontsize',14,'fontweight','bold');
            ylabel('Frequency (Hz)','fontsize',14,'fontweight','bold');
            ylim([f_min f_max]);
            hold on
            plot_peak_values(peakss{1}.adp.F,F,PfdB_skew,'ro');
            plot_peak_values(peakss{1}.isi.F,F,PfdB_skew,'bo');
            grid on;
            
            subplot(1,4,4);
            plot(PfdB_kurt,F,'k');xlabel('kurtosis ','fontsize',14,'fontweight','bold');
            ylabel('Frequency (Hz)','fontsize',14,'fontweight','bold');
            ylim([f_min f_max]);
            hold on
            plot_peak_values(peakss{1}.adp.F,F,PfdB_kurt,'ro');
            plot_peak_values(peakss{1}.isi.F,F,PfdB_kurt,'bo');
            grid on;
            
            yes=menu('Redo?','Yes','No');
            if yes==1
                Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, Batch_type);
                if	isempty(Batch_vars)
                    uiwait(errordlg('Vessel Processing cancelled!'));
                    return;
                end
                
                try
                    sec_avg=str2double(Batch_vars.sec_avg);
                    threshold=str2double(Batch_vars.threshold);
                    df_search=1000*str2double(Batch_vars.df_search);
                    f_min=1000*str2double(Batch_vars.f_min);
                    f_max=1000*str2double(Batch_vars.f_max);
                catch
                    uiwait(errordlg('Vessel Processing cancelled! Bad batch parameters'));
                    return
                end
            else
                yes=0;
            end
            hold off
        end %while yes==1
        
        %Plot bands on spectrogram...
        if ~isempty(peakss{1}.adp.F)
            FF=peakss{1}.adp.F/1000;
            axes(handles.axes1);
            for II=1:length(FF)
                line(handles.sgram.T([1 end]),FF(II)*[1 1],'color','k','linewidth',3);
            end
        end
        
        if ~isempty(peakss{1}.isi.F)
            FF=peakss{1}.isi.F/1000;
            axes(handles.axes1);
            for II=1:length(FF)
                line(handles.sgram.T([1 end]),FF(II)*[1 1],'color','g','linewidth',3);
            end
        end
    case {'Start Bulk Processing File','Start Bulk Processing Folder'}
        
        start_bulk_processing_psd(handles,Batch_mode,Batch_type,sec_avg);
    case 'Load Bulk Processing'
        Batch_vars_bulkload.start_time	=	'folder start';
        
        Batch_desc_bulkload{1}	=	'Time to begin loading (e.g. "here" to use visible start time, "datenum(2011,1,2,0,0,0)", "file start" for current file , "select" to choose file from list, "folder start" for first file in folder)';
        Batch_vars_bulkload.end_time	=	'folder end';
        Batch_desc_bulkload{2}	=	'Time to end loading (e.g. "here" to use GUI end time, "2" for two hours from start time, "datenum(2011,1,2,0,0,0)", "file end" for time at current file end, "select" to choose file from list, "folder end" to import through last file)';
        Batch_vars_bulkload.time_window      =   'file';
        Batch_desc_bulkload{3} =  'hours to display in a single window; "all" means put in a single file, "file" displays one file per figure';
        Batch_vars_bulkload.xlabel_style = 'auto';
        Batch_desc_bulkload{4} =  'Datenumber format to use when plotting x-axis; "auto" will try to decide for you; examples include "mm/dd", "dd", "dd/HH", "HH" or "HH:MM":';
        Batch_vars_bulkload.stationary_time = get(handles.edit_winlen,'string')
        Batch_desc_bulkload{5} =  'Time over which signal is assumed stationary (sec); for fast boat passages this should be lowered:';
        
        Batch_vars_bulkload	=	input_batchparams(Batch_vars_bulkload, Batch_desc_bulkload, Batch_type);
        if	isempty(Batch_vars_bulkload)
            uiwait(errordlg('Processing cancelled!'));
            return;
        end
        
        %%Set stationary averaging (or median) time
        tabs_inc_stationary=datenum(0,0,0,0,0,str2num(Batch_vars_bulkload.stationary_time));
        
        %%%%%%%%Set parameters for plotting images of bulk data
        %   split_windows =1 will display a separate figure for each PSD
        %   file
        split_windows=1;plot_interval=Inf;val=-1;
        
        if ~isempty(strfind(Batch_vars_bulkload.time_window,'file'))
            plot_interval=Inf;
        elseif ~isempty(strfind(Batch_vars_bulkload.time_window,'all'))
            split_windows=0;
        elseif isnumeric(str2num(Batch_vars_bulkload.time_window))
            val=str2num(Batch_vars_bulkload.time_window);
            plot_interval=datenum(0,0,0,val,0,0);
        end
        
        date_tick_chc=get_datetick_style(Batch_vars_bulkload.xlabel_style,val,'hours');
        
        %Load bulk run time and file data
        bulk_params=load_PSD_Bulk_Run(handles, Batch_vars_bulkload);
        
        tabs_folder_start=bulk_params.tabs_folder_start;
        tabs_folder_end=bulk_params.tabs_folder_end;
        tabs_start=bulk_params.tabs_start;
        tabs_end=bulk_params.tabs_end;
        Other_FileNames=bulk_params.Other_FileNames;
        Icurrent_file=bulk_params.Icurrent_file;
        Ifinal_file=bulk_params.Ifinal_file;
        FF=bulk_params.FF;
        
        %Start processing
        tabs_loop_begin=tabs_start;
        boat_detections=[];Tabs_all=[];
        Iplot=1;
        for I=Icurrent_file:Ifinal_file
            fname=Other_FileNames(I).name;
            fprintf('Processing %s...\n',fname);
            [PSD,F,Tsec,Tabs,params]=read_Java_PSD(fname,tabs_loop_begin,Inf);  %Read an entire file into memory
            
            sec_avg=(params.Nfft+params.dn*params.Nsamps)/params.Fs;  %Reconstruct averaging time used in these files.
            
            Istrip=find(Tabs<=tabs_end);
            PSD=PSD(:,Istrip);
            Tabs=Tabs(Istrip);
            PSD_dB=10*log10(abs(PSD));
            
            %Round start down to nearest minute
            tmp=datevec(Tabs(1));
            tmp(end)=0;
            Tabs(1)=datenum(tmp);
            
            %Create time intervals over which to create peaks
            Tabs_stationary=Tabs(1):tabs_inc_stationary:Tabs(end);
            Npeaks=length(Tabs_stationary)-1;
            Maxbands=10;
            
            %Create peaks for entire run
            for Is=1:(length(Tabs_stationary)-1)
                
                Ipass=find(Tabs>=Tabs_stationary(Is)&Tabs<=Tabs_stationary(Is+1));
                PfdB_avg	=	sum(PSD_dB(:,Ipass),2)/length(Ipass);
                peakss=peak_picker_Thode(PfdB_avg,F,df_search,[f_min f_max],threshold,[]);
                boat_detections=merge_structure(boat_detections,peakss{1}.adp,Maxbands,Npeaks);
                
                if isempty(boat_detections)
                    Tabs_all=median(Tabs(Ipass));
                else
                    Tabs_all=[Tabs_all median(Tabs(Ipass))];
                end
                
                
            end  %Is
            
            
            if split_windows
                
                [hprint(Iplot),save_tag{Iplot}]=plot_boat_detections(boat_detections,Tabs_all);
                if Ifinal_file-Icurrent_file>10
                    save_str=sprintf('Boat_%s', save_tag{Iplot});
                    save(save_str,'boat_detections','Tabs_all','params','titlestr','Batch_vars','Batch_vars_bulkload','sec_avg');
                    %uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
                    
                    orient tall
                    print(hprint(II),'-djpeg','-r300',save_str);
                    saveas(hprint(II),save_str,'fig');
                    close(hprint(II));
                    
                end
                Tabs_all=[]; boat_detections=[];
                Iplot=Iplot+1;
            end
            
            %If we move into additional files, start at the beginning of
            %the file...
            tabs_loop_begin=0;
            
        end  %Icurrent_file
        
        if ~split_windows
            [hprint(Iplot),save_tag{Iplot}]=plot_boat_detections(boat_detections,Tabs_all);
            Iplot=Iplot+1;
        end
        
        if Ifinal_file-Icurrent_file>10
            
            yess=menu('Save boat detection Data and figure?','Yes','No');
            if yess==1
                
                Nfigs=sort(get(0,'Child'));
                Nfigs=Nfigs(Nfigs<150);
                
                for II=1:length(Nfigs)
                    
                    
                    save_str=sprintf('Boat_%s', save_tag{II});
                    save(save_str,'boat_detections','Tabs_all','params','titlestr','Batch_vars','Batch_vars_bulkload','sec_avg');
                    uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
                    
                    orient tall
                    print(hprint(II),'-djpeg','-r300',save_str);
                    saveas(hprint(II),save_str,'fig');
                    if length(Nfigs)>10
                        close(hprint(II));
                    end
                    
                end %II
                
            end %%yes
        end %If Ifinal_file
end

%%%%Inner function for plotting statistical values associated with PSD
%%%%peaks
    function plot_peak_values(Fpeaks,F,PP,plotstr)
        %Plot peak detections from isi and adaptive band detectors
        if ~isempty(Fpeaks)
            
            for JJ=1:length(Fpeaks)
                [~,Ipeaks]=min(abs(Fpeaks(JJ)-F));
                plot( PP(Ipeaks),Fpeaks(JJ),plotstr,'markersize',10);
            end
        end
    end

%%Inner function for plotting bulk results

    function [hprint,save_tag]=plot_boat_detections(boat_detections,Tabs_all)
        
        if isempty(Tabs_all)
            hprint=-1;save_tag=-1;
            return
        end
        
        twin1=min(Tabs_all);
        twin2=max(Tabs_all);
        save_tag=[datestr(twin1,30) '_' datestr(twin2,30)];
        
        hprint=figure;
        
        subplot(2,1,1)
        plot(Tabs_all,boat_detections.Npeaks,'ko');grid on
        [AX,H1,H2]=plotyy(Tabs_all,boat_detections.Npeaks,Tabs_all,boat_detections.pdBinttotal(1,:));
        %plot(Tabs_all,boat_detections.Npeaks);
        ylabel(AX(1),'Bands detected');
        ylabel(AX(2),'Int. power, dB re 1uPa');
        
        pmax=ceil(max(boat_detections.pdBinttotal(1,:)));
        pmin=floor(min(boat_detections.pdBinttotal(1,:)));
        ylimm=get(AX(2),'ylim');ylimm(1)=pmin;ylimm(2)=pmax;
        ylim(AX(2),ylimm);
        set(AX(2),'ytick',round(linspace(0,pmax,4)));
        if isempty(date_tick_chc) %auto adjustment has failed..
            date_tick_chc=get_datetick_style('auto',twin2-twin1,'datenumber');
        end
        
        for IJ=1:length(AX)
            xlim(AX(IJ),[twin1 twin2]);
            datetick(AX(IJ),'x',date_tick_chc,'keeplimits');
            set(AX(IJ),'fontweight','bold','fontsize',14);
            grid(AX(IJ),'on');
        end
        
        %sec_avg=params.Nsamps*params.dn/params.Fs;
        titlestr=(sprintf('Start time: %s, End Time: %s, seconds averaged: %6.2f, stationary time: %s sec', ...
            datestr(twin1,30),datestr(twin2,30),sec_avg,Batch_vars_bulkload.stationary_time));
        title(titlestr);
        
        
        subplot(2,1,2)
        colorstr='bgrcmyk';
        for J=1:boat_detections.maxbands
            ff=boat_detections.F(J,:);
            bw=boat_detections.bandwidth(J,:)/2;
            
            XX=[1;1]*Tabs_all;
            YY=[ff+bw;ff-bw];
            Icolor=max([1 rem(J,length(colorstr))]);
            hh=line(XX,YY,'color',colorstr(Icolor),'linewidth',3);
            hold on
        end
        xlim([twin1 twin2]);
        datetick('x',date_tick_chc,'keeplimits');
        set(gca,'fontweight','bold','fontsize',14);
        grid on
        ylabel('Hz');
        title(sprintf('%s-%s kHz bands scanned, adaptive bandwidth %s Hz, threshold %s dB', ...
            Batch_vars.f_min,Batch_vars.f_max,Batch_vars.df_search,Batch_vars.threshold))
        %  threshold=str2double(Batch_vars.threshold);
        %df_search=str2double(Batch_vars.df_search);
        
        ylim([0 str2num(get(handles.edit_fmax,'string'))*1000]);
        
    end
end

function start_bulk_processing_psd(handles,Batch_mode,Batch_type,sec_avg)
if strcmp(Batch_mode,'Start Bulk Processing File')
    uiwait(msgbox(sprintf('Processing all data in file %s\n for %s',fullfile(handles.mydir,handles.myfile), Batch_type)));
    param.exten=['*' handles.myfile];
else
    uiwait(msgbox(sprintf('Processing all data in folder %s\n for %s',fullfile(handles.mydir), Batch_type)));
    
    param.exten=['*' handles.myext];
end

data_folder_name		=	uigetdir('.', 'Select a folder to store analysis results:');
param.dir_out=data_folder_name;
param.file_dir=handles.mydir;
contents=get(handles.popupmenu_Nfft,'String');
param.Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});
contents=get(handles.popupmenu_ovlap,'String');
param.ovlap=round(param.Nfft*str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100);
param.Fs=handles.sgram.Fs;
param.channel=str2num(get(handles.edit_chan,'String'));
param.dumpsize=100000;
param.nstart=0;
param.nsamples=0;
param.f_low=1000*str2num(get(handles.edit_fmin,'String'));
param.f_high=1000*str2num(get(handles.edit_fmax,'String'));
if sec_avg>0
    param.sec_avg=sec_avg;
else
    param.sec_avg=param.Nfft/param.Fs;  %Make a spectrogram--no averaging
end

cd(param.dir_out);
write_Java_script('PSD',param);
! ./masterPSD.scr > outt.txt &
return
end

% --------------------------------------------------------------------
function MenuItem_psd_Callback(hObject, eventdata, handles)
% hObject    handle to MenuItem_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Batch_type	=	get(hObject, 'Label');
Batch_mode	=	input_batchmode(Batch_type);

if	isempty(Batch_mode)
    uiwait(errordlg('Processing cancelled!'));
    return;
end

if isempty(strfind(Batch_mode,'Load'))
    Batch_vars.sec_avg	=	'2';	Batch_desc{1}	=	'Seconds to average PSD for long-term display, if "0" no averaging' ;
    Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, Batch_type);
    if	isempty(Batch_vars)
        uiwait(errordlg('Processing cancelled!'));
        return;
    end
    
    sec_avg=str2double(Batch_vars.sec_avg);
    
end

%%%%%
close_all_figures
switch	Batch_mode
    case 'Process Visible Window'
        uiwait(msgbox(['Processing all data within window for ' Batch_type]));
        
        dT=handles.sgram.T(2)-handles.sgram.T(1);
        Ncol=max([1 floor(sec_avg/dT)]);
        Itime=1:Ncol:length(handles.sgram.T);
        Tnew=handles.sgram.T(Itime);
        Tnew=Tnew(1:(end-1));
        F=handles.sgram.F;
        PSD=zeros(length(F),length(Itime)-1);
        for I=1:(length(Itime)-1)
            PSD(:,I)=mean(handles.sgram.B(:,Itime(I):Itime(I+1)),2);
        end
        PSDdB=10*log10(PSD);
        hprint=figure;
        
        imagesc(Tnew,F/1000,PSDdB);axis('xy')
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (sec)');ylabel(' Frequency (kHz)');
        
        titlestr=(sprintf('PSD Start time: %s, Time processed: %6.2f sec',datestr(handles.tdate_start),handles.tlen));
        title(titlestr);
        
        
        fmin=str2num(get(handles.edit_fmin,'String')); %kHz
        fmax=str2num(get(handles.edit_fmax,'String'));  %kHz
        ylimm=[fmin fmax];
        ylim(ylimm);
        colorbar
        cmin=eval(get(handles.edit_mindB,'string'));
        cmax=cmin+eval(get(handles.edit_dBspread,'string'));
        caxis([cmin cmax]);
        
        %yess=menu('Save PSD Data and figure?','Yes','No');
        ButtonName = questdlg('What do you want to do now?', ...
            'PSD options', ...
            'Leave Me Alone', 'Save Figure', 'Plot percentiles', 'Leave Me Alone');
        switch ButtonName
            case 'Save Figure'
                
                save_str=sprintf('PSD_plot_%4.2f_%s_%s', sec_avg,datestr(handles.tdate_start,30), ...
                    datestr(handles.tdate_start+datenum(0,0,0,0,0,handles.tlen),30));
                save(save_str,'F','titlestr','Tnew','PSD','cmin','cmax','ylimm','sec_avg');
                uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
                
                orient landscape
                print(hprint,'-djpeg','-r300',save_str);
            case 'Plot percentiles'
                
                pms=get_PSD_percentile_params(fmin,fmax);
                pms.label_style=get_datetick_style([]);
                
                pms.y_label='dB re 1uPa';
                Tabs=datenum(get(handles.edit_datestr,'String'))+datenum(0,0,0,0,0,Tnew);
                 
                %%Process PSD matrix..
                for Iff=1:length(pms.fmin)
                    Igood= (F>=pms.fmin(Iff)&F<=pms.fmax(Iff));
                    sumPSD=10*log10((F(2)-F(1))*sum(PSD(Igood,:),1));  %Now power spectral density converted to power
                    pms.title=sprintf('Spectral power between %i and %i Hz, beginning %s',pms.fmin(Iff),pms.fmax(Iff),datestr(Tabs(1)));
                    [temp,hprint]=create_percentile_distributions(Tabs, sumPSD,pms);  %make whisker plot
                    
                    if Iff==1
                        percents.data=zeros(length(pms.fmin),length(pms.percentiles),size(temp.data,2));
                        percents.x=temp.x;
                        percents.percentiles=pms.percentiles;
                        
                        
                    end
                    percents.data(Iff,:,:)=temp.data;
                    percents.title{Iff}=pms.title;
                end
                
                %%Plot contour plots if enough frequency bands
                if pms.plot.image
                    for Ipp=1:length(pms.percentiles)
                        figure
                        
                        yy=squeeze(percents.data(:,Ipp,:));
                        Igood=~all(isnan(yy));
                        imagesc(percents.x(Igood),pms.fmin,yy(:,Igood));
                        axis('xy');
                        datetick('x',pms.label_style);
                        colorbar
                        caxis(pms.y_limits);
                        title(sprintf('%ith percentile, %6.2f s average, starting %s',100*pms.percentiles(Ipp),sec_avg,datestr(Tabs(1))));
                        xlabel('Time','fontweight','bold','fontsize',14);
                        ylabel('Frequency (Hz)','fontweight','bold','fontsize',14);
                    end
                    
                end
                
               
                
        end
        
        return
    case {'Start Bulk Processing File','Start Bulk Processing Folder'}
        
        start_bulk_processing_psd(handles,Batch_mode,Batch_type,sec_avg);
        
        
    case 'Load Bulk Processing'
        close_all_figures;
        Batch_vars_bulkload.start_time	=	'select';
        
        Batch_desc{1}	=	'Time to begin loading (e.g. "here" to use visible start time, "datenum(2011,1,2,0,0,0)", "file start" for current file , "select" to choose file from list, "folder start" for first file in folder)';
        Batch_vars_bulkload.end_time	=	'select';
        Batch_desc{2}	=	'Time to end loading (e.g. "here" to use GUI end time, "2" for two hours from start time, "datenum(2011,1,2,0,0,0)", "file end" for time at current file end, "select" to choose file from list, "folder end" to import through last file)';
        Batch_vars_bulkload.time_window      =   'file';
        Batch_desc{3} =  'hours to display in a single window; "all" means put in a single file, "file" displays one file per figure';
        Batch_vars_bulkload.xlabel_style = 'auto';
        Batch_desc{4} =  'Datenumber format to use when plotting x-axis; "auto" will try to decide for you; examples include "mm/dd", "dd", "dd/HH", "HH" or "HH:MM":';
        
        Batch_vars_bulkload	=	input_batchparams(Batch_vars_bulkload, Batch_desc, Batch_type);
        if	isempty(Batch_vars_bulkload)
            uiwait(errordlg('Processing cancelled!'));
            return;
        end
        
        
        
        %%%%%%%%Set parameters for plotting images of bulk data
        %   split_windows =1 will display a separate figure for each PSD
        %   file
        split_windows=1;plot_interval=Inf;val=-1;
        
        if ~isempty(strfind(Batch_vars_bulkload.time_window,'file'))
            plot_interval=Inf;
        elseif ~isempty(strfind(Batch_vars_bulkload.time_window,'all'))
            split_windows=0;
        elseif isnumeric(str2num(Batch_vars_bulkload.time_window))
            val=str2num(Batch_vars_bulkload.time_window);
            plot_interval=datenum(0,0,0,val,0,0);
        end
        
        date_tick_chc=get_datetick_style(Batch_vars_bulkload.xlabel_style,val,'hours');
        
        %Load bulk run time and file data
        bulk_params=load_PSD_Bulk_Run(handles, Batch_vars_bulkload);
        
        tabs_folder_start=bulk_params.tabs_folder_start;
        tabs_folder_end=bulk_params.tabs_folder_end;
        tabs_start=bulk_params.tabs_start;
        tabs_end=bulk_params.tabs_end;
        Other_FileNames=bulk_params.Other_FileNames;
        Icurrent_file=bulk_params.Icurrent_file;
        Ifinal_file=bulk_params.Ifinal_file;
        FF=bulk_params.FF;
        
        
        %%%%%%%%%Final actions for PSD output displays%%%%%%%%
        PSDButtonName = questdlg('What do you want to do now?', ...
            'PSD options', ...
            'Display and Print Figure', 'Plot percentiles', 'Display and Print Figure');
        
        
        if strcmp(PSDButtonName,'Plot percentiles')
            pms=get_PSD_percentile_params(min(FF)/1000,max(FF)/1000);
            pms.y_label='dB re 1uPa';
            Nf=length(pms.fmin);
            %date_tick_chc=get_datetick_style(date_tick_chc,pms.xlabel_inc,'datenumber');
            
        end
        
        
        %Start processing
        tabs_loop_begin=tabs_start;
        Nfigs_max=10;
        PSD_all=[];Tabs_all=[];
        Iplot=0;
        Nfigs=length(Icurrent_file:Ifinal_file);
        for I=Icurrent_file:Ifinal_file
            fname=Other_FileNames(I).name;
            fprintf('Processing %s...\n',fname);
            [PSD,F,Tsec,Tabs,params]=read_Java_PSD(fname,tabs_loop_begin,Inf);  %Read an entire file into memory
            
            %If we load more than one  file, start at the beginning
            tabs_loop_begin=0;
            
            sec_avg=(params.Nfft+params.dn*params.Nsamps)/params.Fs;
            
            Istrip=find(Tabs<=tabs_end);
            switch PSDButtonName
                case 'Display and Print Figure'
                    PSD_all=[PSD_all PSD(:,Istrip)];
                    Tabs_all=[Tabs_all Tabs(Istrip)];
                    
                    if ~split_windows
                        continue
                    end
                    Iplot=Iplot+1;
                    [~,local_fname{Iplot},~]=fileparts(fname);
                    
                    [hprint(Iplot),save_tag{Iplot}]=image_PSD(min(Tabs_all), max(Tabs_all),local_fname{Iplot});
                    
                    save_str=sprintf('PSD_%s', save_tag{Iplot});
                    if ~exist('pms'),pms=[];end
                    save(save_str,'pms','PSD_all','Tabs_all','params','titlestr','Batch_vars_bulkload','sec_avg');
                    %uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
                    disp([save_str ' mat and jpg file written to ' pwd '\n']);
                    
                    if Nfigs>Nfigs_max
                        
                        orient landscape
                        figure(hprint(Iplot));
                        
                        print(hprint(Iplot),'-djpeg','-r300',save_str);
                        
                        %Plot zooms if desired
                        Tabs_limm=get(gca,'xlim');
                        Tabs_frame=unique([Tabs_limm(1):plot_interval:Tabs_limm(2) Tabs_limm(2)]);
                        Tabs_frame(isnan(Tabs_frame))=[];
                        
                        
                        for JJ=1:(length(Tabs_frame)-1)
                            try
                                xlim([Tabs_frame(JJ) Tabs_frame(JJ+1)]);
                                datetick('x',date_tick_chc,'keeplimits');
                                save_str_zoom=sprintf('PSDzoom_%s_%s', datestr(Tabs_frame(JJ),30), datestr(Tabs_frame(JJ+1),30));
                                
                                titlestr=sprintf('%s\nStart time: %s, End Time: %s, seconds averaged: %6.2f', ...
                                    local_fname{Iplot},datestr(Tabs_frame(JJ),30),datestr(Tabs_frame(JJ+1),30),sec_avg);
                                title(titlestr);
                                print(hprint(Iplot),'-djpeg','-r300',save_str_zoom);
                            catch
                                disp('Failure to plot Tabs_frame: plot_interval is likely bad');
                            end
                            
                        end
                        close;
                        
                    end  %if Nfigs>Nfigs_max
                    PSD_all=[];Tabs_all=[];
                    
                    
                case 'Plot percentiles'
                    
                    %%Process PSD matrix..
                    for Iff=1:length(pms.fmin)
                        Igood= (F>=pms.fmin(Iff)&F<=pms.fmax(Iff));
                        sumPSD=10*log10((F(2)-F(1))*sum(PSD(Igood,Istrip),1));  %Now power spectral density converted to power
                        
                        %%Reserve memory for rest of run..
                        if Iff==1&&I==Icurrent_file
                            PSD_all=zeros(Nf,(Ifinal_file-Icurrent_file+1)*length(sumPSD));
                            Tabs_all=zeros(1,(Ifinal_file-Icurrent_file+1)*length(sumPSD));
                            Icount=1;
                        end
                        PSD_all(Iff,Icount:(Icount+length(sumPSD)-1))=sumPSD;
                        
                    end %pms.fmin
                    
                    Tabs_all(Icount:(Icount+length(Istrip)-1))=Tabs(Istrip);
                    Icount=Icount+length(Istrip);
            end  %switch
            
            
        end  %Icurrent_file
        
        
        %%All files have been loaded and processed.
        
        switch PSDButtonName
            case 'Display and Print Figure'
                
                
                if ~split_windows
                    [hprint(1),save_tag{1}]=image_PSD(min(Tabs_all), max(Tabs_all));
                    save_str=sprintf('PSD_%s', save_tag{1});
                    if ~exist('pms'),pms=[];end
                    save(save_str,'pms','PSD_all','Tabs_all','params','titlestr','Batch_vars_bulkload','sec_avg');
                    disp([save_str ' mat and jpg file written to ' pwd '\n']);
                    
                end
                
                if Nfigs>Nfigs_max
                    return
                end
                yess=menu('Save PSD figures?','Yes','No');
                if yess==1
                    
                    Nfigs=sort(get(0,'Child'));
                    Nfigs=Nfigs(Nfigs<150);
                    
                    for II=1:length(Nfigs)
                        
                        save_str=sprintf('PSD_%s', save_tag{II});
                        
                        %uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
                        orient landscape
                        print(hprint(II),'-djpeg','-r300',save_str);
                        
                        %Plot zooms if desired
                        figure(hprint(II));
                        Tabs_limm=get(gca,'xlim');
                        Tabs_frame=unique([Tabs_limm(1):plot_interval:Tabs_limm(2) Tabs_limm(2)]);
                        Tabs_frame(isnan(Tabs_frame))=[];
                        
                        
                        for JJ=1:(length(Tabs_frame)-1)
                            try
                                xlim([Tabs_frame(JJ) Tabs_frame(JJ+1)]);
                                datetick('x',date_tick_chc,'keeplimits');
                                save_str_zoom=sprintf('PSDzoom_%s_%s', datestr(Tabs_frame(JJ),30), datestr(Tabs_frame(JJ+1),30));
                                %[~,local_fname,~]=fileparts(fname);
                                titlestr=sprintf('%s\nStart time: %s, End Time: %s, seconds averaged: %6.2f', ...
                                    local_fname{II},datestr(Tabs_frame(JJ),30),datestr(Tabs_frame(JJ+1),30),sec_avg);
                                title(titlestr);
                                print(hprint(II),'-djpeg','-r300',save_str_zoom);
                            catch
                                disp('Failure to plot Tabs_frame: plot_interval is likely bad');
                            end
                            
                        end
                        
                    end %II
                    
                end %%yes
                %%%%%%Plot and save percentile files to file
            case 'Plot percentiles'
                %Remove excess storage
                %Icount=Icount-length(Istrip);
                if Icount<length(Tabs_all)
                    Tabs_all=Tabs_all(1:Icount);
                    PSD_all=PSD_all(:,1:Icount);
                    
                end
                
                Igood=find(Tabs_all>0);
                Tabs_all=Tabs_all(Igood);
                PSD_all=PSD_all(:,Igood);
                
                
                pms.y_label='dB re 1uPa';
                yes=1;
                while yes
                    close_all_figures;  %Close all open figures..
                    clear hprint
                    
                    
                    %%Process PSD matrix..
                    for Iff=1:length(pms.fmin)
                        tmp=datevec(pms.x_inc);
                        tmp=3600*tmp(4)+60*tmp(5)+tmp(6);
                        save_tag{Iff}=sprintf('Percentile_PSD_%s_%s_fmin%iHz_fmax%iHz_interval%isec', ...
                            datestr(Tabs_all(1),30),datestr(Tabs_all(end),30),(pms.fmin(Iff)),(pms.fmax(Iff)),tmp);
                        sumPSD=PSD_all(Iff,:);  %Now power spectral density converted to power
                        pms.title=sprintf('Spectral power between %i and %i Hz, averaged %6.2f sec, beginning %s', ...
                            pms.fmin(Iff),pms.fmax(Iff),sec_avg,datestr(Tabs_all(1)));
                        [temp,hprint(Iff)]=create_percentile_distributions(Tabs_all, sumPSD,pms);
                        
                        
                        if Iff==1
                            percents.data=zeros(length(pms.fmin),length(pms.percentiles),size(temp.data,2));
                            percents.x=temp.x;
                            percents.percentiles=pms.percentiles;
%                             if suppress_output==0
%                                 uiwait(msgbox('Please click on screen to print out percentiles used'));
%                                 try
%                                     gtext(sprintf('Percentiles: %s',num2str(pms.percentiles)),'fontweight','bold','fontsize',18)
%                                 end
%                             end
                        end
                        percents.data(Iff,:,:)=temp.data;
                        percents.title{Iff}=pms.title;
                        %if suppress_output==0
                        %  save(save_tag{Iff}, 'pms','Tabs_all','sumPSD','percents');
                        % save(save_tag,'Tabs_all','PSD_all','pms','params','Batch_vars_bulkload','sec_avg','percents');
                        
                        %end
                    end  %length pms.fmin
                    
                    
                    
                    if pms.plot.image
                        
                        %%Plot contour plots if enough frequency bands
                        clear hprint
                        for Ipp=1:length(pms.percentiles)
                            hprint(Ipp)=figure;
                            
                            yy=squeeze(percents.data(:,Ipp,:));
                            Igood=~all(isnan(yy));
                            imagesc(percents.x(Igood),pms.fmin,yy(:,Igood));
                            
                            %imagesc(percents.x,pms.fmin,squeeze(percents.data(:,Ipp,:)));
                            axis('xy');
                            %datetick('x',date_tick_chc);
                            Ibin=get(gca,'xtick');
                            Istep=floor(pms.xlabel_inc/pms.x_inc);
                            Ibin=1:Istep:length(percents.x);
                            
                            %Ibin=percents.x(Ibin);
                            set(gca,'xtick',percents.x(Ibin));
                            
                            datetick('x',pms.label_style,'keepticks','keeplimits')
                            %set(gca,'xticklabel',);
                            
                            colorbar
                            caxis(pms.y_limits);
                            title(sprintf('%ith percentile, %6.2f s average, starting %s',100*pms.percentiles(Ipp),sec_avg,datestr(Tabs(1))));
                            xlabel('Date/Time','fontweight','bold','fontsize',14);
                            ylabel('Frequency (Hz)','fontweight','bold','fontsize',14);
                            grid on;
                            
                            
                            
                        end %for Ipp
                    end  %if plot.image
                    yes=menu('Redo formatting? ','Yes','No');
                    if yes==1
                        pms=get_PSD_percentile_params(pms.fmin/1000,pms.fmax/1000,pms.def);
                        if isempty(pms)
                            yes=0;
                        else
                            pms.y_label='dB re 1uPa';
                            %Nf=length(pms.fmin);
                            date_tick_chc=get_datetick_style([],pms.xlabel_inc,'datenumber');
                        end
                        
                    else
                        yes=0;
                    end
                end  %while yes
                
                %Print figure
                for Iprint=1:length(hprint)
                    tmp=datevec(pms.x_inc);
                    tmp=3600*tmp(4)+60*tmp(5)+tmp(6);
                    
                    if ~pms.plot.image
                        save_tag=sprintf('Percentile_PSD_%s_%s_fmin%iHz_fmax%iHz_interval%isec', ...
                            datestr(Tabs_all(1),30),datestr(Tabs_all(end),30),(pms.fmin(Iprint)),(pms.fmax(Iprint)),tmp);
                    else
                        save_tag=sprintf('Percentile%i_PSD_%s_%s_fmin%iHz_fmax%iHz_interval%isec', ...
                            100*pms.percentiles(Iprint),datestr(Tabs_all(1),30),datestr(Tabs_all(end),30),min(pms.fmin),max(pms.fmax),tmp);
                        
                    end
                    set(hprint(Iprint),'paperpositionmode','auto')
                    print(hprint(Iprint),'-djpeg','-r300',save_tag)
                    saveas(hprint(Iprint), save_tag, 'fig');
                    save(save_tag,'Tabs_all','PSD_all','pms','params','Batch_vars_bulkload','sec_avg','percents');
                    
                end
                
                
        end  %switch
    otherwise
        error('Batch mode not recognized');
end %switch Batch Mode


    function [hprint,save_tag]=image_PSD(twin1,twin2,file_name)
        
        if isempty(Tabs_all)
            hprint=-1;save_tag=-1;
            return
        end
        PSD_all=10*log10(PSD_all);
        
        
        hprint=figure;
        imagesc(Tabs_all,F/1000,real(PSD_all));axis('xy')
        fmin=str2num(get(handles.edit_fmin,'String'));
        fmax=str2num(get(handles.edit_fmax,'String'));
        if fmax>0
            ylimm=[fmin fmax];
            ylim(ylimm);
        end
        colorbar
        
        cmin=eval(get(handles.edit_mindB,'string'));
        cmax=cmin+eval(get(handles.edit_dBspread,'string'));
        climm=[cmin cmax];
        caxis(climm);
        xlim([twin1 twin2]);
        if isempty(date_tick_chc) %auto adjustment has failed..
            date_tick_chc=get_datetick_style('auto',twin2-twin1,'datenumber');
        end
        datetick('x',date_tick_chc,'keeplimits');
        %sec_avg=params.Nsamps*params.dn/params.Fs;
        %[~,local_fname,~]=fileparts(fname);
        titlestr=(sprintf('%s\n Start time: %s, End Time: %s, seconds averaged: %6.2f',file_name,datestr(twin1,30),datestr(twin2,30),sec_avg));
        title(titlestr);
        save_tag=[datestr(twin1,30) '_' datestr(twin2,30)];
    end
end

function date_tick_chc=get_datetick_style(input_style,val,val_units)
%Set tick date format
% input_style: string with initial style; most often "auto"
% val: number represenintg typical time interval to be plotted.  If empty
%   or less than one, return empty matrix
% val_units:  units of val number

if isempty(input_style)
    prompt = {'Datenumber format to use when plotting x-axis; examples include "mm/dd", "dd", "dd/HH", "HH" or "HH:MM":'};
    dlg_title = 'Datetick style';
    num_lines = 1;
    def = {'HH'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    input_style=answer{1};
end

if isempty(strfind(input_style,'auto'))
    date_tick_chc=input_style;
    return
end

if isempty(val)||val<0
    %uiwait(warndlg('Need to define an input number for get_datetick_style'));
    date_tick_chc=[];
    return
end

%%Convert time units to hours
switch val_units
    case 'datenumber'
        tmp=datevec(val);
        val=tmp(6)/3600+tmp(5)/60+tmp(4)+24*tmp(3);
    case 'hours'
    otherwise
        uiwait(errordlg('datetick style val units not recognized'));
        return
end

if val<1/6
    date_tick_chc='HH:MM:SS';
elseif val<=1  %less than an hour
    date_tick_chc='HH:MM'; %HH:MM
elseif val>1&&val<3
    date_tick_chc='HH';
elseif val>3
    date_tick_chc='dd';
end


end

%
function pms=get_PSD_percentile_params(fmin,fmax,def_old)
%
% pms.plot.line=false;
% pms.plot.boxplot=false;
% pms.plot.none=false;
% pms.plot.image=false;
% pms.x_label=answer{1};
% pms.x_inc=datenum(0,0,0,0,0,str2num(answer{2}));
% pms.fmin=eval(answer{3});
% pms.fmax=eval(answer{4});
% pms.y_limits=eval(answer{5});
% pms.xlabel_inc=datenum(0,0,0,0,0,str2num(answer{6}));
% pms.percentiles=eval(answer{7});
% pms.plot_chc=answer{8};
% pms.label_style
 
if exist('def_old','var')
    def=def_old;
    def{3}=num2str(1000*fmin);
    def{4}=num2str(1000*fmax);
    
else
    def = {'Date/Time','4*3600',num2str(1000*fmin),num2str(1000*fmax),'[90 140]','48*3600','[0.1 0.5 0.9]','line'};
    
end
yes=1;
while yes
    prompt = {'Time unit','Time increment (sec):','Min Frequencies (Hz) [Examples: ''[100 200 300]" or "100:10:300"]:', ...
        'Max Frequencies(Hz) [Examples: ''[100 200 300]" or "100:10:300"]', ...
        'dB limits [min max]','tick increment (sec)','percentile','Plot type (boxplot, line, image, none):'};
    dlg_title = 'Input for PSD statistics';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        pms=[];
        return
    end
    
    pms.x_label=answer{1};
    pms.x_inc=datenum(0,0,0,0,0,str2num(answer{2}));
    pms.fmin=eval(answer{3});
    pms.fmax=eval(answer{4});
    pms.y_limits=eval(answer{5});
    pms.xlabel_inc=datenum(0,0,0,0,0,str2num(answer{6}));
    pms.percentiles=eval(answer{7});
    pms.plot_chc=answer{8};
    
    if rem(length(pms.percentiles),2)==0
        uiwait(msgbox('Must be an odd number of percentiles'));
    else
        yes=0;
    end
    date_tick_chc=get_datetick_style('auto',pms.xlabel_inc,'datenumber');
    
    pms.label_style=date_tick_chc;
    pms.def=answer;
    
    prompt = {'Datenumber Format, examples include "mm/dd", "dd", "dd/HH", "HH" or "HH:MM":'};
    dlg_title = 'Check the datenumber format...';
    num_lines = 1;
    def2{1}=date_tick_chc;
    answer = inputdlg(prompt,dlg_title,num_lines,def2);
    pms.label_style=answer{1};
    
    pms.plot.line=false;
    pms.plot.boxplot=false;
    pms.plot.none=false;
    pms.plot.image=false;
    
    switch pms.plot_chc
        case 'line'
            pms.plot.line=true;
        case 'image'
            pms.plot.image=true;
        case 'boxplot'
            pms.plot.boxplot=true;
        case 'none'
            pms.plot.none=true;
        otherwise
            uiwait(msgbox('Invalid plot format selected, will just plot nothing'));
            pms.plot.none=true;
    end
              
    %%Determine what type of plot to use, if any
    if length(pms.fmin)>5&&~pms.plot.none
        ButtonName = questdlg(sprintf('You have chosen over %i frequency bands.  Use image plot instead?',length(pms.fmin)), ...
            'Plot percentile image?', ...
            'Yes', 'No', 'No');
        if strcmp(ButtonName,'Yes')
            pms.plot.image=true;
            pms.plot.line=false;
            pms.plot.boxplot=false;
        end
        
    end
                
               
    
    % if pms.xlabel_inc>=datenum(0,0,0,12,0,0) %label spacing on order of days
    %     pms.label_style='dd';
    % elseif pms.xlabel_inc>=datenum(0,0,0,1,0,0);  % If label spacing is over one hour
    %     pms.label_style='HH';
    % else
    %     pms.label_style='HH:MM';
    % end
end
end

% --------------------------------------------------------------------
function MenuItem_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to MenuItem_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Batch_type	=	get(hObject, 'Label');
Batch_mode	=	input_batchmode(Batch_type);

if	isempty(Batch_mode)
    uiwait(errordlg('Processing cancelled!'));
    return;
end

if isempty(strfind(Batch_mode,'Start Bulk'))
    Batch_vars.mean	=	'yes';	Batch_desc{1}	=	'Enter yes for averaged linear spectrum (instead of just percentiles)';
    Batch_vars.percentile	=	'[10 50 90]';	Batch_desc{2}	=	'Percentile display';
    
    Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, Batch_type);
    percentile=str2num(Batch_vars.percentile)/100;
    
    if	isempty(Batch_vars)
        uiwait(errordlg('Processing cancelled!'));
        return;
    end
end
% %



%Graphic parameters
line_width=2;


switch	Batch_mode
    case 'Process Visible Window'
        uiwait(msgbox(['Processing all data within window for ' Batch_type]));
        
        %%Make a new figure and plot spectrum.
        Ileg=1;
        SS_mean=[];
        hprint=figure;
        
        if strcmp(Batch_vars.mean,'yes')
            SS_mean=mean(handles.sgram.B');
            plot(handles.sgram.F,10*log10(SS_mean),'k','linewidth',line_width);grid on
            hold on
            leg_str{Ileg}='mean';Ileg=Ileg+1;
            
        end
        
        if ~isempty(percentile)
            xsort=sort(handles.sgram.B,2)';
            Iout=round(percentile*size(xsort,1));
            SS_percentile=xsort(Iout,:);
            plot(handles.sgram.F,10*log10(SS_percentile),'linewidth',line_width);
            
            for II=1:length(percentile)
                leg_str{Ileg}=num2str(100*percentile(II));
                Ileg=Ileg+1;
            end
            
        end
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Frequency(Hz)');ylabel(' PSD (dB re 1 uPa^2/Hz)');
        
        legend(leg_str);
        titlestr=(sprintf('Start time: %s, Time processed: %6.2f sec',datestr(handles.tdate_start),handles.tlen));
        title(titlestr);
        ylimm=ylim;
        ylimm(1)=0;
        ylim(ylimm);
        
        yess=menu('Save Spectrum Data and figure?','Yes','No');
        if yess==1
            cd(handles.outputdir);
            F=handles.sgram.F;
            save_str=sprintf('Spectrum_%s_%s', datestr(handles.tdate_start,30),datestr(handles.tdate_start+datenum(0,0,0,0,0,handles.tlen),30));
            save(save_str,'F','percentile','leg_str','SS_percentile','SS_mean','titlestr');
            uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
            
            orient landscape
            print(hprint,'-djpeg','-r300',save_str);
        end
        
        return
    case {'Start Bulk Processing File','Start Bulk Processing Folder'}
        
        if strcmp(Batch_mode,'Start Bulk Processing File')
            uiwait(msgbox(sprintf('Processing all data in file %s\n for %s',fullfile(handles.mydir,handles.myfile), Batch_type)));
            param.exten=['*' handles.myfile];
        else
            uiwait(msgbox(sprintf('Processing all data in folder %s\n for %s',fullfile(handles.mydir), Batch_type)));
            
            param.exten=['*' handles.myext];
        end
        data_folder_name		=	uigetdir('.', 'Select a folder to store analysis results:');
        param.dir_out=data_folder_name;
        
        param.file_dir=handles.mydir;
        contents=get(handles.popupmenu_Nfft,'String');
        param.Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});
        contents=get(handles.popupmenu_ovlap,'String');
        param.ovlap=round(param.Nfft*str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100);
        param.Fs=handles.sgram.Fs;
        param.channel=str2num(get(handles.edit_chan,'String'));
        param.dumpsize=100000;
        param.nstart=0;
        param.nsamples=0;
        param.f_low=1000*str2num(get(handles.edit_fmin,'String'));
        param.f_high=1000*str2num(get(handles.edit_fmax,'String'));
        param.sec_avg=param.Nfft/param.Fs;
        
        cd(param.dir_out);
        write_Java_script('PSD',param);
        ! ./masterPSD.scr > outt.txt &
        return
    case 'Load Bulk Processing'
        
        
        Batch_vars_bulkload.start_time	=	'select';
        Batch_desc{1}	=	'Time to begin loading (e.g. "here" to use visible start time, "datenum(2011,1,2,0,0,0)", "file start" for current file , "select" to choose file from list, "folder start" for first file in folder)';
        Batch_vars_bulkload.end_time	=	'select';
        Batch_desc{2}	=	'Time to end loading (e.g. "here" to use GUI end time, "2" for two hours from start time, "datenum(2011,1,2,0,0,0)", "file end" for time at current file end, "folder end" to import through end of folder)';
        
        
        Batch_vars_bulkload.dB_bins = '30:2:140';        Batch_desc{3} = 'vector of dB levels to build long-term histograms';
        Batch_vars_bulkload.duty_cycle='0';        Batch_desc{4} = 'Type "1" to process a duty cycle (e.g. diel or weekend analysis)';
        Batch_vars_bulkload.cycle_duration='0';     Batch_desc{5}='Duty cycle:  Hours in a complete cycle.  Cycle begins at date/time in first row above.  Example: "24*7" defines a week';
        Batch_vars_bulkload.process_duration='0';   Batch_desc{6}='Duty cycle:  Hours to process, beginning at the start of a cycle.  Example: "24*2" defines a weekend';
        
        Batch_vars_bulkload	=	input_batchparams(Batch_vars_bulkload, Batch_desc, Batch_type);
        duty_cycle_chc=logical(eval(Batch_vars_bulkload.duty_cycle));
        
        if	isempty(Batch_vars_bulkload)
            uiwait(errordlg('Processing cancelled!'));
            return;
        end
        
        
        % Load bulk run information
        %[tabs_folder_start,tabs_folder_end,tabs_start,tabs_end,Other_FileNames,Icurrent_file,Ifinal_file,FF]=load_PSD_Bulk_Run(handles, Batch_vars_bulkload);
        bulk_params=load_PSD_Bulk_Run(handles, Batch_vars_bulkload);
        
        tabs_folder_start=bulk_params.tabs_folder_start;
        tabs_folder_end=bulk_params.tabs_folder_end;
        tabs_start=bulk_params.tabs_start;
        tabs_end=bulk_params.tabs_end;
        Other_FileNames=bulk_params.Other_FileNames;
        Icurrent_file=bulk_params.Icurrent_file;
        Ifinal_file=bulk_params.Ifinal_file;
        FF=bulk_params.FF;
        
        %%Translate duty cycle...
        if duty_cycle_chc
            cycle_duration=datenum(0,0,0,eval(Batch_vars_bulkload.cycle_duration),0,0);
            process_duration=datenum(0,0,0,eval(Batch_vars_bulkload.process_duration),0,0);
            
            current_cycle_start=tabs_start;
            %next_cycle_start=tabs_start+cycle_duration;
            process_end=tabs_start+process_duration;
            
        end
        
        
        %Start processing
        dB_bins=str2num(Batch_vars_bulkload.dB_bins);
        hist_total=zeros(length(FF),length(dB_bins));
        tabs_loop_begin=tabs_start;
        for I=Icurrent_file:Ifinal_file
            fname=Other_FileNames(I).name;
            fprintf('Processing %s...\n',fname);
            [PSD,F,Tsec,Tabs,params]=read_Java_PSD(fname,tabs_loop_begin,Inf);  %Read an entire file into memory
            Istrip=find(Tabs<=tabs_end);
            PSD=PSD(:,Istrip);Tabs=Tabs(Istrip);
            
            
            %%Extract processing times.
            
            if ~duty_cycle_chc
                Igood=1:length(Tabs);
            else
                Igood=[];
                while current_cycle_start<max(Tabs)
                    Igood=[Igood find(Tabs>=current_cycle_start&Tabs<=process_end)];  %This has to be true for first file
                    current_cycle_start=current_cycle_start+cycle_duration;
                    process_end=process_end+cycle_duration;
                end
                
                %%Back up a step in case desired process time runs into
                %%next file
                current_cycle_start=current_cycle_start-cycle_duration;
                process_end=process_end-cycle_duration;
            end
            
            %If we move to next file, start at the beginning
            tabs_loop_begin=0;
            
            if isempty(Igood)
                continue
            end
            
            PSD=10*log10(PSD);
            
            for If=1:length(FF)
                Ncount=histc(PSD(If,Igood),dB_bins);
                hist_total(If,:)=hist_total(If,:)+Ncount;
            end
            
            
            
        end  %Icurrent_file
        
        %%%%%Convert bulk histogram into percentile distribution.
        cum_dist=cumsum(hist_total');
        cum_dist=cum_dist./(ones(length(dB_bins),1)*max(cum_dist));
        
        CC=contourc(FF,dB_bins,cum_dist,percentile);
        
        SS_percentile=zeros(length(percentile),length(FF));
        for Ip=1:length(percentile)
            Npt=CC(2,1);
            p_check{Ip}=num2str(CC(1,1));
            SS_raw{Ip}=CC(2,2:(Npt+1));
            Freq{Ip}=CC(1,2:(Npt+1));
            CC=CC(:,(Npt+2):end);
            SS_percentile(Ip,:)=interp1(Freq{Ip},SS_raw{Ip},FF);
            
        end
        
        %%Debug check
        hprint=figure;
        subplot(2,1,1)
        imagesc(FF,dB_bins,cum_dist);colorbar;
        hold on
        for Ip=1:length(percentile)
            plot(Freq{Ip},SS_raw{Ip},'y');
            
        end
        plot(FF,SS_percentile,'o');
        titlestr=(sprintf('Start time: %s, End Time: %s ',datestr(tabs_start,30),datestr(tabs_end,30)));
        title(titlestr);
        
        grid on
        
        subplot(2,1,2)
        plot(FF,SS_percentile,'linewidth',line_width);
        grid on
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Frequency(Hz)');ylabel(' PSD (dB re 1 uPa^2/Hz)');
        legend(p_check);
        
        yess=menu('Save Spectrum Data and figure?','Yes','No');
        if yess==1
            cd(handles.outputdir);
            F=FF;
            leg_str=p_check;
            save_str=sprintf('Spectrum_%s_%s', datestr(tabs_start,30),datestr(tabs_end,30));
            save(save_str,'F','percentile','leg_str','SS_percentile','titlestr','Batch_vars','Batch_vars_bulkload');
            uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
            
            orient landscape
            print(hprint,'-djpeg','-r300',save_str);
        end
    otherwise
        error('Batch mode not recognized');
end



end

function bulk_params=load_PSD_Bulk_Run(handles,Batch_vars_bulkload)
%Note that a file does not have to be loaded for this to work...
% Output parameters:
%       tabs_folder_start:  datenumber of time of start of first psd file
%           in folder.
%       tabs_folder_end:  datenumber of time of end of last psd file
%           in folder.
%       tabs_start:  datenumber of time at start of file loaded
%       tabs_end:  datenumber of time at end of file loaded
%       Other_FileNames: cell array of PSD files that share common
%           processing parameters in a folder.
%       Icurrent_file:  Index of Other_FileNames that lists the PSD file
%           associated with the loaded file.  Is set to '1' if all files
%           desired.
%       Ifinal_file:  number of PSD files that occur after the current loaded
%           file.
%       FF: vector of frequencies associated with power spectral density.

%

bulk_params.tabs_folder_start=[];
bulk_params.tabs_folder_end=[];
bulk_params.tabs_start=[];
bulk_params.tabs_end=[];
bulk_params.Other_FileNames=[];
bulk_params.Icurrent_file=[];
bulk_params.Ifinal_file=[];
bulk_params.FF=[];

if isfield(handles,'myfile')
    [~,token,extt] = fileparts(handles.myfile);
else
    token=[];
end

%Select first example file
if ~strcmp(Batch_vars_bulkload.start_time,'select')
    
    dialog_title	=	'Select an example Bulk file to load: Note that I am looking in analysis directory, not data directory';
else
    dialog_title	=	'Select the first Bulk file to load: Note that I am looking in analysis directory, not data directory';
    token=[];
end
[FileName,DirectoryName] = uigetfile([token '*.psd'],dialog_title);
if isnumeric(FileName)
    disp('No file selected');
    return;
end

%Select second example file if needed
if strcmp(Batch_vars_bulkload.end_time,'select')
    dialog_title	=	'Select the final Bulk file to load: Note that I am looking in analysis directory, not data directory';
    
    [FileNameEnd,DirectoryNameEnd] = uigetfile([token '*.psd'],dialog_title);
    if isnumeric(FileName)
        disp('No second file selected');
        return;
    end
else
    FileNameEnd='';
end

fprintf('Changing working directory to %s\n',DirectoryName);
cd(DirectoryName);
Iscore=min(strfind(FileName,'_chan'));

token2=FileName((Iscore+1):end);  %PSD files having the same input parameters..

Other_FileNames=dir(['*_' token2]);

%Sort by numerical order, not ASCII order.


if isempty(Other_FileNames)
    uiwait(errordlg('load_PSD_Bulk_Run cannot find matches to file you selected:','Files not found!'));
    return;
end

tabs_select_start=[];
tabs_select_end=[];
tabs_selectEnd_start=[];
tabs_selectEnd_end=[];
for II=1:length(Other_FileNames)
    [~,~,~,~,params]=read_Java_PSD(Other_FileNames(II).name,0,Inf,1); %Read header only
    tabs_file_start(II)=params.tstart_file;
    tabs_file_end(II)=params.tend_file;
    
    
end

%Sort by start time
[tabs_file_start,Isort]=sort(tabs_file_start);
tabs_file_end=tabs_file_end(Isort);
Other_FileNames=Other_FileNames(Isort);

for II=1:length(Other_FileNames)
    if strcmp(Other_FileNames(II).name,FileName)
        tabs_select_start=tabs_file_start(II);
        tabs_select_end=tabs_file_end(II);
        Iselect_file=II;
    end
    
    if strcmp(Other_FileNames(II).name,FileNameEnd)
        tabs_selectEnd_start=tabs_file_start(II);
        tabs_selectEnd_end=tabs_file_end(II);
        IselectEnd_file=II;
        
    end
    
end


tabs_folder_start=tabs_file_start(1);  %datenumber of start of data in folder (assuming all filenames have chronological order)
tabs_folder_end=tabs_file_end(end);  %datenumber of end of data in folder (assuming all filenames have chronological order)


[~,FF,~,~,~]=read_Java_PSD(Other_FileNames(1).name,0,10);  %Get FF vector...

% [~,FF,~,~,params]=read_Java_PSD(Other_FileNames(Icurrent_file).name,0,Inf,1);
% tabs_end=params.tend_file;  %datenumber of end of data in focal file (assuming all filenames have chronological order)
% tabs_start=params.tstart_file;  %datenumber of start of data in focal file (assuming all filenames have chronological order)


%%Translate start time.  Need to define tabs_start and Icurrent_file
switch lower(Batch_vars_bulkload.start_time)
    case 'select'
        if ~isempty(tabs_select_start)
            tabs_start=tabs_select_start;
        else
            uiwait(msgbox('Start file not selected'));
            return
        end
        %Reset other variables
        Icurrent_file=Iselect_file;
    case 'here'  %Use edit window
        tabs_start=datenum(get(handles.edit_datestr,'String'));
        Icurrent_file=find(tabs_start>=tabs_file_start, 1, 'last' );
    case 'file start'  %start of file
        tabs_start=datenum(get(handles.edit_datestr,'String'));
        Icurrent_file=find(tabs_start>=tabs_file_start, 1, 'last' );
        tabs_start=tabs_file_start(Icurrent_file);
        
    case 'folder start'  %start of folder
        tabs_start=tabs_folder_start;
        %Reset other variables
        Icurrent_file=1;
        %Ifinal_file=length(Other_FileNames);
    otherwise %check for datenum
        if strfind(Batch_vars_bulkload.start_time,'datenum')
            try
                tabs_start_want=eval(Batch_vars_bulkload.start_time);
                if tabs_start_want<tabs_folder_start || tabs_start_want>tabs_folder_end
                    uiwait(errordlg('Can''t request a time outside selected folders''s start time; start time set to folder start','Bulk Processing Menu Error'));
                    tabs_start_want=tabs_folder_start;
                    tabs_start=tabs_start_want;
                    Icurrent_file=1;
                    
                else  %location appropriate start time...
                    %Identify current file...
                    Icurrent_file=find(tabs_start_want>=tabs_file_end, 1, 'last' )+1; %Do this because of duty-cycled files
                    tabs_start=max([tabs_start_want tabs_file_start(Icurrent_file)]);
                end
            catch
                uiwait(errordlg('Can''t process start_time','Bulk Processing Menu Error'));
                return;
            end
            
            
        end
end


%%Translate end time;  Need to define tabs_end and Ifinal_file
switch lower(Batch_vars_bulkload.end_time)
    case 'select'
        if ~isempty(tabs_selectEnd_end)
            tabs_end=tabs_selectEnd_end;
            
        else
            uiwait(msgbox('End file not selected'));
            return
        end
        Ifinal_file=IselectEnd_file;
        
    case 'here' %Use edit window
        tabs_end=datenum(get(handles.edit_datestr,'String'))+datenum(0,0,0,0,0,str2num(get(handles.edit_winlen,'string')));
        Ifinal_file=Icurrent_file;  %Load one file only
    case 'folder end'
        tabs_end=tabs_folder_end;
        Ifinal_file=length(Other_FileNames);
    case 'file end'  %end of file
        %tabs_end=Inf;
        Ifinal_file=Icurrent_file;  %Load one file only.
        tabs_end=tabs_file_end(Ifinal_file);
    otherwise %check for datenum
        if strfind(Batch_vars_bulkload.end_time,'datenum')
            try
                tabs_end_want=eval(Batch_vars_bulkload.end_time);
                if tabs_end_want<tabs_folder_start || tabs_end_want>tabs_folder_end
                    uiwait(errordlg('Can''t request a time outside selected folder''s end time; end time set to folder end','Bulk Processing Menu Error'));
                    tabs_end_want=tabs_folder_end;
                    tabs_end=tabs_end_want;
                    Ifinal_file=length(Other_FileNames);
                else  %location appropriate start time...
                    %Identify current file...
                    
                    Ifinal_file=find(tabs_end_want<=tabs_file_start, 1, 'first' )-1;  %Do this because of duty-cycled files
                    tabs_end=min([tabs_end_want tabs_file_end(Ifinal_file)]);
                end
                
            catch
                uiwait(errordlg('Can''t process end_time','Bulk Processing Menu Error'));
                return;
                
            end
        elseif isnumeric(str2num(Batch_vars_bulkload.end_time))  %%Assume end time is duration in hours after start time
            tabs_end=tabs_start+datenum(0,0,0,str2num(Batch_vars_bulkload.end_time),0,0);
            Ifinal_file=find(tabs_end<=tabs_file_end, 1, 'first' );
            if isempty(Ifinal_file)
                uiwait(errordlg('Can''t process end_time','Bulk Processing Menu Error'));
                return;
            end
            
        end
end

fprintf('Folder start time: %s, File start time: %s, File end time: %s, Folder end time: %s\n', ...
    datestr(tabs_folder_start,30),datestr(tabs_start,30),datestr(tabs_end,30),datestr(tabs_folder_end,30));


bulk_params.tabs_folder_start=tabs_folder_start;
bulk_params.tabs_folder_end=tabs_folder_end;
bulk_params.tabs_start=tabs_start;
bulk_params.tabs_end=tabs_end;
bulk_params.Other_FileNames=Other_FileNames;
bulk_params.Icurrent_file=Icurrent_file;
bulk_params.Ifinal_file=Ifinal_file;
bulk_params.FF=FF;

end

% --------------------------------------------------------------------
function MenuItem_eventdet_Callback(hObject, eventdata, handles)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Batch_type	=	get(hObject, 'Label');
Batch_mode	=	input_batchmode(Batch_type);

if	isempty(Batch_mode)
    uiwait(errordlg('Processing cancelled!'));
    return;
end

switch	Batch_mode
    case 'Process Visible Window'
        %%%%%%%%Energy_detector_review.m%%%%%%%%%%%%%%%
        %function Energy_detector_review(param,fname,tstart,twant,tlen,tview)
        %
        %  Aaron Thode
        %  December 9. 2008
        %  Run a short segment of time through Energy detector for insight and
        %  adjustment.
        %
        %  Input:
        %       param: structure of parameters for energy detector; must contain
        %           field param.energy.
        %       fname: complete pathname to raw data file
        %       tstart: datenumber of start time of file
        %       twant:  datenumber of desired start time energy detector processing
        %       tlen:   number of seconds to process
        %       tview: two element vector giving bounds for viewing output (sec)
        %            This permits viewing of a short transient window while
        %            permitting a long "start up" time for the equalization.
        
        
        burn_in_time=1; %minutes
        
        [param,param_desc]=load_energy_parameters;
        if isempty(param)
            return
        end
        param.energy=param;
        
        
        %tel=datevec(datenum(param.nstart)-datenum(handles.tdate_min));
        %sec=tel(:,6)+60*tel(:,5)+3600*tel(:,4)+24*3600*tel(:,3);
        tlen=(param.nsamples);
        tstart=datenum(param.energy.nstart);
        %param.nstart=round(param.Fs*sec);
        
        climm=str2num(get(handles.edit_mindB,'String'))+[0 str2num(get(handles.edit_dBspread,'String'))];
        
        data_all=Energy_detector_review(param,fullfile(handles.mydir,handles.myfile),  ...
            handles.tdate_min,tstart,tlen,climm,handles.axes1,[60*burn_in_time tlen]);
        
        %%Return 'param' to a form that can be used by the input_batch
        %%  function
        fields=fieldnames(param.energy);
        param.energy.f_low=param.energy.f_low/1000;
        param.energy.f_high=param.energy.f_high/1000;
        param.energy.bandwidth=param.energy.bandwidth/1000;
        param.energy.ovlap=round(100*param.energy.ovlap/param.Nfft);
        
        for II=1:length(fields)
            tmp=num2str(param.energy.(fields{II}));
            if ~isempty(tmp)
                param.energy.(fields{II})=tmp;
            end
        end
        
        handles.energy_parameters_default=param.energy;
        handles.energy_parameters_default_desc=param_desc;
        !rm *sel *scr *snips *detsum out*txt
        guidata(hObject, handles);
        
    case {'Start Bulk Processing File','Start Bulk Processing Folder'}
        
        burn_in_time=0;
        param=load_energy_parameters;
        
        if isempty(param)
            return
        end
        
        script_name=writeEnergyDetectorCshell(param);
        pause(1);
        eval(sprintf('!./%s > outt_click.txt &',script_name));
        pause(1);
        
        handles.energy_parameters_default=param;
        
        guidata(hObject, handles);
        return
    case 'Load Bulk Processing'
        
        %info_str='Select the appropriate detsum file using "Select Directory" button in the Annotations Section';
        %disp(info_str);
        %h	=	msgbox(info_str);
        
        dialog_title	=	'Select location for Event Detector import file(s):';
        start_path		=	handles.mydir;
        folder_name		=	uigetdir(start_path, dialog_title);
        
        if isnumeric(folder_name)
            return;
        end
        mydir=handles.mydir;
        
        handles.mydir=folder_name;
        cd(handles.mydir);  %Change location to be inside this folder, since we assume it is a data analysis folder.
        
        %%Make sure all files with detsum and snips in this folder can be
        %%imported.
        
        [filename, pathname] = uigetfile('*.detsum', 'Pick a representative detection file (I''ll look at all files that share same parameter set:');
        if isequal(filename,0) || isequal(pathname,0)
            disp('User pressed cancel')
            return
        end
        
        
        %%Grab template for detsum file; look for a 'xxtoyyHz' substring;
        Istart=strfind(filename,'_');
        Idot=max(strfind(filename,'.'))-1;
        if isfield(handles,'myfile')
            myfile_org=handles.myfile;
        else
            myfile_org=[];
        end
        handles.myfile=['*' filename(Istart(7):Idot)];
        handles		=	load_notes_file(handles, folder_name);  %Bulk Load Event detector
        handles.myfile=myfile_org;
        if ~isfield(handles,'tdate_start')
            handles.tdate_start=[];
        end
        
        %guidata(hObject, handles);
        
        file_path	=	handles.notes.file_path;
        if	isempty(file_path)
            warning('Notes file/folder not set, file not saved');
            return;
        end
        
        %	Save whole data structure, allows implicit expansion of named variables
        %	Only save if changed
        GUI_params		=	save_gui_params(handles);
        Data=handles.notes.Data;
        save(file_path, 'Data', 'GUI_params');
        handles.notes.saved	=	true;
        set(handles.pushbutton_notes_stats,'Enable','on');
        handles.mydir=mydir;
        cd(handles.mydir);
        
        guidata(hObject, handles);
        
        
    otherwise
        error('Batch mode not recognized');
end

    function [param,param_desc]=load_energy_parameters
        
        if isfield(handles,'energy_parameters_default')&&isfield(handles,'energy_parameters_default_desc')&&~isempty(strfind(Batch_mode,'Bulk'))
            %%Logic assumes that a normal user has used "Process Visible
            %%Window and is now starting a bulk run
            param=handles.energy_parameters_default;
            param_desc=handles.energy_parameters_default_desc;
            %param_desc=load_energy_parameter_description;
            param.nstart='0';
            param.nsamples='0';
            if strcmp(Batch_mode,'Start Bulk Processing File')
                uiwait(msgbox(sprintf('Processing all data in file %s\n for %s',fullfile(handles.mydir,handles.myfile), Batch_type)) );
                param.exten=['*' handles.myfile];
                
            else
                uiwait(msgbox(sprintf('Processing all data in folder %s\n for %s',fullfile(handles.mydir), Batch_type)) );
                param.exten=['*' handles.myext];
            end
            
        else
            K=1;
            if strcmp(Batch_mode,'Start Bulk Processing File')
                uiwait(msgbox(sprintf('Processing all data in file %s\n for %s',fullfile(handles.mydir,handles.myfile), Batch_type)) );
                param.exten=['*' handles.myfile];
            else
                uiwait(msgbox(sprintf('Processing all data in folder %s\n for %s',fullfile(handles.mydir), Batch_type)) );
                param.exten=['*' handles.myext];
            end
            param_desc{K}='Filename or file extension to process';K=K+1;
            data_folder_name		=	uigetdir('.', 'Select a folder to store analysis results:');
            
            param.dir_out=data_folder_name;         param_desc{K}='Output directory';K=K+1;
            param.file_dir=handles.mydir;   param_desc{K}='Directory containing files';K=K+1;
            contents=get(handles.popupmenu_Nfft,'String');
            param.Nfft=(contents{get(handles.popupmenu_Nfft,'Value')});param_desc{K}='FFT size';K=K+1;
            
            contents=get(handles.popupmenu_ovlap,'String');
            param.ovlap=(contents{get(handles.popupmenu_ovlap,'Value')});
            param_desc{K}='percent overlap';K=K+1;
            
            try
                param.Fs=num2str(handles.sgram.Fs);param_desc{K}='Sampling rate (will be overidden for most files)';K=K+1;
            catch
                uiwait(msgbox(sprintf('Need to load spectrogram first')));
                return
            end
            param.channel=(get(handles.edit_chan,'String'));param_desc{K}='Data channel (starting from 1):';K=K+1;
            param.dumpsize='100000';param_desc{K}='JAVA dump size (points):';K=K+1;
            
            if ~isempty(strfind(Batch_mode,'Start Bulk Processing'))
                param.nstart='0';param_desc{K}='number of samples to skip before processing:';K=K+1;
                param.nsamples='0';param_desc{K}='number of samples to process:';K=K+1;
            else
                test_time=datenum(get(handles.edit_datestr,'String'))-datenum(0,0,0,0,burn_in_time,0);
                if test_time<handles.tdate_min
                    uiwait(msgbox(sprintf('You need to start at least one minute into file to allow processor to burn in')));
                    return
                end
                param.nstart=datestr(test_time);param_desc{K}='Time to begin processing:';K=K+1;
                param.nsamples=num2str(str2num(get(handles.edit_winlen,'String'))+60*burn_in_time);param_desc{K}='seconds to process:';K=K+1;
                
            end
            param.f_low=(get(handles.edit_fmin,'String'));param_desc{K}='minimum frequency (kHz)';K=K+1;
            param.f_high=(get(handles.edit_fmax,'String'));param_desc{K}='maximum frequency (kHz)';K=K+1;
            %param.bandwidth=param.f_high-param.f_low;
            
            
            param.eq_time='10';   param_desc{K}='Equalization time (s): should be roughly twice the duration of signal of interest';K=K+1;
            param.bandwidth='.05';     param_desc{K}='Bandwidth of detector in kHz';K=K+1;
            param.threshold='10';  param_desc{K}='threshold in dB to accept a detection';K=K+1;
            param.snips_chc='1';  param_desc{K}='0 for no snips file, 1 for snips file of one channel, 2 for snips file of all channels';K=K+1;
            param.bufferTime='0.5'; param_desc{K}='buffer Time in seconds to store before and after each detection snip, -1 suppress snips file';K=K+1;
            param.TolTime='1e-4';  param_desc{K}='Minimum time in seconds that must elapse for two detections to be listed as separate';K=K+1;
            param.MinTime='0';     param_desc{K}='Minimum time in seconds a required for a detection to be logged';K=K+1;
            param.MaxTime='3';     param_desc{K}= 'Maximum time in seconds a detection is permitted to have';K=K+1;
            param.debug='0';       param_desc{K}= '0: do not write out debug information. 1:  SEL output.  2:  equalized background noise. 3: SNR.';K=K+1;
            
            
            
        end
        
        
        param= 	input_batchparams(param, param_desc, Batch_type);
        if isempty(param)
            uiwait(msgbox(sprintf('Energy Detector aborted','replace')));
            return
        end
        fields=fieldnames(param);
        for II=1:length(fields)
            tmp=str2num(param.(fields{II}));
            if ~isempty(tmp)
                param.(fields{II})=tmp;
            end
        end
        param.f_low=param.f_low*1000;
        param.f_high=param.f_high*1000;
        param.bandwidth=1000*param.bandwidth;
        %param.bandwidth=round(2*(param.f_high-param.f_low)/(1+param.Nbands));  %50% overlap between bands
        param.ovlap=round(param.ovlap*param.Nfft/100);
    end

    function param_desc=load_energy_parameter_description
        K=1;
        param_desc{K}='Filename or file extension to process';K=K+1;
        
        param_desc{K}='Output directory';K=K+1;
        param_desc{K}='Directory containing files';K=K+1;
        param_desc{K}='FFT size';K=K+1;
        
        param_desc{K}='percent overlap';K=K+1;
        
        param_desc{K}='Sampling rate (will be overidden for most files)';K=K+1;
        
        param_desc{K}='Data channel (starting from 1):';K=K+1;
        param_desc{K}='JAVA dump size (points):';K=K+1;
        
        param_desc{K}='Time to begin processing:';K=K+1;
        param_desc{K}='seconds to process:';K=K+1;
        
        param_desc{K}='minimum frequency (kHz)';K=K+1;
        param_desc{K}='maximum frequency (kHz)';K=K+1;
        
        param_desc{K}='Equalization time (s): should be roughly twice the duration of signal of interest';K=K+1;
        param_desc{K}='Bandwidth of detector in kHz';K=K+1;
        param_desc{K}='threshold in dB to accept a detection';K=K+1;
        param_desc{K}='0 for no snips file, 1 for snips file of one channel, 2 for snips file of all channels';K=K+1;
        param_desc{K}='buffer Time in seconds to store before and after each detection snip, -1 suppress snips file';K=K+1;
        param_desc{K}='Minimum time in seconds that must elapse for two detections to be listed as separate';K=K+1;
        param_desc{K}='Minimum time in seconds a required for a detection to be logged';K=K+1;
        param_desc{K}= 'Maximum time in seconds a detection is permitted to have';K=K+1;
        param_desc{K}= '0: do not write out debug information. 1:  SEL output.  2:  equalized background noise. 3: SNR.';K=K+1;
        
    end  %load_energy_parameter_description
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask user for mode of processing operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return string description processing choice...
function	[Batch_mode,Mode_options]	=	input_batchmode(Batch_type)
%Mode_options	=	{'Process Visible Window', 'Start Bulk Processing File', 'Start Bulk Processing Folder','Load Bulk Processing'};
Mode_options	=	{'Process Visible Window', 'Start Bulk Processing', 'Load Bulk Processing'};
qstring	=	'Which operation do you want?';
title	=	Batch_type;

button	=	questdlg(qstring, title ,...
    Mode_options{1}, Mode_options{2}, Mode_options{3},...
    Mode_options{1});

Batch_mode	=	button;

switch Batch_mode
    case 'Start Bulk Processing'
        Mode_options	=	{'File','Folder'};
        title	=	Batch_type;
        
        button	=	questdlg(qstring, title ,...
            Mode_options{1}, Mode_options{2}, ...
            Mode_options{1});
        
        Batch_mode	=	[Batch_mode ' ' button];
        
end


end

% Asks user for specified batch processing parameters-return strings
function	[Batch_vars, Batch_names]	=	input_batchparams(Batch_vars, Batch_desc, Batch_type, num_out)

if	isstruct(Batch_vars)
    names		=	fieldnames(Batch_vars);
elseif	iscell(Batch_vars)
    names		=	Batch_vars;
    Batch_vars	=	struct;
else
    error('Batch_vars must be either a structure, or cell array of variable names');
end

N_fields	=	length(names);
if	length(Batch_desc) ~= N_fields
    error('# of Descriptors differs from # of fields');
end

%	Default values
defaults	=	cell(N_fields,1);
if	isstruct(Batch_vars)
    for	ii	=	1:N_fields
        value	=	Batch_vars.(names{ii});
        defaults{ii}	=	num2str(value);
    end
end

%	Size of each prompt field
num_lines		=	1;

%	Title for window
dlgTitle	=	['Parameters for  ' Batch_type];


%	Create input dialogbox
answer		=	inputdlg(Batch_desc, dlgTitle, num_lines, defaults,'on');


%	Put variables into output structure
if	isempty(answer)
    Batch_vars	=	[];
    return;
end

for	ii	=	1:N_fields
    Batch_vars.(names{ii})	=	answer{ii};
    if nargin>3
        try
            Batch_vars.(names{ii})	=	eval(answer{ii});
        end
        
    end
end

Batch_names=fieldnames(Batch_vars);

end

% Asks user for specified batch processing parameters-return numerics

%%	Button callbacks

% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	Stop audio playback if it's running
if	~isempty(handles.audioplayer) && isplaying(handles.audioplayer)
    stop(handles.audioplayer);
end

%	Disable update buttons while loading/processing
set(handles.pushbutton_update, 'Enable', 'off');
set(handles.pushbutton_next, 'Enable', 'off');
set(handles.pushbutton_prev, 'Enable', 'off');

%	Call actual function to produce figure
handles		=	load_and_display_spectrogram(handles);

%	Call helper function to overlay event notes
handles				=	plot_events(handles);

%	Update static date text
set(handles.text_datestr_demo,'String',datestr(handles.tdate_start));

%	Audio object is no longer valid
handles.audio_stale		=	true;
handles.fig_updated		=	true;

%	Turn on buttons that required Update
toggle_initial_buttons(handles, 'on');

%	Renable buttons
set(handles.pushbutton_update, 'Enable', 'on');
set(handles.pushbutton_next, 'Enable', 'on');
set(handles.pushbutton_prev, 'Enable', 'on');


guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_datestr_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
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

function edit_datestr_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to edit_datestr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%tmin=datenum(get(handles.text_mintime,'String'));
%tmax=datenum(get(handles.text_maxtime,'string'));
tmin	=	handles.tdate_min;
tmax	=	handles.tdate_max;
if isempty(get(hObject,'String'))
    return
end
try
    new_date	=	datenum(get(hObject,'String'));
catch %#ok<*CTCH>
    uiwait(errordlg('Incorrect datestr'));
    new_date =	[];
end

if	~isempty(new_date)
    if		new_date < tmin
        new_date	=	tmin;
    elseif	new_date > tmax
        new_date	=	tmax;
    end
    newval	=	(new_date-tmin)/(tmax-tmin);
    handles.tdate_start		=	new_date;
    set(hObject,'String', datestr(new_date, 'dd-mmm-yyyy HH:MM:SS.FFF'));
    set(handles.slider_datestr,'Value',newval);
    guidata(hObject, handles);
end



% Hints: get(hObject,'String') returns contents of edit_datestr as text
%        str2double(get(hObject,'String')) returns contents of edit_datestr as a double

end

% --- Executes during object creation, after setting all properties.
function edit_winlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

handles.tlen	=	[];
guidata(hObject, handles);
end

function edit_winlen_Callback(hObject, eventdata, handles)
% hObject    handle to edit_winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_winlen as text
%        str2double(get(hObject,'String')) returns contents of edit_winlen as a double
tlen	=	str2double(get(hObject,'String'));
if ~(strcmpi(handles.filetype,'psd'))
    tlen_max=	5*60;
    
    if tlen > tlen_max
        %	Make sure user actually intended to set a window this long
        tlen_opts{1}	=	['Given, ' num2str(tlen) 's'];
        tlen_opts{2}	=	['Suggested, ' num2str(tlen_max) 's'];
        tlen_opts{3}	=	['Original, ' num2str(handles.tlen) 's'];
        choice	=	questdlg('Are you sure you want to use a window size this big? Select...', ...
            'Window length sanity check', ...
            tlen_opts{1}, tlen_opts{2}, tlen_opts{3},...
            tlen_opts{3});
        % Handle response
        switch choice
            case tlen_opts{1}
                %tlen	=	tlen;
            case tlen_opts{2}
                tlen	=	tlen_max;
            case tlen_opts{3}
                tlen	=	handles.tlen;
        end
    end
end  %not PSD

set(hObject, 'String', num2str(tlen));
handles.tlen	=	tlen;
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
%contents=get(hObject,'String');
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
%contents=get(hObject,'String');
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
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor')); %#ok<UNRCH>
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
myval	=	get(hObject,'Value');

%%Following commented out because tmax fails to read in WAV data correctly.
%tmin=datenum(get(handles.text_mintime,'String'));
%tmax=datenum(get(handles.text_maxtime,'string'));

tmin	=	handles.tdate_min;
tmax	=	handles.tdate_max;

datenumm	=	tmin + myval*(tmax-tmin);
handles.tdate_start		=	datenumm;
set(handles.edit_datestr,'String',datestr(datenumm,'dd-mmm-yyyy HH:MM:SS.FFF'));

guidata(hObject, handles);
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over slider_datestr.
function	slider_datestr_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to slider_datestr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
climm=str2double(get(hObject,'String'))+[0 diff(climm)];
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
climm(2)=str2double(get(hObject,'String'))+climm(1);
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
handles.filter.f_min	=	str2double(get(hObject,'String'));
handles.filter.changed	=	true;
guidata(hObject, handles);

end
function edit_minfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_minfreq as a double

handles.filter.f_min	=	str2double(get(hObject,'String'));
handles.filter.changed	=	true;
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
handles.filter.f_max	=	str2double(get(hObject,'String'));
guidata(hObject, handles);

end
function edit_maxfreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxfreq as text
%        str2double(get(hObject,'String')) returns contents of edit_maxfreq as a double

handles.filter.f_max	=	str2double(get(hObject,'String'));
handles.filter.changed	=	true;
guidata(hObject, handles);
end
% --- Executes on button press in checkbox_filter.
function checkbox_filter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_filter
handles.audio_stale		=	true;
guidata(hObject,handles);
end
% --- Executes on selection change in popupmenu_scalesound.
function popupmenu_scalesound_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_scalesound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_scalesound contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_scalesound

handles.audio_scale		=	get(hObject,'Value');
handles.audio_stale		=	true;

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_scalesound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_scalesound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.audio_scale		=	1;
guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_playsound.
function pushbutton_playsound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playsound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if	isfield(handles, 'audioplayer')
    player	=	handles.audioplayer;
else
    player	=	[];
end

%	If existing audio player is running, stop it
if	~isempty(player) && ~handles.audio_stale
    ispaused	=	(player.CurrentSample ~= 1) &&...
        (player.CurrentSample ~= player.TotalSamples);
    if	isplaying(player) || ispaused
        stop(player);
        audio_stop(player,[],handles);
        return;
    end
end
%	Only need to recreate player if, file has changed, filter is enabled,
%	or filter has changed.
use_filter	=	get(handles.checkbox_filter,'Value');

audio_stale	=	handles.audio_stale | (use_filter & handles.filter.changed);

if	audio_stale
    Fs	=	handles.Fs;
    x	=	handles.x;
    
    %	Perform filtering if requested
    if	use_filter
        if	handles.filter.changed
            f_max	=	handles.filter.f_max;
            f_min	=	handles.filter.f_min;
            
            %	check frequency ranges
            if	f_min < 0 || f_min >= Fs/2
                f_min	=	0;
            end
            if	f_max <= 0 || f_max > Fs/2
                f_max	=	Fs/2;
            end
            if f_min > f_max
                f_temp	=	f_max;
                f_max	=	f_min;
                f_min	=	f_temp;
                clear	f_temp;
            end
            handles.filter.f_min	=	f_min;
            set(handles.edit_minfreq, 'String', num2str(f_min));
            handles.filter.f_max	=	f_max;
            set(handles.edit_maxfreq, 'String', num2str(f_max));
            
            
            %	Skip filtering if range covers whole band
            if	f_min <=0 && f_max >= Fs/2
                use_filter = false;
            else
                %	default band-pass filter
                df	=	0.1*[f_min f_max]; df = df(df>0);
                df	=	min(df);
                f	=	[f_min-[df 0] f_max+[0 df]];
                a	=	[0 1 0];
                rp	=	0.1;
                rs	=	40;
                dev	= [10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)];
                
                %	low-pass filter
                if	min(f) <= 0
                    f	=	f(end-1:end);
                    a	=	a(end-1:end);
                    dev	=	dev(end-1:end);
                    
                    %	high-pass fiter
                elseif	max(f) >= Fs/2
                    f	=	f(1:2);
                    a	=	a(1:2);
                    dev	=	dev(1:2);
                end
                
                [n,fo,ao,w] =	firpmord(f, a, dev, Fs);
                b			=	firpm(n,fo,ao,w);
                
                handles.filter.b		=	b;
                handles.filter.changed	=	false;
            end
        end
    end
    
    %	Apply filter if needed
    if	use_filter
        b	=	handles.filter.b;
        x	=	filter(b, 1, x);
    end
    
    %	Resample if not within generic sound-card limits
    valid_audio_Fs	=	[5000,8000,11025,22050,44100,48000];
    if ~ismember(Fs, valid_audio_Fs)
        [~, ii]		=	min(abs(valid_audio_Fs - Fs));
        ii	=	ii(1);
        newFs	=	valid_audio_Fs(ii);
        
        disp('Resampling sound');
        x	=	resample(x, newFs, Fs);
        Fs	=	newFs;
    end
    
    %	Floating point data must be scaled to -1:+1 for playback
    x	=	x - mean(x);
    x	=	x ./ max(abs(x));
    
    %	Create audio player object
    Fs_play		=	Fs * handles.audio_scale;
    player		=	audioplayer(x, Fs_play);
    handles.audio_stale		=	false;
    handles.audioplayer		=	player;
    handles.audioFs			=	Fs;
    
    %	Draw initial marker line
    YL		=	ylim;					% get the y-axis limits
    hold on;
    hline	=	plot([0 0], YL, 'r');	% plot the marker
    hold off;
    try	delete(handles.hline);	end
    handles.hline			=	hline;
    
    %	set callbacks for audioplayer
    player.StopFcn = {@audio_stop, handles};
    player.TimerFcn = {@audio_timer, handles}; % timer callback function (defined below)
    player.TimerPeriod = 0.01; % period of the timer in seconds
    
end

%	start playback
play(player);

%	Set controls
set(hObject, 'String', 'Stop');
set(handles.pushbutton_pausesound, 'String', 'Pause');
set(handles.pushbutton_pausesound, 'Enable', 'on');

guidata(hObject, handles);

end

% --- Executes on button press in pushbutton_pausesound.
function pushbutton_pausesound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pausesound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	If existing audio player is running, pause or resume it
if	isfield(handles, 'audioplayer')
    player	=	handles.audioplayer;
    if	isplaying(player)
        pause(player);
        set(hObject, 'String', 'Resume');
    else
        resume(player);
        set(hObject, 'String', 'Pause');
    end
end


end

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;
try
    chan=eval(get(handles.edit_chan,'String'));
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

if ~isempty(strfind(lower(computer),'mac'))
    Islash=strfind(handles.mydir,'/');
else
    Islash=strfind(handles.mydir,'\');
end
chan=get(handles.edit_chan,'String');

save_name=sprintf('soundsamp%s_%s',handles.mydir((Islash(end)+1):end),datestr(tdate_start,30));
disp(['Saving ...' save_name]);

save_path	=	fullfile(handles.outputdir,[ save_name ]);  %AARON: save to local directory, not server

try
    if size(x,2)>size(x,1)
        x=x';
    end
    for Ichan=1:size(x,2)
        xfilt(:,Ichan)=filter(handles.b,1,x(:,Ichan));
    end
    audiowrite(save_path,xfilt/(1.1*max(max(abs(xfilt)))),Fs);
    
catch
    disp('No filtering desired... exists; saving raw acoustic data to WAV');
    xfilt=[];
    audiowrite([save_path '.wav'],(x/(1.1*max(max(abs(x))))),Fs);
    
end
save([save_path '.mat'],'x','xfilt','Fs','tdate_start','hdr');

end

% --- Executes on button press in pushbutton_print.
function pushbutton_print_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(handles.filetype,'MDAT')
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
            
            [~,FF,TT,B] = spectrogram(x(Ichan,:),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
            %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
            %axes(handles.axes1);
            imagesc(TT,FF,10*log10(B));%
            axis('xy')
            fmax=str2double(get(handles.edit_fmax,'String'));
            fmin=str2double(get(handles.edit_fmin,'String'));
            if fmax==0,
                ylim([0 Fs/2]);
                set(handles.edit_fmax,'String',num2str(Fs/2));
            else
                ylim([fmin fmax]*1000);
            end
            %ylim([0 1]);axis('xy')
            climm(1)=str2double(get(handles.edit_mindB,'String'));
            climm(2)=climm(1)+str2double(get(handles.edit_dBspread,'String'));
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
end  %MDAT

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if ~strcmp(handles.display_view,'New Fig')
    figchc=[];
    axes(handles.axes1);
else
    figure(gcf);
    chcc=get(0,'child');
    chcc=[chcc.Number];
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
save_name=sprintf('soundsamp_%s_%s_%s',handles.myfile(1:Idot),datestr(tdate_start,30),chan);
save_path	=	fullfile(handles.outputdir, save_name);
disp(['Printing %s ...' save_name]);

orient landscape
print(gcf,'-djpeg','-r300',sprintf('%s.jpg',save_path));
%print(figchc,'-dtiff',save_path);

%print('-depsc',save_path);

end

% --- Executes on button press in pushbutton_selectpoints.
function pushbutton_selectpoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt={'Number of selections'};
def={'2'};
start_time=handles.tdate_start;
lineNo=1;
answer=inputdlg(prompt,'',lineNo,def);
Np=str2num(answer{1});
tmp=ginput(Np);
msg=[];
for If=1:Np
    msg=sprintf('%s Point %i Time: %8.6f Frequency: %6.2f Hz \n',msg,If,tmp(If,1),1000*tmp(If,2));
end

msg=sprintf('%s\n Absolute time at min time: %s \n' ,msg, datestr(start_time+datenum(0,0,0,0,0,min(tmp(:,1))),0));

duration=max(tmp(:,1))-min(tmp(:,1));
bandwidth=max(tmp(:,2))-min(tmp(:,2));
msg=sprintf('%s Duration: %6.2f sec \n Bandwidth: %6.2f Hz \n Slope: %6.2f Hz/sec\n', ...
    msg,duration,1000*bandwidth,1000*bandwidth/duration);
uiwait(msgbox(msg,'Modal'));

disp(sprintf('Slope is %6.2f Hz/sec',1000*bandwidth/duration))
% disp('Absolute times:')
% disp(datestr(start_time+datenum(0,0,0,0,0,tmp(:,1)),0));
% disp('Relative times and frequencies:');
% disp(tmp(:,1)');
% disp(tmp(:,2)');
% disp('Time elaspsed from first point:');
% disp(tmp(:,1)'-tmp(1,1));


end

% --- Executes on button press in pushbutton_fileread.
function pushbutton_fileread_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fileread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Assumes datenumber is first part of line
tline	=	fgetl(handles.fid);
%Remove leading blanks;
tline	=	fliplr(deblank(fliplr(tline)));
disp(tline);
Idash	=	strfind('-',tline);
Icolon	=	strfind(':',tline);
newtime	=	tline(1:(max(Icolon)+2));


if	tline ~= -1
    tmin	=	handles.tdate_min;
    tmax	=	handles.tdate_max;
    try
        new_date	=	datenum(newtime);
    catch
        uiwait(errordlg('Incorrect datestr'));
        new_date =	[];
    end
    
    if	~isempty(new_date)
        if		new_date < tmin
            new_date	=	tmin;
        elseif	new_date > tmax
            new_date	=	tmax;
        end
        newval	=	(new_date-tmin)/(tmax-tmin);
        handles.tdate_start		=	new_date;
        set(hObject,'String', datestr(new_date,0));
        set(handles.slider_datestr,'Value',newval);
        guidata(hObject, handles);
    end
else
    disp('End of file reached');
    fclose(handles.fid);
end

end

% --- Executes on button press in checkbox_restart.
function checkbox_restart_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_restart
if get(hObject,'Value')==1,
    try
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
fmin=str2double(get(handles.edit_fmin,'String'));
fmax=str2double(get(hObject,'String'));
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

fmin=str2double(get(hObject,'String'));
fmax=str2double(get(handles.edit_fmax,'String'));
ylim([fmin fmax]);

%%%%%%%%%gui_startup_information.m%%%%%%%%%%%%%%%%
% Aaron Thode
% November 22, 2004
% Contains directory locations and names of various files
%   used by the Schultz GUI system

end

% --- Executes on button press in pushbutton_binary.
function pushbutton_binary_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
Ichan=str2double(get(handles.edit_chan,'String'));
[x,~,Fs,tstart]=load_data(handles.filetype,tdate_start,tlen,Ichan,handles);

Nmax=(2^16)-1;
Fs_want=125000;  %Actual playback rate

if Fs<Fs_want
    [q,p]=rat(Fs/Fs_want,0.0001);
    x=resample(x,p,q);
end

y_scale=(x-min(x))*(Nmax/(max(x)-min(x)));
Istart=input('Enter an integer for filename: ');

save_name	=	['PLAY.' int2str(Istart)];
save_path	=	fullfile(handles.notes.folder_name, save_name);

fid=fopen(save_path,'w','ieee-be');
COUNT = fwrite(fid,y_scale,'uint16');
fclose(fid);

end

% --- Executes during object creation, after setting all properties.
function pushbutton_update_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.fig_updated		=	false;
guidata(hObject, handles);
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
handles.audioplayer		=	[];
guidata(hObject, handles);

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

%WARNING! AARON Changes on 11/16 to be compatible
if verLessThan('matlab','8.4.0')
    set(hObject,'SelectionChangeFcn',@uipanel_type_SelectionChangeFcn);
    % execute code for R2014a or earlier
else
    % execute code for R2014b or later
end
%


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

% --- Executes on button press in pushbutton_GSIbearing.
function pushbutton_GSIbearing_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_GSIbearing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get_GSI_bearing(hObject,eventdata,handles);

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

if ~strcmp(handles.mydir(end-6),'S')
    fprintf('Directory %s does not have form S***** \n',handles.mydir);
    return
end
for Idasar=1:NDasar
    %%Assume final directoryname containing files is of the form 'S510G0/'
    
    
    handles.mydir(end-2)=strr(Idasar);
    %%S510G0T20100831T000000.gsi form
    handles.myfile(5)=strr(Idasar);
    fprintf('Directory %s contains %s\n',handles.mydir,handles.myfile);
    try
        handles		=	load_and_display_spectrogram(handles);
        [theta(Idasar),kappa(Idasar),tsec]=get_GSI_bearing(hObject,eventdata,handles);
    catch
        disp('Directory does not exist');
    end
    
    
end

faill=false;
if isfield(handles,'notes')
    try
        locc=handles.notes.Data.Events(2).localization.station_position;
        strr=handles.notes.Data.Events(2).link_names(:,5);
        indicies=double(strr)-64;
        DASAR_coords=zeros(7,2);
        DASAR_coords(indicies,:)=[locc.easting locc.northing];
        VA_cords=[];
    catch
        faill=true;
    end
else
    faill=true;
end

if faill
    prompt1={'File of locations','Site','UTM location to compute range'};
    dlgTitle1='Parameters for GSI localization...';
    %def1={'/Volumes/ThodePortable2/2010_Beaufort_Shell_DASAR/DASAR_locations_2010.mat', '5','[4.174860699660919e+05 7.817274204098196e+06]'};
    %'/Volumes/Data/Shell2010_GSI_Data/DASARlocations/DASAR_locations_2010.mat',
    def1={[filesep fullfile('Volumes','Data','Shell2010_GSI_Data','DASARlocations','DASAR_locations_2010.mat')], '5','[4.174860699660919e+05 7.817274204098196e+06]'};
    answer=inputdlg(prompt1,dlgTitle1,1,def1);
    locs=load(answer{1});
    Isite=str2double(answer{2});
    VA_cords=str2num(answer{3})/1000;
    
    DASAR_coords=[locs.Site{Isite}.easting locs.Site{Isite}.northing];
    
end

Ikeep=find(~isnan(theta)&theta>0);  %Remove absent or unselected DASARs
Igood=find(~isnan(theta)&theta>-2);  %Remove absent DASARS (used to compute center of array)
%theta=theta(Ikeep);
%kappa=kappa(Ikeep);
%Note that the Site contains all DASAR locations,
%  whether they collected data or not

%Using kappa gives too precise an estimate
%[VM,Qhat,~,outcome] = vmmle_r(theta(Ikeep)',DASAR_coords(Ikeep,:),'h',kappa(Ikeep)');
[VM,Qhat,~,outcome] = vmmle_r(theta(Ikeep)',DASAR_coords(Ikeep,:),'h');
mean_coords=mean(DASAR_coords(Igood,:));
CRITVAL=4.60517; %chi2inv(0.90,2);

[~,A,B,ANG,Baxis] = ellipsparms(Qhat,CRITVAL,mean_coords,VM);

figure
Ikeep=find(~isnan(theta(Igood))&theta(Igood)>0);
[~,xg,yg]=plot_location(DASAR_coords(Igood,:),theta(Igood),Ikeep,VM,A,B,ANG);
hold on

if ~isempty(VA_cords)
    tmp=(VA_cords-[xg yg]);
    plot(tmp(1),tmp(2),'o');
    range=sqrt(sum((VA_cords-VM/1000).^2));
    title(sprintf('Range of source from chosen location: %6.2f +/- %3.2f km, minor axis %3.2f km',range,A/1000,B/1000));
end
fprintf('Position of location (can copy for annotation): %s\n',mat2str(VM,10));
yes=menu('Print and Save?','Yes','No');
if yes==1
    orient landscape
    tstart=datestr(datenum(0,0,0,0,0,tsec)+handles.tdate_start,30);
    figure(1);
    print(1,'-djpeg','-r300',sprintf('Localization_G_%s.jpg',tstart));
    print('-djpeg','-r300',sprintf('Spectrogram_G_%s.jpg',tstart));
    save(sprintf('DASAR_localization_%s',tstart),'DASAR_coords','Igood','theta','Ikeep','VM','A','B','ANG','range');
end
close(1)
end

% --- Executes on button press in pushbutton_CSDM.
function pushbutton_CSDM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_CSDM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdate_start=handles.tdate_start;
tlen=handles.tlen;
mydir=pwd;
Ichan='all';  %Hardwire first channel
%set(handles.togglebutton_ChannelBeam,'String','Channel');
[x,t,Fs,~,~,head]=load_data(handles.filetype, tdate_start,tlen,Ichan,handles);

%If not multichannel data, return
if isempty(x)||~head.multichannel
    uiwait(msgbox('CSDM requires multichannel data or non-zero x value'));
    return
end
%Will need to have columns be channels, rows time
if floor(tlen*Fs)~=size(x,1)
    x=x.';
end
%%For MDAT file, channel 1 is already shallowest channel, so remove data below...

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
ovlap=min([1-1/Nfft ovlap]);

fmin=str2double(get(handles.edit_fmin,'String'));
fmax=str2double(get(handles.edit_fmax,'String'));
%chann=input('Enter channel indicies (remember that channels have been flipped, 1 is now shallowest) (1:Nchan):');
%if isempty(chann)
chann=1:head.Nchan; %We now assume that any bad channels are marked in the LOG files.
%end

Ichc=menu('Enter type of frequency processing:','Contour','All (averaged)','Individual Frames (and ray tracing)','Time delay and sum');

delay_and_sum_flag=false;  %Boolean variable is true if delay and sum option selected
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
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);
    
    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);
    
    
    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Threshold (dB)', };
    def1={'-Inf'};
    answer=inputdlg(prompt1,'CSDM threshold?',1,def1);
    
    threshold=eval(answer{1});
    
    [Ksout.Kstot,Ksout.freq,Ksout.VV,Ksout.EE]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,-1,Nfft,0,threshold);
    if ~isempty(Ksout.EE)
        %Ksout.SNR=Ksout.EE(1,:)./Ksout.EE(2,:);
        Ksout.SNR=Ksout.EE(1,:)./sum(Ksout.EE(2:end,:));
    else
        Ksout.SNR=Inf*ones(size(Ksout.freq));
    end
    
    figure
    plot(Ksout.freq,10*log10(Ksout.SNR),'-x');ylabel('dB SNR');xlabel('Freq (Hz)')
    grid on;title('estimated SNR of CSDM frequency components: first eignvalue/sum(rest)')
    
    %%%%%%%%%Extract individual frames and ray trace%%%%%%%%%
elseif Ichc==3 %%Extract frames and ray trace
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);
    
    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);
    
    
    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Threshold (dB)', 'optional FFT for longer FFT of small segments'};
    def1={'-Inf',''};
    answer=inputdlg(prompt1,'CSDM threshold?',1,def1);
    
    threshold=eval(answer{1});
    if isempty(answer{2})
        Nfft2=Nfft;
    else
        Nfft2=eval(answer{2});
    end
    
    [Ksout.Kstot,Ksout.freq,Ksout.t]=extractKsframes(x,ovlap,Nfft2,chann,frange,Fs,Nfft,0,threshold);
    
elseif Ichc==4 %% simple time and delay on filtered signal...
    
    delay_and_sum_flag=true;
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);
    
    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);
    
    
    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Min Freq (Hz)', 'Max Freq (Hz)'};
    def1={num2str(frange(1)),num2str(frange(2))};
    answer=inputdlg(prompt1,'Check frequency selection',1,def1);
    
    frange=[eval(answer{1}) eval(answer{2})];
    
    bpFilt = designfilt('bandpassiir', 'FilterOrder', 20, ...
             'HalfPowerFrequency1', frange(1), 'HalfPowerFrequency2', frange(2),...
             'SampleRate', Fs);
  
    %fvtool(bpFilt);
    x=filter(bpFilt,x);
    
    
end


if ~delay_and_sum_flag
    yes=menu('Beamform?','No','Conventional','MV','Both','Reflection Coefficient Estimation','MFP');
else
    yes=2;
end
yes_eigen=0;
while yes>1
    
    if yes<6
        
        prompt1={'Vector of angles (deg) [(-20:0.2:20)]','Vector of sin angles (deg) [(-0.2:0.005:0.2)]','hydrophone indicies [all]','sound speed (m/sec)', ...
            };
        def1={'-10:0.2:10', '-0.3:0.005:0.3',sprintf('[1:%i]',length(chann)),'1480'};
        
        answer=inputdlg(prompt1,'Beamforming parameters',1,def1);
        try
            if isempty(answer{2})
                angles=eval(answer{1});
            else
                angles=180*asin(eval(answer{2}))/pi;
                fprintf('Angular range is %6.2f to %6.2f degrees\n',min(angles),max(angles));
            end
            Igood_el=eval(answer{3});
            cc=eval(answer{4});
            
        catch
            errdlg('Could not understand your beamforming parameters');
            return
            
        end
    end
    switch yes
        case 2
            beam_str='CV';
            
            if delay_and_sum_flag
                beamchc='delaynsum';
            else
                Ndim=ndims(Ksout.Kstot);
                
                if Ndim==3
                    beamchc='freqvspwr';
                elseif Ndim==4
                    beamchc='angvstime';  %migration
                end
            end
            if strcmp(beamchc,'freqvspwr')
                B=conventional_beamforming(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
                
                figure(20);
                imagesc(Ksout.freq,[],10*log10(abs(Ksout.EE)));xlabel('frequency (Hz)');ylabel('Eigenvalue');colorbar
                
                Ieig=input('Pick eigenvector [return yields none]:');
                if ~isempty(Ieig)
                    V1=squeeze(Ksout.VV(:,Ieig,:));
                    for If=1:length(Ksout.freq)
                        Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';
                    end
                    B_eig=conventional_beamforming(Ksout.Kstot_eig(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
                    yes_eigen=1;
                else
                    close(20);
                end
                
                ray_angles=plot_beamforming_results(true); % peakpicking
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%Individual snapshots and ray tracing%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif strcmp(beamchc,'angvstime') %Individual snapshots
                Bsum=zeros(length(angles),length(Ksout.t));
                df=Ksout.freq(2)-Ksout.freq(1);
                %plot(Bsum,angles,'k');
                
                for Isnap=1:length(Ksout.t)
                    if rem(Isnap,10)==0,disp(Isnap);end
                    B=conventional_beamforming(squeeze(Ksout.Kstot(Igood_el,Igood_el,:,Isnap)),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
                    %Bsum(:,Isnap)=10*log10(sum(abs(B)))/length(Ksout.freq);
                    Bsum(:,Isnap)=(sum(10*log10(abs(B))))/length(Ksout.freq);
                end
                
                %plot vs sin angle
                
                figure
                imagesc((Ksout.t-min(Ksout.t))*1000,sin(angles*pi/180),Bsum)
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('sine of Elevation angle');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                set(gca,'xtick',0:5:xlimm(2));
                
                %plot migration angle...
                figure(21)
                imagesc((Ksout.t-min(Ksout.t))*1000,angles,Bsum)
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('Elevation angle (deg)');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                set(gca,'xtick',0:5:xlimm(2));
                printstr=sprintf('AngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
                saveas(gcf,[printstr '.fig'],'fig')
                
                
                
                
                %%Permit ray tracing if desired
                prep_ray_trace;
            
            
            elseif strcmp(beamchc,'delaynsum') %Individual snapshots
            
           
                Bsum=zeros(length(angles),size(x,1));
                space=head.geom.rd(Igood_el);
                space=space(2)-space(1);
                for Isnap=1:length(angles)
                    if rem(Isnap,10)==0,disp(Isnap);end
                    xtot=delaynsum(x,-angles(Isnap),space,Fs,Igood_el,cc);
                   % Bsum(Isnap,:)=abs(hilbert(xtot'));
                    Bsum(Isnap,:)=xtot';
                end
                
                Bsum_org=Bsum;
                Bsum=abs(hilbert(Bsum.'))';
                Bsum=Bsum/max(max(Bsum));
                Bsum=20*log10(Bsum);
                
                %plot vs sin angle
                tt=(1:length(xtot))/Fs;
                figure(20);clf;
                imagesc(tt*1000,sin(angles*pi/180),Bsum);
                caxis([-20 0]);colorbar
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('sine of Elevation angle');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                set(gca,'xtick',0:5:xlimm(2));
                printstr=sprintf('SinAngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
                saveas(gcf,[printstr '.fig'],'fig')
                
                %plot migration angle...
                figure(21);clf;
                imagesc(tt*1000,angles,Bsum)
                caxis([-20 0]);colorbar
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('Elevation angle (deg)');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                if diff(xlimm)<100
                    set(gca,'xtick',0:5:xlimm(2));
                else
                    set(gca,'xtick',0:100:xlimm(2));
                    
                end
                printstr=sprintf('AngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
                saveas(gcf,[printstr '.fig'],'fig')
                
                %Repeat using cross-correlation result (to remove
                %   complexities in time waveform)
                figure(21)
                uiwait(msgbox('Please select a box to cross-correlate with other beams (matched filter), hit return to skip:'));
                tmp=ginput(2);
               % tmp=[];
                if isempty(tmp)
                    prep_ray_trace;
                    return
                end
                
                Bsum_corr=zeros(length(angles),size(Bsum,2));
                angle_want=mean(tmp(:,2));
                indd=round(tmp(1,1)*Fs):round(Fs*tmp(2,1));
                [~,Iang]=min(abs(angle_want-angles));
                x_match=fliplr((Bsum_org(Iang,:)));Nx=length(x_match);
                
                for Iang=1:length(angles)
                   Bsum_corr(Iang,:)= conv(x_match,Bsum_org(Iang,:),'same');
                    
                end
                Bsum_corr=20*log10(abs(hilbert(Bsum_corr.')))';
                Bsum_corr=Bsum_corr-max(max(Bsum_corr));
                %Bsum_conv=Bsum_conv./max(max(Bsum_conv));
                %plot vs sin angle
                 figure(22)
                imagesc(1000*(1:size(Bsum_corr,2))/Fs,(angles),Bsum_corr);
                caxis([-20 0]);colorbar
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('sine of Elevation angle');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                 if diff(xlimm)<100
                    set(gca,'xtick',0:5:xlimm(2));
                else
                    set(gca,'xtick',0:100:xlimm(2));
                    
                end
                printstr=sprintf('MatchSinAngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
               % print(gcf,'-djpeg',[printstr '.jpg']);
               % saveas(gcf,[printstr '.fig'],'fig')
               %%Permit ray tracing if desired
               prep_ray_trace;
               
            return
            end
            
        case 3
            beam_str='MV';
            
            B=MV_beamforming(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
        case 4
            beam_str='CVnMV';
            
            B=conventional_beamforming(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
            
            B2=MV_beamforming(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
        case 5
            beam_str='RC';
            
            R=derive_reflection_coefficient2(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
        case 6
            beam_str='MFP';
            
            prompt1={'Model file..','tilt offset between top and bottom phone (m)','ranges (m)', 'depths (m):','plot intermediate images?', ...
                'Nfft for plotting:','Frequency SNR_cutoff (dB), or array with specific frequencies to process..','Primary Eigenvector only?','receiver depths (m):'};
            
            dlgTitle1='Parameters for matched field processing...';
            
            [MFP_replica_file,range_str,default_tilt,default_SNR]=load_MFP_scenario;
            
            if length(MFP_replica_file)>1
                error('Cannot select multiple environmental models during MFP');
            end
            
            %head.geom.rd=get_VA_2010_depths_and_offset(2010,handles.myfile,MFP_replica_file{1});
            head.geom.rd=head.geom.rd;
            
            def1={MFP_replica_file{1}, default_tilt{1},range_str{1},'1:1:55','0','2048',['[' num2str(default_SNR{1}) ']'],'yes',num2str(head.geom.rd)};
            answer=inputdlg(prompt1,dlgTitle1,1,def1);
            model_name=answer{1};
            tilt_offset=str2num(answer{2});
            ranges=str2num(answer{3});
            depths=str2num(answer{4});
            plot_chc=str2num(answer{5});
            Nfft2=str2num(answer{6});
            SNRmin=str2num(answer{7});  %May also be frequency
            head.geom.rd=str2num(answer{9});
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
                save_result=matched_field_processor(model_name,tilt_offset,ranges,depths,Ksout.Kstot_eig,Ksout.freq,10*log10(Ksout.SNR),head.geom.rd,SNRmin,srcname);
            else
                save_result=matched_field_processor(model_name,tilt_offset,ranges,depths,Ksout.Kstot,Ksout.freq,10*log10(Ksout.SNR),head.geom.rd,SNRmin,srcname);
            end
            
            
            if ~isempty(save_result)
                save(srcname,'Ksout','head','tdate_start','tlen');
                write_covmat(Ksout.freq,Ksout.Kstot,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
                if ~isempty(Ksout.VV)
                    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
                    write_covmat(Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
                end
            end
            return
    end
    
    plot_beamforming_results(true); %ppeak picking
    
    
    yes=menu('Beamform?','No','Conventional','MV','Both','Reflection Coefficient Estimation','MFP');
    
end  %while yes

yes=menu('Save CSDM?','Yes','No');
if yes==1
    
    Ksout.Nfft=Nfft;
    Ksout.ovlap=ovlap;
    V1=squeeze(Ksout.VV(:,1,:));
    for If=1:length(Ksout.freq)
        Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';
        
    end
    
    
    %     srcname=sprintf('CSDM_Nfft%i_%s_flipped',Nfft,datestr(tdate_start,30));
    %     write_covmat(Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
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
    srcname=sprintf('CSDM_Nfft%i_%s',Nfft,datestr(tdate_start,30));
    save(srcname,'Ksout','head','tdate_start','tlen');
    write_covmat(Ksout.freq,Ksout.Kstot,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    
    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
    write_covmat(Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    
end

%%Inner function for prepping a ray trace and/or invariant trace
    function prep_ray_trace(ray_angles2,dt_data)
        
        repeat=1;
        while repeat==1
        ButtonName = questdlg('What to do?', 'Migration Analysis','Ray Tracing', 'Invariant Tracing', 'Nada', 'Nada');
        switch ButtonName,
            case 'Nada'
                return
            case 'Ray Tracing'
                scenarios=raytrace_MATLAB;
                [Selection,OK]=listdlg('liststring',scenarios,'SelectionMode','single');
                if ~OK
                    return
                end
                
                %figure(21);pause(0.5);
                %%Interpret results depending on plot
                if strcmp(beamchc,'angvstime')||strcmp(beamchc,'delaynsum')
                    
                    prompt1={'Number of paths to pick?'};
                    def1={'0'};
                    %Other defaults: MURI_Feb2014_20140218T140039,
                    answer=inputdlg(prompt1,'Ray picks',1,def1);
                    Nrays=str2num(answer{1});
                    
                    tmp=ginput(Nrays);
                    
                    
                    dt=max(tmp(:,1))-tmp(:,1);  %Rays that arrive first will have positive numbers
                    [dt_meas,Isort]=sort(dt);  %Arrange by first-arriving rays
                    ang_meas=tmp(Isort,2);
                elseif strcmp(beamchc,'freqvspwr')
      
                        ang_meas=(ray_angles2);
                        dt_meas=[dt_data];  %No relative arrival information
                    
                        dt_meas=dt_meas-min(dt_meas);
                        [dt_meas,Isort]=sort(dt_meas);
                        ang_meas=ang_meas(Isort);
                end
                
                [x,dt,dz]=raytrace_MATLAB(ang_meas,dt_meas,scenarios{Selection});
            case 'Invariant Tracing'
                
                prompt1={'Ranges (km)','Beta','sound speed at receiver (m/sec)'};
                
                dlgTitle1='Parameters for invariant tracing...';
                
                def1={'[40:10:60 90:10:120]','-9','1481'};
                answer=inputdlg(prompt1,dlgTitle1,1,def1);
                
                try
                    figure(20);
                    pause(0.5);
                    invariant_range_estimate(eval(answer{1}),eval(answer{2}),eval(answer{3}));
                catch
                    uiwait(errordlg('Invariant Processing failed!')); 
                end
                
                
        end % switch
        
        repeat=menu('Another analysis?','Yes','No');
        end
        
        
    end

%%Inner function for plot_beamforming_results
    function ray_angles=plot_beamforming_results(peak_picking)
        ray_angles=[];
        %Plot beamforming output
        try
            close(1);
        end
        figure(1);clf
        if yes==4||yes_eigen*yes==2
            subplot(2,1,1)
        end
        
        %Image of beamform output vs look angle and frequency.
        imagesc(Ksout.freq/1000,angles,10*log10(B'));
        colorbar
        cmap=colormap;
        caxis([60 100])
        caxis('auto');
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Frequency (kHz)');ylabel('Angle from horizontal (deg)');grid on;
        title(sprintf('%s, %i FFT, %i elements',datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
        set(gcf,'colormap',cmap(1:4:64,:));
        
        figure(2);clf
        if yes==4||yes_eigen*yes==2
            subplot(2,1,1)
        end
        
        %%%Plot summed beampattern
        df=Ksout.freq(2)-Ksout.freq(1);
        Bsum=sum(10*log10((B)))/length(Ksout.freq);
        plot(Bsum,angles,'k');
        set(gca,'fontweight','bold','fontsize',14);axis('ij');
        xlabel('Mean dB Beampower ');ylabel('Angle from horizontal (deg)');grid on;
        title(sprintf('%s, %i FFT, %i elements',datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
        
        if ~peak_picking
            return
        end
        %%%%%%%%%%%
        %%%Look for local maxima:
        not_happy=true;
        hplot=[];
        while not_happy
            if ~isempty(hplot)
                set(hplot,'vis','off');
            end
            prompt1={'Automatic peak pick?','Half-beamwidth (deg)', 'threshold SNR'};
            def1={'no','3', '3'};
            
            answer=inputdlg(prompt1,'Peakpicking parameters',1,def1);
            if isempty(answer)
                return
            end
            try
                peak_pick_chc=strcmpi(answer{1},'yes');
                halfbeamwidth=eval(answer{2});
                thresholddB=eval(answer{3});
                
            catch
                errdlg('Could not understand your peakpicking parameters');
                return
                
            end
            
            if peak_pick_chc
                peaks=peak_picker_Thode(Bsum,angles,halfbeamwidth,[min(angles) max(angles)],thresholddB);
                figure(2)
                hold on
                hplot=plot(peaks{1}.adp.PdB,peaks{1}.adp.F,'go');
                for I=1:length(peaks{1}.adp.PdB)
                    fprintf('Adp Path %i: %6.2f deg beampower %6.2f dB\n',I,peaks{1}.adp.F(I),peaks{1}.adp.PdB(I));
                    %try
                    %fprintf('ISI Path %i: %6.2f deg beampower %6.2f dB\n',I,peaks{1}.isi.F(I),peaks{1}.isi.PdB(I));
                    %end
                end
                [PdB,Isort]=sort(peaks{1}.adp.PdB,'descend');
                ray_angles=peaks{1}.adp.F(Isort);
            else
                prompt1={'Number of picks'};
                def1={'3'};
                try
                    answer=inputdlg(prompt1,'Peakpicking parameters',1,def1);
                    
                    tmp=ginput(eval(answer{1}));
                catch
                    return
                end
                PdB=tmp(:,1);
                ray_angles=tmp(:,2);
                
                %Find local maxima near each selection
                fprintf('Original angles selected: %s\n',mat2str(ray_angles,3));
                sigma_ang=2;
                for Iray=1:length(ray_angles)
                    Iseg=find(angles>ray_angles(Iray)-sigma_ang&angles<ray_angles(Iray)+sigma_ang);
                    [~,Imax]=max(Bsum(Iseg));
                    ray_angles(Iray)=angles(Iseg(Imax));
                    
                end
                fprintf('Revised angles: %s\n',mat2str(ray_angles,3));
                
                %[PdB,Isort]=sort(PdB,'descend');
                %ray_angles=ray_angles(Isort);
            end
            not_happy = questdlg('Redo Peak Pick?', ...
                ' Question', ...
                'Yes', 'No', 'No');
            not_happy=strcmp(not_happy,'Yes');
        end %while not_happy
        %%save ray angles, sorted by power
        
        if length(ray_angles)>1
            [Bshift,dt_data]=conventional_beam_timedelay(Ksout.Kstot(Igood_el,Igood_el,:),ray_angles,Ksout.freq,head.geom.rd(Igood_el),cc,Nfft,Fs);
            prep_ray_trace(ray_angles,dt_data);
        end
        
        
        if yes==4
            figure(1)
            subplot(2,1,2);
            imagesc(Ksout.freq,angles,10*log10(abs(B2)'));
            colorbar
            %cmap=colormap;
            caxis([60 100])
            caxis('auto');
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');ylabel('Angle from horizontal (deg)');grid on;
            title('MV processor');
            %set(gcf,'colormap',cmap(1:4:64,:));
            
            figure(2);
            subplot(2,1,2)
            plot(angles,sum(10*log10((B2)))/length(Ksout.freq),'k');
            set(gca,'fontweight','bold','fontsize',14);
            ylabel('Mean Beampower (dB)');xlabel('Angle from horizontal (deg)');grid on;
            title(sprintf(' %s, %i FFT, %i elements',datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
            
        elseif yes_eigen*yes==2
            figure(1)
            subplot(2,1,2)
            imagesc(Ksout.freq,angles,10*log10(B_eig'));
            colorbar
            %cmap=colormap;
            caxis([60 100])
            caxis('auto');
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');ylabel('Angle from horizontal (deg)');grid on;
            title(sprintf('Eigenvector %i only: %s, %i FFT, %i elements',Ieig,datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
            %set(gcf,'colormap',cmap(1:4:64,:));
            
            figure(2);
            subplot(2,1,2)
            plot(angles,sum(10*log10((B_eig)))/length(Ksout.freq),'k');
            set(gca,'fontweight','bold','fontsize',14);
            ylabel('Mean dB Beampower ');xlabel('Angle from horizontal (deg)');grid on;
            title(sprintf('Eigenvector %i only: %s, %i FFT, %i elements',Ieig,datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
            
        end
        
        yes_print=menu('Print?','Yes','No');
        if yes_print==1
            for III=1:3
                figure(III)
                orient tall
                
                 ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                printstr=sprintf('Beamforming_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
               
                print(gcf,'-djpeg','-r300',sprintf('Beamforming%s_%s_%ito%ideg_%4.2fres_%i.jpg', ...
                    beam_str,datestr(ttt,'yyyymmddTHHMMSS.FFF'),min(angles),max(angles),angles(2)-angles(1),III));
            end
        end
    end %function plot_beamforming_results
end

% --- Executes on button press in pushbutton_Mode.
function pushbutton_Mode_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%Only called if MDAT exists...

tmp=ginput(1);
twant=tmp(1);
%fwant=tmp(2)*1000; %in Hz


tdate_start=handles.tdate_start;
tlen=handles.tlen;
%for Ichan=1:8
[x,~,Fs,tstart,junk,head]=load_data(handles.filetype, tdate_start,tlen,'all',handles); %#ok<*ASGLU>

%If not multichannel data, return
if isempty(x)||~head.multichannel
    uiwait(msgbox('Mode filtering requires multichannel data or non-zero x value'));
    return
end
%AARON: removed lines below because multichannel data will always get
%flipped
%if size(x,2)>1
%    x=x';
%end
%xall{Ichan}=x;
%end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
ovlap=min([1-1/Nfft ovlap]);

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');
for Ichan=1:head.Nchan
    
    %%Test 1
    [S,FF,TT,B] = spectrogram(x(:,Ichan)-mean(x(:,Ichan)),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
    [junk,Iwant]=min(abs(twant-TT));
    Bslice(:,Ichan)=B(:,Iwant); %#ok<*AGROW>
    
    
    %Sign test
    index=round(twant*Fs)+(1:Nfft);
    Bslice2(:,Ichan)=fft((x(index,Ichan)-mean(x(index,Ichan))).*(hanning(Nfft)));
    
end
Bslice2=Bslice2((1:(Nfft/2+1)),:);
figure(1)
for I=1:2
    subplot(2,1,I)
    if I==1
        contourf(FF,head.geom.rd,10*log10(abs(Bslice2')), 'linestyle','none');axis('ij')
        %caxis([80 120])
        title(sprintf('log magnitude Raw FFT %s: %s plus %6.2f sec',get(handles.text_filename,'String'),datestr(tdate_start),twant));
    else
        imagesc(FF,head.geom.rd,abs(Bslice2'));  title('magnitude of FFT');
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
        contourf(FF,head.geom.rd,10*log10(Bslice'), 'linestyle','none');axis('ij')
        caxis([80 120])
        title(sprintf('Log spectrogram %s: %s plus %6.2f sec',get(handles.text_filename,'String'),datestr(tdate_start),twant));
    else
        imagesc(FF,head.geom.rd,(Bslice'));title('Spectrogram slice');
        %contourf(FF,head.geom.rd,(Bslice'),1e10:1e10:10e10,'linestyle','none');
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
plot(head.geom.rd,mode/max(abs(mode)),'-x');grid on;title(sprintf('Depth: %6.2f Frequency: %6.2f', head.D,freq))

subplot(3,1,2)
angg=180*(unwrap(angle(Bslice2(Ibest,:))))/pi;

plot(head.geom.rd,angg,'-x');grid on;
title('unwrapped angle in degrees');

subplot(3,1,3)
angg=diff(unwrap(angle(Bslice2(Ibest,:))));
angg=[0 angg];
%%since angg=kL =2*pi*f*L/c in case of tilt
tilt_offset=angg*1490./(2*pi*freq);
plot(head.geom.rd,cumsum(tilt_offset),'-x');grid on;
title('horizontal meter offset from top element');

yes=input('Save a mode example?','s');
if ~isempty(yes)
    modeno=input('Enter mode number:');
    save(deblank(sprintf('Mode%iExample%3.0fHz_%s',modeno,FF(Ibest),datestr(tdate_start+datenum(0,0,0,0,0,twant),30))),'mode','head','freq','angg','tilt_offset');
end

close(1:3)


end

% --- Executes on button press in pushbutton_tilt.
function pushbutton_tilt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tdate_start=handles.tdate_start;
tlen=handles.tlen;
%for Ichan=1:8
[x,t,Fs,tstart,junk,head]=load_data(handles.filetype,tdate_start,tlen,'all',handles);

%If not multichannel data, return
if isempty(x)||~head.multichannel
    uiwait(msgbox('Tilt estimate requires multichannel data or non-zero x value'));
    return
end
%x(1) is shallowest element...

%AARON: make sure channels are columns
if size(x,2)~=head.Nchan
    %if size(x,2)>1
    x=x';
end

%xall{Ichan}=x;
%end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
ovlap=min([1-1/Nfft ovlap]);

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
plot(cumsum(offsets),head.geom.rd,'-x',offsets_cum,head.geom.rd,'-ro');
set(gca,'fontweight','bold','fontsize',14);
grid on;ylabel('receiver depth (m)');xlabel('tilt offset (m)');
legend('offset accumulated from adjacent elements','direct comparison to shallowest element','Location','best');
xlim([-5 5]);
hold on;axis('ij')
title(sprintf('Depth: %6.2f m Times: %s, freq range: %i to %i Hz', head.geom.D,datestr(tdate_start,30),round(minfreq),round(maxfreq)));

yes=menu('Save a tilt example?','Yes','No');
if yes==1
    figure(1)
    orient landscape
    print(1,'-djpeg','-r300',sprintf('TiltExample_%s_%ito%iHz',datestr(tdate_start,30),round(minfreq),round(maxfreq)));
end
%close(1);
%%%Test to understand xcorr lags
%x1=[0 0 0 1 0 0 0 0 0 0];
%x2=[0 0 0 0 0 0 1 0 0 0];  %x2 lags x1 (signal arrives later on x2)
%[cc,lags]=xcov(x1,x2); %lags max at -4

%Thus if first vector leads second, lags are negative
%  If first vector lags second, lags are positive

end

function pushbutton_modalfiltering_Callback(hObject, eventdata, handles)
% Estimate distance and range by isolating modes on a vertical arrray...
% handles    structure with handles and user data (see GUIDATA)

figss=get(0,'child');
Iclose=find(figss-floor(figss)==0);
if ~isempty(figss(Iclose))
    close(figss(Iclose));
end
tdate_start=handles.tdate_start;
tlen=handles.tlen;
%for Ichan=1:8
[data.x,t,Fs,tstart,junk,head]=load_data(handles.filetype, tdate_start,tlen,'all',handles);

%If not multichannel data, return
if isempty(x)||~head.multichannel
    uiwait(msgbox('Modal filtering requires multichannel data or non-zero x value'));
    return
end
%xall{Ichan}=x;
%end

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
ovlap=min([1-1/Nfft ovlap]);


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
        Islash=max(strfind(MFP_scenario{J},'/'))+1;
        fprintf('tilt: %s File: %s \n',default_tilt{J},MFP_scenario{J}(Islash:end));
    end
    dlgTitle1='Parameters for modal beamforming...';
    def1={range_str{1},'0.1','0.025','0','2048','1:1:55',default_tilt{1},'3'};
    answer=inputdlg(prompt1,dlgTitle1,1,def1);
    range_guess=str2num(answer{1});
    maxdelay(1)=str2double(answer{2});
    maxdelay(2)=maxdelay(1);
    maxdelay(3)=str2double(answer{3});
    plot_chc=str2double(answer{4});
    Nfft2=str2double(answer{5});
    depths=str2num(answer{6});
    tilt_offset=str2num(answer{7});
    Maxmodes=str2double(answer{8});
    
    
end


for Imodel=1:length(MFP_scenario)
    
    if length(MFP_scenario)==1
        
        prompt1={'Model file..','tilt offsets (deg)','range guesses (m)', ...
            'maxdelay for mode 1:2 range estimation (s):', ...
            'maxdelay for mode 1:N range estimation (s):', ...
            'plot intermediate images?', ...
            'Nfft for plotting:','depths','max modes'};
        
        dlgTitle1='Parameters for modal beamforming...';
        
        def1={MFP_scenario{Imodel}, default_tilt{Imodel},range_str{Imodel},'0.1','0.025','0','2048','1:0.5:55','2'};
        
        answer=inputdlg(prompt1,dlgTitle1,1,def1);
        model=load(answer{1});
        disp(['Loaded ',answer{1}]);
        tilt_offset=str2num(answer{2});
        range_guess=str2num(answer{3});
        maxdelay(1)=str2double(answer{4});
        maxdelay(2)=maxdelay(1);
        maxdelay(3)=str2double(answer{5});
        plot_chc=str2double(answer{6});
        Nfft2=str2double(answer{7});
        depths=str2num(answer{8});
        Maxmodes=str2double(answer{9});
    else
        
        model=load(MFP_scenario{Imodel});
        disp(sprintf('Loaded %s',MFP_scenario{Imodel})); %#ok<*DSPS>
        %tilt_offset=str2double(default_tilt{Imodel});
        
    end
    maxdelay(4:Maxmodes)=maxdelay(3);
    
    %disp('Save signal sample if you wish...');
    % keyboard
    %following operation is now done in load_data...
    %data.x=flipud(data.x); %Shallowest data now first...
    
    %%AARON:  REMOVE THIS!
    %head.geom.rd=head.geom.rd+2;
    %disp('WARNING WARNING WARNING I added 2 m depth to array elements!!!!!!');
    
    Igood	=	zeros(1,length(head.geom.rd));
    for I = 1:length(head.geom.rd)
        [junk, Igood(I)]		=	min(abs(model.rd-head.geom.rd(I)));
    end
    
    %%depths are desired modeled source depths.
    Igood_source=zeros(1,length(depths));
    for I=1:length(depths) %#ok<*FXUP>
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
                % frequencies.  This makes a convenient time shift..
                
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
        %compute_depth_estimate;
        
    end %Itilt
    %%Plot summary of range estimation time lags...
    plot_range_estimates;
    
    if plot_chc==1
        if exist('MM', 'var');
            movie2avi(MM,'tiltTest');
        end
        keyboard
        for I=1:12
            if any(get(0,'child')==I)
                figure(I)
                orient tall
                print(I,'-djpeg','-r300',sprintf('ModalBeamform%i',I));
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
                    
                    
                    set(h1,'markeredgecolor',[plotstr_tilt(Itilt2) ]) %#ok<*NBRAK>
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
                    print(gcf,'-djpeg','-r300',sprintf('RangeEstimate_%s.jpg',savename_base));
                    
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
                print(gcf,'-djpeg','-r300',sprintf('RangeError1_%s.jpg',savename_base));
                
                figure(2)
                orient landscape
                legend('full model','Pekeris model','Averaged Pekeris')
                print(gcf,'-djpeg','-r300',sprintf('RangeError2_%s.jpg',savename_base));
                
                
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
            print(gcf,'-djpeg','-r300',sprintf('DepthEstimate2D_%s.jpg',savename_base));
            
            
            %Plot frequency-averaged depth estimator%
            
            figure(100+Itilt+Ir_best);
            plot1=10*log10(sum(amb{Imodel}'))/length(Iwant); %#ok<*UDIM>
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
            print(gcf,'-djpeg','-r300',sprintf('DepthEstimate_%s.jpg',savename_base));
            
            
            
            
            
        end  %Inner function plot_depth_estimates
        
    end  %inner function compute_depth_estimates
%%%%%%%%%%%%%%%%%


%%Inner function for plotting
    function tltchk=plot_mode_filtering
        persistent timmes_signal timmes_noise pwrr yes %#ok<NUSED>
        
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
            
            for Imm=1:Maxmodes
                [Bbeam,FF1,TT] = specgram(yout(:,Imm),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
                figure(10*Itilt+1)
                set(gcf,'units','norm','pos',[ 0.1589    0.3345    0.2216    0.6102]);
                
                subplot(Maxmodes,1,Imm);
                imagesc(TT,FF1,20*log10(abs(Bbeam).*normm));%
                orient landscape
                grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
                axis('xy');ylim([minfreq maxfreq]);
                caxis([70 120]-30);
                %caxis([-150 -70])
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
            
            
            
            for Imm=1:Maxmodes
                [Bbeam1,FF1,TT] = specgram(yout_corrected(:,Imm),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
                figure(2)
                subplot(Maxmodes,1,Imm);
                imagesc(TT,FF1,20*log10(abs(Bbeam1).*normm));%
                orient landscape
                grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
                ylabel('Hz'); xlabel('Time (s)');
                
                axis('xy');ylim([minfreq maxfreq]);
                caxis([70 120]-30);grid on
                %caxis([-150 -70])
                if Imm==1
                    %title(sprintf('Mode %i, modal beamforming and dispersion correction: vertical tilt: %6.2f deg',Imm, tilt));
                end
                if Imm<Maxmodes
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
            [~,Imax]=max(abs(tmp));
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
        subplot(Maxmodes-1,1,Im-1)
        plot(lags/Fs,xc);xlim(10*[-0.05 0.05]);
        grid on
        title(sprintf('Mode %i, leads mode 1 by %10.6f ms',Im,1000*talign{Imodel}(Itilt,Ir,Im)));
        xlabel('Relative delay (s)');
        ylabel('Cross-correlation');
        %end
    end
end  %function modal_filtering

% --- Executes on button press in checkbox_teager.
function checkbox_teager_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_teager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_teager
end

% --- Executes on button press in pushbutton_next.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	increment start date
new_date		=	handles.tdate_start;
tlen			=	handles.tlen;
tlen			=	tlen/60/60/24;	%	convert to fractional days
new_date		=	new_date + tlen;
%	Trigger callback to update/check other GUI elements
set(handles.edit_datestr, 'String', datestr(new_date,'dd-mmm-yyyy HH:MM:SS.FFF'));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

pushbutton_update_Callback(handles.pushbutton_update,eventdata,handles);

end

% --- Executes on button press in pushbutton_prev.
function pushbutton_prev_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	decrement start date
new_date		=	handles.tdate_start;
tlen			=	handles.tlen;
tlen			=	tlen/60/60/24;	%	convert to fractional days
new_date		=	new_date - tlen;
%	Trigger callback to update/check other GUI elements
set(handles.edit_datestr, 'String', datestr(new_date,'dd-mmm-yyyy HH:MM:SS.FFF'));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

pushbutton_update_Callback(handles.pushbutton_update,eventdata,handles);

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Beginning of annotation GUI subroutines %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	new annotation stuff

% --- Executes on button press in pushbutton_notes_last.
function pushbutton_notes_last_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N		=	length(handles.notes.Data.Events);
if	(N == 0)
    error('What happened?');
end

handles.notes.i_sel		=	N;

handles		=	update_events(handles);

%	this might be redundant
guidata(hObject, handles);

end

% --- Executes on button press in pushbutton_notes_first.
function pushbutton_notes_first_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N		=	length(handles.notes.Data.Events);
if	(N == 0)
    error('What happened?');
end

handles.notes.i_sel		=	1;

handles		=	update_events(handles);

%	this might be redundant
guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_notes_next.
function pushbutton_notes_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i_sel	=	handles.notes.i_sel;
N		=	length(handles.notes.Data.Events);

if	(N == 0)
    error('What happened?');
end

if	isempty(i_sel)
    i_sel	=	1;
else
    i_sel	=	i_sel + 1;
end

if	i_sel > N
    i_sel	=	1;
end
handles.notes.i_sel		=	i_sel;

handles		=	update_events(handles);

%	this might be redundant
guidata(hObject, handles);

end

% --- Executes on button press in pushbutton_notes_prev.
function pushbutton_notes_prev_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i_sel	=	handles.notes.i_sel;
N		=	length(handles.notes.Data.Events);

if	isempty(i_sel) || (i_sel == 0) || (N == 0)
    warning('Selection index is not valid');
    return;
end

i_sel	=	i_sel - 1;
if	i_sel < 1
    i_sel	=	N;
end
handles.notes.i_sel		=	i_sel;

handles		=	update_events(handles);

%	this might be redundant
guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function pushbutton_notes_save_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.notes.saved	=	true;
guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_notes_save.
function pushbutton_notes_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	Shouldn't be able to click this button anyway, but can't hurt to check
readonly	=	handles.notes.readonly;
if readonly
    warning('Readonly flag set, file not saved');
    return;
end

%	The following should never be empty, but...
Data		=	handles.notes.Data;
if	isempty(Data)
    warning('Notes data empty, nothing to save');
    return;
end

file_path	=	handles.notes.file_path;
if	isempty(file_path)
    warning('Notes file/folder not set, file not saved');
    return;
end

%	Save whole data structure, allows implicit expansion of named variables
%	Only save if changed
if	~handles.notes.saved
    GUI_params		=	save_gui_params(handles);
    save(file_path, 'Data', 'GUI_params');
    set(hObject,'Enable','off');
    handles.notes.saved	=	true;
    guidata(hObject, handles)
end

end

% --- Executes during object creation, after setting all properties.
function pushbutton_notes_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.notes.folder_name	=	[];
handles.notes.file_name		=	[];
handles.notes.file_path		=	[];

guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_notes_select.
%  Used to select either an annotations or an automated detector file.
function pushbutton_notes_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%

dialog_title	=	'Select location for Annotation files';
start_path		=	handles.notes.folder_name;
folder_name		=	uigetdir(start_path, dialog_title);

if isnumeric(folder_name)
    return;
end

handles		=	load_notes_file(handles, folder_name);  %notes_select callback
guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function checkbox_notes_readonly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_notes_readonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.notes.readonly	=	false;
guidata(hObject, handles);
end

% --- Executes on button press in checkbox_notes_readonly.
function checkbox_notes_readonly_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_notes_readonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_notes_readonly
read_only	=	get(hObject, 'Value');
saved		=	handles.notes.saved;
if	read_only || saved
    opt		=	'off';
else
    opt		=	'on';
end
set(handles.pushbutton_notes_save, 'Enable', opt);

handles.notes.readonly	=	read_only;

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function checkbox_notes_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_notes_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.notes.show	=	false;
handles.notes.i_show=	[];
handles.notes.h_show=	[];
handles.notes.i_sel =	[];
handles.notes.h_info=	[];

guidata(hObject, handles);
end

% --- Executes on button press in checkbox_notes_show.
function checkbox_notes_show_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_notes_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_notes_show

handles.notes.show	=	get(hObject, 'Value');
handles				=	plot_events(handles);
guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_notes_new.
function handles=pushbutton_notes_new_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Click on two points defining a rectangle on the figure: ');
[Times, Freq, Buttons]		=	ginput(2);

%	Check if user tried to terminate input early
if	(length(Buttons) < 2) || any(Buttons > 3)
    disp('Input cancelled');
    return;
end

%	Parameters taken from ginput
start_time	=	handles.tdate_start...
    +	datenum(0,0,0,0,0,min(Times));
min_freq	=	(1000*min(Freq));
max_freq	=	(1000*max(Freq));
duration	=	(abs(Times(2) - Times(1)));

%	Draw square on plot
x		=	min(Times);
y		=	min_freq/1000;
width	=	duration;
height	=	(max_freq - min_freq)/1000;
axes(handles.axes1);
hrec	=	rectangle('Position',[x,y,width,height],...
    'Curvature',[0.3],...
    'LineWidth',2,'LineStyle','-',...
    'EdgeColor','r');

%	Initial signal type selection--can change when editing event.
% sig_types	=	{'Pulsive','FM'};
% choice		=	menu('Signal type?',sig_types);
% if	choice == 0
%     disp('Input cancelled');
%     delete(hrec);
%     return;
% end
% sig_type	=	sig_types{choice};
sig_type = 'FM';

%	Get automated parameters from basic information
params_extract	=	extract_automated_fields(Times, Freq, handles);
if isempty(params_extract)
    uiwait(errordlg({'SNR and level could not be extracted:';...
        'Time window is too small';...
        'Show larger time window and retry'}));
    delete(hrec);
    return
end


%	Get existing notes
Data	=	handles.notes.Data;
i_sel	=	handles.notes.i_sel;

%	Default values from last (selected) event or template
if	~isempty(Data.Events)
    if	~isempty(i_sel)
        Event	=	Data.Events(i_sel);
        
    else
        Event	=	Data.Events(end);
    end
    previous_start_time=Event.start_time;
else
    Event	=	Data.Template;
    previous_start_time=0;
end

%	Insert current values
Event.start_time	=	start_time;
Event.sig_type		=	sig_type;
Event.min_freq		=	min_freq;
Event.max_freq		=	max_freq;
Event.duration		=	duration;



Event.noise_se_dB		=	params_extract.noise_se_dB;
Event.noise_rms_dB		=	params_extract.noise_rms_dB;
Event.noise_peakpsd_dB	=	params_extract.noise_peakpsd_dB;
Event.signal_se_dB		=	params_extract.signal_se_dB;
Event.signal_rms_dB		=	params_extract.signal_rms_dB;
Event.signal_peakpsd_dB	=	params_extract.signal_peakpsd_dB;
Event.SNR_rms_dB		=	params_extract.SNR_rms_dB;



%	modified defualts for other signal types
switch sig_type
    case	'Pulsive'
        Event.call_type		=	'S1';
        Event.num_pulses	=	10;
        Event.num_harmonics	=	-1;
    case	'FM'
        Event.call_type		=	'moan';
        Event.num_pulses	=	-1;
        Event.num_harmonics	=	1;
    otherwise
        warning('Signal type not recognized for defaults');
end

%Finally, determine hash tag...


%if linked file, offer additional capabilities bearing of annotation
linked_file=isfield(handles.file_flags,'linked')&&handles.file_flags.linked&&isfield(Event,'Istation');

%switch handles.filetype
%case 'GSI'
if ~linked_file
    
    %If event is within 2 sec of previously-selected annotation, offer a chance to use
    %  same hashtags (useful for harmonics)
    tmp=(abs(start_time-previous_start_time));
    yes_independent_selection=tmp>datenum(0,0,0,0,0,2);
    if ~yes_independent_selection
        ButtonName = questdlg('You are close to previously-selected annotation!', ...
            'Annotation status', ...
            'Completely New', 'Adding harmonic to previous', 'Completely New');
    else
        ButtonName='Completely New';
    end
    if strcmp(ButtonName,'Completely New')||~isfield(Event,'hash_tag')||isempty(Event.hash_tag)
        
        Event.hash_tag      =   2*datenum(1970,1,1,0,0,0)-now;  %manual hashtags go back into past, automated hashtags into future.
        %If this is a brand new creation with empty link file, zero out
        %hashtags and position
        
    else
        disp('copying hashtag...');
        %Keep current event hashtag...
    end
else
    
    ButtonName = questdlg('What is this annonation status?', ...
        'Annotation status', ...
        'Completely New', 'Adding harmonic', 'Replacing','Completely New');
    
    if strcmp(ButtonName,'Completely New')||~isfield(Event,'hash_tag')||isempty(Event.hash_tag)
        Event.hash_tag      =   2*datenum(1970,1,1,0,0,0)-now;  %manual hashtags go back into past, automated hashtags into future.
        Event.link_hashtags=-1*ones(size(Event.link_hashtags,1),1);
        Istation=str2num(Event.Istation);
        Event.link_hashtags(Istation)=Event.hash_tag;
        Event.link_hashtags=num2str(Event.link_hashtags);
    end
    
    
    if isfield(Event,'bearing')
        
        switch handles.filetype
            case 'GSI'
                [bearing,kappa,tsec]=get_GSI_bearing(hObject,eventdata,handles,[Times Freq]);
                %Event.bearing=num2str(bearing);
                %if isfield(Event,'localization')&&strcmp(ButtonName,'Replacing')
                if isfield(Event,'localization')&&~isempty(getfield(Event,'localization'))
                    Event=update_GSI_localization(Event,str2num(Event.Istation),bearing,kappa);
                    
                end
        end
    end
end


%	JIT: Unnecssary, removed recomputation, it's done above.
Event	=	edit_event(Event, Data.Description, handles.notes.edit_fields);

try %#ok<*TRYNC>
    delete(hrec);
end

if	~isempty(Event)
    %	Add new event to data store
    [Data.Events, ii]	=	add_event(Data.Events, Event);
    handles.notes.Data	=	Data;
    
    handles.notes.saved	=	false;
    checkbox_notes_readonly_Callback(handles.checkbox_notes_readonly, [], handles);
    set(handles.checkbox_notes_show, 'Enable', 'on');
    
    %	Set new note as currently selected one and replot if enabled
    handles.notes.i_sel	=	ii;
    handles		=	plot_events(handles);
    guidata(hObject, handles);
end
end

% --- Executes on button press in pushbutton_notes_edit.
function pushbutton_notes_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if	isempty(handles.notes.i_sel)
    warning('No event selected');
    return;
else
    i_sel	=	handles.notes.i_sel;
end

if	isempty(handles.notes.Data) || isempty(handles.notes.Data.Events)
    warning('No events data');
    return;
end

%	Delete entry if enabled
delete_on	=	get(handles.checkbox_notes_delete, 'Value');
if	delete_on
    
    current_time=handles.notes.Data.Events(i_sel).start_time;
    handles.notes.Data.Events(i_sel)	=	[];
    N	=	length(handles.notes.Data.Events);
    start_times=[handles.notes.Data.Events.start_time];
    if	N == 0
        i_sel	=	[];
    elseif	i_sel > N
        i_sel	=	N;  %Closest by default
    else
        [~,i_sel]=min(abs(current_time-start_times));
    end
    handles.notes.i_sel	=	i_sel;
    
    %Remove delete as a safety...
    set(handles.checkbox_notes_delete, 'Value',0);
    set(handles.pushbutton_notes_edit,'String','Edit');
    %	Otherwise open edit window with existing values
else
    Event	=	handles.notes.Data.Events(i_sel);
    
    %	Recalculate noise and signal stats
    %	JIT: datevec components unnecessary
    tmp			=	Event.start_time - handles.tdate_start;
    Times(1)	=	tmp*24*60*60;
    Times(2)	=	Times(1) + str2double(Event.duration);
    Freq(1)		=	str2double(Event.min_freq)/1000;
    Freq(2)		=	str2double(Event.max_freq)/1000;
    
    params_extract	=	extract_automated_fields(Times, Freq, handles);
    if isempty(params_extract)
        uiwait(errordlg('SNR and level could not be extracted: Time window is too small'));
        NewEvent=[];
        return
    else
        names	=	fieldnames(params_extract);
        for	ii	=	1:length(names)
            if	~isfield(Event,names{ii})
                error('Event does not contain field for this parameter');
            end
            Event.(names{ii})	=	params_extract.(names{ii});
        end
    end
    
    Event	=	edit_event(Event, handles.notes.Data.Description, handles.notes.edit_fields);
    
    %Check if need to relocalize...
    if isfield(Event,'bearing')&&~isempty(getfield(Event,'localization'))
        
        switch handles.filetype
            case 'GSI'
                [bearing,kappa,tsec]=get_GSI_bearing(hObject,eventdata,handles,[Times' Freq']);
                %Event.bearing=num2str(bearing);
                %if isfield(Event,'localization')&&strcmp(ButtonName,'Replacing')
                if isfield(Event,'localization')
                    Event=update_GSI_localization(Event,str2num(Event.Istation),bearing,kappa);
                    
                end
        end
    end
    
    if	isempty(Event)
        return;
    else
        handles.notes.Data.Events(i_sel)	=	Event;
    end
end

handles.notes.saved	=	false;
checkbox_notes_readonly_Callback(handles.checkbox_notes_readonly, [], handles);
%	Replot events, next one should be highlighted if it's on same screen
handles		=	plot_events(handles);
guidata(hObject,handles);

end

% --- Executes on button press in checkbox_notes_delete.
function checkbox_notes_delete_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_notes_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_notes_delete
delete_on	=	get(hObject, 'Value');

if	delete_on
    set(handles.pushbutton_notes_edit, 'String', 'Delete');
else
    set(handles.pushbutton_notes_edit, 'String', 'Edit');
end

end

function edit_folder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_folder as text
%        str2double(get(hObject,'String')) returns contents of edit_folder as a double
end

% --- Executes during object creation, after setting all properties.
function edit_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton_notes_screen.
function i_show=pushbutton_notes_screen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_screen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h_axes	=	handles.axes1;
%	Window limits
Times	=	mean(xlim(h_axes));
Times	=	handles.tdate_start + datenum(0,0,0,0,0,Times);

%	Event times
Start_Times	=	cell2mat({handles.notes.Data.Events.start_time});

%	AARON Find closest event
[junk,i_show]	=	min(abs(Start_Times-Times));

if ~isempty(i_show)
    handles.notes.i_sel=i_show;
end

handles		=	update_events(handles);

%	this might be redundant
guidata(hObject, handles);

end

function pushbutton_previous_linked_annotation_Callback(hObject, eventdata, handles)
shift_linked_annotation(hObject,eventdata,handles,'prevlink');
end
% --- Executes on button press in pushbutton_next_linked_annotation.
function pushbutton_next_linked_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next_linked_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shift_linked_annotation(hObject,eventdata,handles,'nextlink');

end

% --- Executes on button press in pushbutton_next_station.
function pushbutton_next_station_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shift_linked_annotation(hObject,eventdata,handles,'nextstation');
end

% --- Executes on button press in pushbutton_prev_station.
function pushbutton_prev_station_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_prev_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shift_linked_annotation(hObject,eventdata,handles,'prevstation');
end

function [newdir,newfile,new_annotation_file_name,Iarray,Iarray_org]=get_new_filename(current_file_name, ...
    mydir, myfile,filetype,stepp_type,link_names,link_hashtags)

% Returns a new filename, new directory filename, and new annotation file name to load when a link pushbutton is pressed
newdir=[];
newfile=[];
Iarray=[];
Iarray_org=[];
new_annotation_file_name=[];

if ~exist('link_names')
    link_names=[];
end

if ~exist('link_hashtags')
    link_hashtags=[];
end

if ~isempty(strfind(stepp_type,'prev'))
    stepp=-1;
else
    stepp=1;
end

if ~isempty(strfind(stepp_type,'link'))
    yes_nextlink_only=true;
else
    yes_nextlink_only=false;
end

switch filetype
    case 'GSI'
        if isempty(link_names)
            return
        end
        current_letter=current_file_name(5);
        DASAR_letters=link_names(:,5);
        NDASAR=length(DASAR_letters);
        Iarray_org=strmatch(current_letter,DASAR_letters);  %Position of current annotationfile in link_names;
        
        %Check that station number matches Iarray_org (derived from
        %strings)
        
        
        if stepp==1
            Iarray=Iarray_org+1;
            if Iarray>NDASAR
                Iarray=1; %Return to letter 'A'
            end
        else
            Iarray=Iarray_org-1;
            if Iarray==0
                Iarray=NDASAR; %Return to letter 'A'
            end
            
        end
        
        %Keep going until a linked annotation is found
        if yes_nextlink_only&~isempty(link_names)&~isempty(link_hashtags)
            next_link_tag=str2num(link_hashtags(Iarray,:));
            while(next_link_tag<0)
                if stepp==1
                    Iarray=Iarray+1;
                    if Iarray>NDASAR
                        Iarray=1; %Return to letter 'A'
                    end
                else
                    Iarray=Iarray-1;
                    if Iarray==0
                        Iarray=NDASAR; %Return to letter 'A'
                    end
                    
                end
                next_link_tag=str2num(link_hashtags(Iarray,:));
            end
            
            
        end
        next_letter=DASAR_letters(Iarray);
        
        if ~strcmp(mydir(end-6),'S')
            fprintf('Directory %s does not have form S***** \n',mydir);
            return
        end
        mydir(end-2)=next_letter;
        %  Assume S510G0T20100831T000000.gsi form
        myfile(5)=next_letter;
        fprintf('Directory %s contains %s\n',mydir,myfile);
        
        newdir=mydir;
        newfile=myfile;
        
        %get annotation file name
        
        %Check that everything makes sense
        if ~strcmp(current_file_name,link_names(Iarray_org,:))
            fprintf('%s in notes files does not match %s in link file names \n', file_name,link_names(Iarray_org,:));
            
        end
        
        new_annotation_file_name=link_names(Iarray,:);
end
end

function handles=shift_linked_annotation(hObject,eventdata,handles,stepp_type)

annotation_exists=true;
switch handles.filetype
    case 'GSI'
        
        %Determine next DASAR to load, keep skipping until a hashtag is
        %  discovered...
        if ~isfield(handles,'notes')||isempty(handles.notes)
            return
        end
        
        i_sel=handles.notes.i_sel;
        if isempty(i_sel) %We have no annotation
            annotation_exists=false;
            i_sel=1;
            handles.notes.i_sel=1;
            stepp_type=[stepp_type(1:4) 'station'];  %Can't link if no current annotation
        end
        currentEvent=handles.notes.Data.Events(i_sel);
        
        % If no annotation selected, get link_names through other means.
        
        [mydir,myfile,new_annotation_file_name,Iarray,Iarray_org]=get_new_filename(handles.notes.file_name, ...
            handles.mydir, handles.myfile, handles.filetype,stepp_type, ...
            currentEvent.link_names,currentEvent.link_hashtags);
        
        if isempty(mydir)||isempty(myfile)
            return
        end
        
        handles.mydir=mydir;
        handles.myfile=myfile;
        
        if str2num(currentEvent.Istation)~=Iarray_org
            disp('Warning!! Events.Istation not equal to Iarray_org');
        end
        
        try
            handles		=	load_and_display_spectrogram(handles);
        catch
            disp(sprintf('Could not load %s',fullfile(handles.mydir,handles.myfile)));
        end
        
        
        %%Now that event has been copied, delete annotation if desired
        delete_on	=	get(handles.checkbox_notes_delete, 'Value');
        if delete_on & annotation_exists
            %Double check
            ButtonName = questdlg('Delete annotation before changing link?', ...
                'The Delete checkbox is checked!', ...
                'Yes','No', 'No');
            switch ButtonName
                case 'Yes'
                    %Remove existance of annotation from CurrentEvent
                    %hashtables
                    
                    currentEvent=update_GSI_localization(currentEvent,Iarray_org, NaN,-1);
                    
                    %Delete original annotation
                    handles.notes.Data.Events(i_sel)	=	[];
                    N	=	length(handles.notes.Data.Events);
                    %Keep same index, unless at end
                    if	N == 0
                        i_sel	=	[];
                    elseif	i_sel > N
                        i_sel	=	N;
                    end
                    handles.notes.i_sel	=	i_sel;
                    handles.notes.saved=false;
                    %Explicitly turn off delete checkbox to eliminate
                    %inadvertant deletions.
                    set(handles.checkbox_notes_delete,'Value',false)
                    checkbox_notes_delete_Callback(handles.checkbox_notes_delete, eventdata, handles);
                    
            end % switch
        end %delete_on
        
        % Load and replace current notes with new file
        fprintf('I am file %s, station %i just before load_notes_file\n',handles.notes.file_name,str2num(handles.notes.Data.Events(handles.notes.i_sel).Istation));
        
        if annotation_exists
            pushbutton_notes_save_Callback(handles.pushbutton_notes_save, eventdata, handles);  %Force a save
            handles.notes.saved=true;
        end
        
        handles=load_notes_file(handles, [],new_annotation_file_name);
        handles.notes.saved=true;
        
        %Create hashtag array
        hashtag=-1*ones(length(handles.notes.Data.Events),1);
        for I=1:length(handles.notes.Data.Events)
            hashtag(I)=str2num(handles.notes.Data.Events(I).hash_tag);
        end
        
        
        %Determine whether a linked annotation already exists in this
        %station
        
        if annotation_exists
            next_link_tag=str2num(currentEvent.link_hashtags(Iarray,:));
            handles.notes.i_sel=find(next_link_tag==hashtag, 1 );
            
            handles	=	plot_events(handles);
        else % We are at a station with no annotation
            handles	=	plot_events(handles);
            handles.notes.i_sel=min(handles.notes.i_show);
            guidata(hObject, handles);
            return
        end
        
        flag_fix_event=false;
        if ~isempty(findstr(stepp_type,'link'))
            
            if isempty(handles.notes.i_sel) %The data in link_hashtags is out of date!  Annotation has been deleted
                msgbox('Warning!  The annotation does not exist here!');
                return
            end
            
        else  %If I'm switching channel
            
            ButtonName = questdlg('What Next?', ...
                'Linking choice', ...
                'Do Nothing', 'Add link', 'Replace','Do Nothing');
            axes(handles.axes1)
            
            switch ButtonName
                case 'Do Nothing'
                    %If no linked annotation exists in this station,
                    %return
                    if isempty(handles.notes.i_sel)
                        % i_sel=pushbutton_notes_screen_Callback(handles.pushbutton_notes_screen,eventdata,handles);
                        % handles.notes.i_sel=i_sel;
                        guidata(hObject, handles);
                        return
                    end
                    
                    %Otherwise, continue as if next_link button has
                    %been pressed
                case 'Add link'
                    handles=pushbutton_notes_new_Callback(handles.pushbutton_notes_new, eventdata, handles);
                    
                    % handles.notes.i_sel should be updated
                    % automatically
                case 'Replace'
                    uiwait(msgbox('After hitting OK, please click on left side of annotation you wish to link'));
                    tmp=ginput(1);
                    [~,ichc]=min(abs(tmp(1)-handles.notes.x_show));
                    handles.notes.i_sel=handles.notes.i_show(ichc);
                    handles	=	plot_events(handles);
                    
            end
            
            flag_fix_event=strcmp(ButtonName,'Add link')|strcmp(ButtonName,'Replace');
            
            
        end
        
        %Create the new Event
        newEvent=handles.notes.Data.Events(handles.notes.i_sel);  %Event that is on new station
        
        fprintf('I am file %s, station %i\n',handles.notes.file_name,str2num(newEvent.Istation));
        
        %As a courtesy, copy over author information
        if ~strcmp(newEvent.author,currentEvent.author)
            newEvent.author=currentEvent.author;
            newEvent.call_type=currentEvent.call_type;
            
            handles	=	plot_events(handles);
        end
        %Check that localization information has not changed.  If it has,
        %it means the position has been recalculated, or that original
        % calculation had an incorrect error ellipse.
        
        flag_positions_changed=compare_localizations(currentEvent,newEvent);
        
        
        %If bearings have changed, copy localization object to newEvent and
        %update related fields
        handles.notes.saved=true;
        
        if flag_positions_changed
            
            %Copy over link_names and link_hashtags
            newEvent.link_names=currentEvent.link_names;
            newEvent.link_hashtags=currentEvent.link_hashtags;
            
            %Add current hashtag to link_hashtags
            %Warning, need to make sure not to cut off hashtag...
            if isnumeric(newEvent.hash_tag)
                newEvent.hash_tag=num2str(newEvent.hash_tag);
            end
            
            temp=char(newEvent.link_hashtags,newEvent.hash_tag);  %Ensures blanks are correct
            temp(Iarray,:)=temp(end,:);
            newEvent.link_hashtags=temp(1:(end-1),:);
            %Nchar=min([length(newEvent.hash_tag) size(newEvent.link_hashtags,2)]);
            % newEvent.link_hashtags(Iarray,:)=blanks(size(newEvent.link_hashtags,2));
            % newEvent.link_hashtags(Iarray,1:Nchar)=newEvent.hash_tag(1:Nchar);
            
            %copy localization fields, this ensures that new bearings
            %  are propagated throught linked files
            
            if flag_fix_event
                bearing=newEvent.localization.bearings_all(Iarray);
                kappa=newEvent.localization.kappa(Iarray);
                
                newEvent.localization=currentEvent.localization;
                
                newEvent.localization.bearings_all(Iarray)=bearing;
                newEvent.localization.kappa(Iarray)=kappa;
                
            else
                newEvent.localization=currentEvent.localization;
                
                
            end
            %Update range, bearing, and position, and error ellipse fields
            
            newEvent=update_GSI_localization(newEvent);
            
            %,Iarray_org, str2num(currentEvent.bearing),currentEvent.localization.kappa(Iarray_org));
            
            handles.notes.Data.Events(handles.notes.i_sel)=newEvent;
            handles.notes.saved=false;  %position now changed
            handles=plot_events(handles);
        end
        %handles.notes.i_show=handles.notes.i_sel;
        
        
        guidata(hObject, handles);
    otherwise
        return
end

end

function flag_positions_changed=compare_localizations(currentEvent,newEvent)
%check if bearings consistent
bearings1=currentEvent.localization.bearings_all;
bearings2=newEvent.localization.bearings_all;

Inan1=find(isnan(bearings1));
Inan2=find(isnan(bearings2));
if length(Inan1)~=length(Inan2)
    test_change1=true;
    test_change2=true;
else
    test_change1=~all(Inan1==Inan2);
    Inan1=find(~isnan(bearings1));
    Inan2=find(~isnan(bearings2));
    test_change2=~all(bearings1(Inan1)==bearings2(Inan2));
    
end

flag_positions_changed=test_change2||test_change1;

%%Next, check the error ellipse of the new Event
DASAR_coords=[newEvent.localization.station_position.easting newEvent.localization.station_position.northing];
Igood=find(~isnan(bearings2));
[VM,Qhat,~,outcome] = vmmle_r(bearings2(Igood),DASAR_coords(Igood,:),'h');
mean_coords=mean(DASAR_coords(Igood,:));
CRITVAL=4.60517; %chi2inv(0.90,2);
[~,A,B,ANG,Baxis] = ellipsparms(Qhat,CRITVAL,mean_coords,VM);

%If localization info different than what is stored, flag for a
%change

Astored=newEvent.localization.major;
Bstored=newEvent.localization.minor;
Estored=newEvent.localization.ellipse_ang;

if ~strcmp(outcome,'successful')||isnan(A)||isnan(B)
    flag_positions_changed=true;
    return
end

if abs((A-Astored)./Astored)>0.1||abs((B-Bstored)./Bstored)>0.1
    flag_positions_changed=true;
    return
end
% if isnan(A)&~isnan(

end

%	Supporting functions, i.e. not auto-generated callbacks
function	status	=	dependency_check()

status	=	true;

end

%% load,edit, display, and save annotation files to Ulysses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   ANNOTATIONS annotations %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Checks Notes folder for existing files and loads them
%   Where new annotations are created...

%	disables (or enables) note controls
function	disable_notes_nav(handles, opt)

if	nargin < 2 || isempty(opt)
    opt		=	'off';
end

set(handles.pushbutton_notes_first, 'Enable', opt);
set(handles.pushbutton_notes_last, 'Enable', opt);
set(handles.pushbutton_notes_next, 'Enable', opt);
set(handles.pushbutton_notes_prev, 'Enable', opt);
set(handles.pushbutton_notes_screen, 'Enable', opt);
if strcmpi(handles.filetype,'gsi')
    set(handles.pushbutton_next_linked_annotation, 'Enable', opt);
    set(handles.pushbutton_previous_linked_annotation, 'Enable', opt);
end

end  %disable_notes_nav


function	handles		=	load_notes_file(handles, new_folder,annotation_file_name)
%new_folder can be empty if want to use annotation_file_name

%%	First check that current notes are saved before proceeding
if	~handles.notes.saved && ~handles.notes.readonly
    qtitle	=	'Unsaved changes';
    qstring	=	{'Warning: current changes to notes not saved!'...
        'Write to file before continuing?'};
    str1	=	'Yes';
    str2	=	'No';
    button	=	questdlg(qstring, qtitle, str1, str2, str1);
    
    switch	button
        case	str1
            %	Save file
            pushbutton_notes_save_Callback(handles.pushbutton_notes_save,...
                [], handles);
            handles		=	guidata(handles.pushbutton_notes_save);
        case	str2
            %	Don't bother
        otherwise
            error('This should not be reached');
    end
end


%%	New notes file disables all notes controls, until otherwise activated
opt		=	'off';
disable_notes_nav(handles, opt);
set(handles.pushbutton_notes_new, 'Enable', opt);
set(handles.pushbutton_notes_save, 'Enable', opt);
set(handles.pushbutton_notes_edit, 'Enable', opt);
set(handles.checkbox_notes_show, 'Enable', opt);
set(handles.checkbox_notes_delete, 'Enable', opt);

set(handles.pushbutton_next_linked_annotation, 'Enable', opt);
set(handles.pushbutton_previous_linked_annotation, 'Enable', opt);


%%	Switch to new folder and check for existing files
if	exist('new_folder','var') && ~isempty(new_folder)
    folder_name	=	new_folder;
else
    folder_name	=	handles.notes.folder_name;
end
[~,fname,extt]		=	fileparts(handles.myfile);

%	Default output file name
user_name	=	getusername();
multiple_times_present=0;


if strfind(handles.myfile,'*')  %If inputs indicates that multiple files from different times are to be merged
    Istart=3;
    file_name	=	['Multiple' fname(Istart:end) '-notes-' user_name '.mat'];
    multiple_times_present=1;
    
else
    Istart=1;
    file_name	=	[fname(Istart:end) '-notes-' user_name '.mat'];
    
end


%	Jit: this should be done after finding all relevant files
%	no point showing user this menu if there are no files to select
%% AARON:  Give user option of loading automated detector
%  or existing notes file
file_opts	=	{'Manual Annotation', 'Automated Detector'};
file_exts	=	{'.mat', '.detsum'};
file_flag	=	{true, false};

Nf			=	length(file_exts);
file_found	=	true(Nf,1);

for ii	=	1:Nf  %For every file type...
    %	Find all existing files
    file_listings{ii}	=	dir(fullfile(folder_name, [fname '*' file_exts{ii}]));
    if isempty(file_listings{ii})
        file_found(ii)	=	false;
    end
end
%	Limit options to available files
file_opts		=	file_opts(file_found);
file_exts		=	file_exts(file_found);
file_flag		=	file_flag(file_found);
file_listings	=	file_listings(file_found);

choice	=	0;
if	~isempty(file_listings)
    if length(file_listings) == 1
        choice	=	1;
    else
        %	Ask user to choose which list
        choice	=	menu('Select a file type:', file_opts);
    end
end

if	choice	==	0
    listing		=	[];
    manual_flag	=	true;
else
    listing		=	file_listings{choice};
    manual_flag	=	file_flag{choice};
end

%Check to see if a specific annotation file has already been requested...

if nargin>2 && ~isempty(annotation_file_name)  %if annotation_file_name has been fed in...
    %cd(handles.outputdir)  %Should already be here, but just in case
    mydirr=pwd;
    cd(folder_name)
    if exist(annotation_file_name,'file')==2
        file_name	=	annotation_file_name;
        sel_names{1}=file_name;
        Sel=1;
    else
        fprintf('Requested annotation file %s does not exist!\n',annotation_file_name);
        sel_names=[];
    end
    cd(mydirr);
else  %have user select annotation file interactively
    sel_names	=	[];
    if	~isempty(listing)
        %	Ask user which file(s) to load
        N_files		=	length(listing);
        list_files	=	cell(N_files,1);
        for	ii	=	1:N_files
            list_names{ii}	=	listing(ii).name;
        end
        
        %AARON note: example of listdlg
        [Sel, OK]	=	listdlg('ListString', list_names,...
            'Name', 'Available notes',...
            'PromptString', 'Select file(s) to load:',...
            'OKString', 'Load',...
            'CancelString', 'New');
        if	(OK == 0) || isempty(Sel)
            sel_names	=	[];
        else
            sel_names	=	list_names(Sel);
            
            %	If just one file, then make it the target of future saves
            if	length(Sel) == 1 && manual_flag
                
                file_name	=	sel_names{1};
            elseif length(Sel) == 1 && ~manual_flag
                %%keep file name
                
                
            elseif	length(Sel) > 1
                %	Otherwise new merged file
                user_name	=	'merged';
                file_name	=	[fname(Istart:end) '-' user_name '.mat'];
            end
        end
        
        %	Make sure file name for new files is unique
        if	length(Sel) ~= 1	%i.e. we're not working with a specific file
            ii	=	0;
            while	any(strcmp(file_name, list_names))
                ii	=	ii + 1;
                file_name	=	[fname(Istart:end) '-' user_name '-' num2str(ii) '.mat'];
            end
        end
        
    end  %if listing is not empty (i.e. an annotation file already exists)
end  %annotation_file_name

%%	New file with default data template
[Defaults.Description, Defaults.Template, edit_fields]	=	load_default_annotation_template();
Defaults.Events		=	Defaults.Template;
GUI_params	=	[];
handles.notes.edit_fields	=	edit_fields;


if isempty(listing) || isempty(sel_names)
    %	Create defaults if none already present
    %	Default prompt
    
    Data				=	Defaults;
    Data.Events			=	[];
    
    handles.notes.Data	=	Data;
    handles.notes.show	=	false;
    opt					=	'off';
    
    %	Load selected files and merge data
else
    Data	=	[];
    for	ii	=	1:length(sel_names)  %For every file selected
        file_path		=	fullfile(folder_name, sel_names{ii});
        
        %%AARON changes
        if manual_flag
            LSfile		=	load(file_path);
        else  %import automated file
            %import_JAVA_Energy_into_notes(file_path,sel_names,Defaults,Template,get_dialog)
            LSfile.Data=import_JAVA_Energy_into_notes(file_path,sel_names{ii},Defaults,ii==1);
            %modify saved file name to include time range
            if multiple_times_present
                t1=datestr(LSfile.Data.Events(1).start_time,30);
                t2=datestr(LSfile.Data.Events(end).start_time,30);
                file_name	=	['Multiple_'  t1 '_' t2  fname(2:end) '-notes-' user_name '.mat'];
                
            end
        end
        
        %	Check and merge event data
        if	isempty(Data)
            Data	=	check_notes(LSfile.Data);
        else
            nData	=	check_notes(LSfile.Data);
            Data	=	merge_data(Data, nData);
        end
        
        %	Load gui parameters if present, and only from first encountered
        if	isfield(LSfile, 'GUI_params')
            if isempty(GUI_params)
                GUI_params	=	LSfile.GUI_params;
            end
        end
        
        %AARON:  if annotation contants own edit field matrix, use it
        if	isfield(LSfile.Data, 'edit_fields')
            if ~isempty(LSfile.Data.edit_fields)
                handles.notes.edit_fields	=	LSfile.Data.edit_fields;
            end
        end
        
    end  %ii loop through all names.
    
    %	Always merge with current template data, to update old notes with
    %	new fields
    Data	=	merge_data(Defaults, Data);
    Data.Events(1)	=	[];
    
    
    if	length(Sel) > 1 && manual_flag
        %	Since this is a new merged file, turn save on, and read_only off
        handles.notes.saved	=	false;
        set(handles.checkbox_notes_readonly, 'Value', 0);
        checkbox_notes_readonly_Callback(handles.checkbox_notes_readonly, [], handles);
        set(handles.pushbutton_notes_edit,'Enable','on');  %edit!
        
    elseif ~manual_flag %If automated data loaded, freeze ability to edit but allow saves
        handles.notes.saved	=	false;
        set(handles.checkbox_notes_readonly, 'Value', 0);
        checkbox_notes_readonly_Callback(handles.checkbox_notes_readonly, [], handles);
        set(handles.pushbutton_notes_edit,'Enable','off');  %No editing!
        
    end
    
    handles.notes.Data	=	Data;
    handles.notes.show	=	true;
    opt				=	'on';
end



%%	enable relevant buttons
disable_notes_nav(handles,opt);
set(handles.checkbox_notes_show, 'Value', handles.notes.show);
set(handles.checkbox_notes_show, 'Enable', opt);
if	handles.fig_updated
    set(handles.pushbutton_notes_new, 'Enable', 'on');
end

handles.notes.folder_name	=	folder_name;
handles.notes.file_name		=	file_name;
handles.notes.file_path		=	fullfile(folder_name, file_name);


h_axes	=	handles.axes1;
%	Window limits
Times	=	mean(xlim(h_axes));
handles.notes.i_sel		=	[];

if isfield(handles,'tdate_start')&~isempty(handles.tdate_start)
    
    Times	=	handles.tdate_start + datenum(0,0,0,0,0,Times);
    
    %	Jit, only do this if existing Event data is loaded
    if	~isempty(Data.Events)
        %	Event times
        Start_Times	=	cell2mat({Data.Events.start_time});
        
        %	AARON Find closest event
        [~, i_show]	=	min(abs(Start_Times - Times));
        
        if ~isempty(i_show)
            handles.notes.i_sel	=	i_show;
        end
    end
end
%	Load previous window settings, if present
if ~exist('annotation_file_name','var')  %don't load if I am skipping across files
    handles		=	load_gui_params(handles, GUI_params);
end
%	Set folder text box to selected directory + file_name
set(handles.edit_folder, 'String', handles.notes.file_path);

%	Set note index text accordingly
i_sel		=	handles.notes.i_sel;
N			=	length(handles.notes.Data.Events);
msg			=	{[num2str(i_sel) ' of']; num2str(N)};
set(handles.text_notenum, 'String', msg);

end  %load_notes_file
% --------------------------------------------------------------------

%	Load and apply ancilary window information from notes if present
function	handles		=	load_gui_params(handles, GUI_params)

if	isempty(GUI_params)
    return;
end

%	FFT_parameters
set(handles.popupmenu_Nfft,'Value', GUI_params.iNfft);
set(handles.popupmenu_ovlap,'Value', GUI_params.iOvlap);
set(handles.edit_mindB,'String', num2str(GUI_params.mindB));
set(handles.edit_dBspread,'String', num2str(GUI_params.dBspread));
handles.tlen	=	GUI_params.tlen;
set(handles.edit_winlen,'String', num2str(GUI_params.tlen));

%	Playback parameters
handles.filter.f_min	=	GUI_params.filter.f_min;
handles.filter.f_max	=	GUI_params.filter.f_max;
handles.filter.changed	=	true;
set(handles.edit_minfreq,'String', num2str(GUI_params.filter.f_min));
set(handles.edit_maxfreq,'String', num2str(GUI_params.filter.f_max));

%	Plot parameters
set(handles.edit_fmin,'String', num2str(GUI_params.fmin));
set(handles.edit_fmax,'String', num2str(GUI_params.fmax));
set(handles.edit_datestr, 'String', datestr(GUI_params.tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

end  %load_gui_params
% --------------------------------------------------------------------

%	Saves GUI parameters
function	GUI_params		=	save_gui_params(handles)

if	isempty(handles)
    GUI_params	=	[];
    return;
end

%	FFT_parameters
GUI_params.iNfft	=	get(handles.popupmenu_Nfft,'Value');
GUI_params.iOvlap	=	get(handles.popupmenu_ovlap,'Value');
GUI_params.mindB	=	str2double(get(handles.edit_mindB,'String'));
GUI_params.dBspread	=	str2double(get(handles.edit_dBspread,'String'));
GUI_params.tlen		=	handles.tlen;

%	Playback parameters
GUI_params.filter.f_min	=	handles.filter.f_min;
GUI_params.filter.f_max	=	handles.filter.f_max;

%	Plot parameters
GUI_params.fmin			=	str2double(get(handles.edit_fmin,'String'));
GUI_params.fmax			=	str2double(get(handles.edit_fmax,'String'));
GUI_params.tdate_start	=	handles.tdate_start;

end  %save_guid_params
% --------------------------------------------------------------------

%	Pops up window to edit event data
function	NewEvent	=	edit_event(Event, Description, edit_fields)

if	nargin < 3
    edit_fields	=	[];
end

%	Start time included in window title, should not be part of input
OldEvent		=	Event;
Event			=	rmfield(Event, 'start_time');
Description(1)	=	[];

%	only keep those fields we want to edit
if	~isempty(edit_fields)
    keep	=	[];
    names	=	fieldnames(Event);
    for	ii	=	1:length(names)
        if	any(strcmp(names{ii},edit_fields))
            keep	=	[keep; ii];
        else
            Event	=	rmfield(Event, names{ii});
        end
    end
    Description		=	Description(keep);
end

%	Generate default strings from field values
names		=	fieldnames(Event);
N_fields	=	length(names);
defaults	=	cell(N_fields,1);

%	Size of each prompt field
num_lines		=	ones(N_fields,1);
for	ii	=	1:N_fields
    value	=	Event.(names{ii});
    if iscell(value)
        defaults{ii}=cell2mat(value');
        
    else
        defaults{ii}	=	num2str(value);
        
    end
    num_lines(ii)=max([1 size(defaults{ii},1)]);
end
num_lines(end)	=	2; %for comments




%	title for window
dlgTitle	=	['Annotation for event at ' datestr(OldEvent.start_time)];

if	length(Description) ~= N_fields
    error('# of Descriptors differs from # of fields');
end

%Convert structures into inputsdlg

%	Create input dialogbox
%options.Resize	=	'on';
%answer		=	inputdlg(Description, dlgTitle, num_lines, defaults, options);
%answer		=	newid(Description, dlgTitle, num_lines, defaults, options);

prep_inputsdlg;
[answer,Cancelled] = inputsdlg(Prompt,dlgTitle,Formats,DefAns,Options);
%

if	isempty(answer)||Cancelled
    NewEvent	=	[];
    return;
end

NewEvent	=	OldEvent;
for	ii	=	1:N_fields
    NewEvent.(names{ii})	=	answer.(names{ii});
end


end %edit_event
% --------------------------------------------------------------------

function	params_extract	=	extract_automated_fields(Times,Freq,handles)

duration	=	(abs(Times(2) - Times(1)));


T		=	handles.sgram.T;
F		=	handles.sgram.F/1e3;
dF		=	F(2) - F(1);
dT		=	T(2) - T(1);
B		=	handles.sgram.B;
i_time	=	(min(Times) <= T) & (T <= max(Times));
i_freq	=	(min(Freq) <= F) & (F <= max(Freq));
PSD		=	(B(i_freq, i_time));

%Estimating the signal-to-noise ratio in terms of rms power/rms power
% AARON THODE flag, changed June 16, 2013
% Compute noise over window just before and after detection,
%  with total duration equal  to detection length

% Let X(f) be raw FFT
% Then PSD(f) =|X(f)|^2/(Nfft*Fs)
% From parsevals theorem and dt*X(f)=X(w)_analytic
%    rms power over FFT sample~sum(df*PSD) in a single column
%       Thus mean rms is mean(sum(df*PSD));
%       Note that this handles overlap between samples
%    SE across FFT sample:  sum(PSD) per column
%    SE_total for entire signal across all spectrogram columns:
%       Simplest solution is mutliply rms power by duration
%       se=rms*duration;

signal_rms		=	mean(sum(dF*PSD)); %note not technically root mean square, actually mean square
signal_se		=	signal_rms*duration;
signal_peakpsd	=	max(PSD(:));

signal_se_dB		=	10*log10(signal_se);
signal_rms_dB		=	10*log10(signal_rms);  %not 20 log because actually mean sqaure
signal_peakpsd_dB	=	10*log10(signal_peakpsd);

%Take noise sample from just before and after current detection
Tmin	=	min(Times);
Tmax	=	max(Times);
duration_noise	=	0.5*(Tmax-Tmin);
i_time_noise1	=	(Tmin-duration_noise <= T) & (T < Tmin);
i_time_noise2	=	(Tmax < T) & (T <= Tmax+duration_noise);%Shift one column earlier to avoid picking up first signal FFT bin
%if (sum(i_time_noise1)*dT < 0.9*duration_noise) ||...
%        (sum(i_time_noise2)*dT < 0.9*duration_noise)

if Tmin-duration_noise<=min(T) || Tmax+duration_noise>=max(T)
    uiwait(msgbox('extract_automated_fields: Selected annotation too close to start or end of spectrogram window; automated SNR cannot be computed','modal'));
    params_extract	=	[];
    return
end

if sum(i_time_noise1)<2||sum(i_time_noise2)<2
    uiwait(msgbox('extract_automated_fields: Selected annotation too short for chosen FFT size; automated SNR cannot be computed','modal'));
    params_extract	=	[];
    return
end

PSD_noise   =  [B(i_freq,i_time_noise1) B(i_freq,i_time_noise2)];

%noise duration should equal signal duration from above procedure.
noise_rms		=	mean(sum(dF*PSD_noise)); %note not technically root mean square, actually mean square
noise_se		=	noise_rms*duration;
noise_peakpsd	=	max(PSD_noise(:));


noise_se_dB			=	10*log10(noise_se);
noise_rms_dB		=	10*log10(noise_rms); %not 20 log because actually mean sqaure
noise_peakpsd_dB	=	10*log10(noise_peakpsd);

SNR_rms		=	signal_rms/noise_rms;
SNR_rms_dB	=	10*log10(SNR_rms);  % not 20log because ratio of mean-square, not root-mean-square

params_extract.noise_se_dB		=	noise_se_dB;
params_extract.noise_rms_dB		=	noise_rms_dB; %not 20 log because actually mean sqaure
params_extract.noise_peakpsd_dB	=	noise_peakpsd_dB;

params_extract.signal_se_dB		=	signal_se_dB;
params_extract.signal_rms_dB	=	signal_rms_dB;  %not 20 log because actually mean sqaure
params_extract.signal_peakpsd_dB=	signal_peakpsd_dB;

%params_extract.SNR_rms=SNR_rms;
params_extract.SNR_rms_dB		=	SNR_rms_dB;


end  %extract_automated_fields
% --------------------------------------------------------------------

%%load JAVA energy detector results into annotation format
function Data=import_JAVA_Energy_into_notes(file_path,sel_names,Defaults,get_dialog)


persistent filter_transition_band filter_passband_frequencies run_options

[auto,head]		=	readEnergySummary(file_path, Inf);

for Ifea=1:length(auto.names)
    if strcmp(auto.names{Ifea},'max_freq')
        Imaxx=Ifea;
    elseif strcmp(auto.names{Ifea},'min_freq')
        Iminn=Ifea;
    end
end


Data			=	Defaults;
Data.param	=	head;
hh	=	waitbar(0,sprintf('Importing %s...',sel_names));
for JJ = 1:length(auto.ctime)
    if rem(JJ,500) == 0
        waitbar(JJ/length(auto.ctime),hh);
    end
    Data.Events(JJ)			=	Defaults.Template;
    Data.Events(JJ).start_time=	datenum(1970,1,1,0,0,auto.ctime(JJ));
    Data.Events(JJ).author	=	'JAVA Energy Processor';
    Data.Events(JJ).duration	=	num2str(auto.features(end,JJ));
    Data.Events(JJ).min_freq	=	num2str(auto.features(Iminn,JJ));
    Data.Events(JJ).max_freq	=	num2str(auto.features(Imaxx,JJ));
end
close(hh);

%%If snips file exists, load level information.
[pathstr,fname,extt] = fileparts(file_path);
full_snips_name=fullfile(pathstr,[fname '.snips']);
snips_name=dir(full_snips_name);
if isempty(snips_name)
    return
    
end

min_freq_events=min(auto.features(Iminn,:));
max_freq_events=max(auto.features(Imaxx,:));
freq_spacing=min(diff(head.flow));
%Making sure that min FIR passband frequency is below event
%   detected.
min_freq_events=max([10 min_freq_events-freq_spacing]);
max_freq_events=min([head.Fs/2 max_freq_events+freq_spacing+1]);

if get_dialog || isempty(filter_transition_band)  %If this is the first file loaded in a sequence
    prompt={'Number of snips to import into RAM:', ...
        'filter transition band (Hz):',...
        'minimum passband frequency for FIR filters (Hz):', ...
        'maximum passband frequency for FIR filters (Hz):', ...
        'passband frequency spacing (Hz)', ...
        'Debug plotting (1=yes; 0=no)'};
    name='Parameters for importing JAVA snips file';
    numlines=1;
    defaultanswer={'500','2',num2str(min_freq_events),num2str(max_freq_events),num2str(freq_spacing),'0'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    
    run_options.Ncalls_to_sample=str2double(answer{1});
    filter_transition_band=str2double(answer{2});
    min_freq_events=str2double(answer{3});
    max_freq_events=str2double(answer{4})-filter_transition_band;
    freq_spacing=str2double(answer{5});
    run_options.debug=str2double(answer{6});
    filter_passband_frequencies=min_freq_events:freq_spacing:max_freq_events;
    filter_passband_frequencies=unique([filter_passband_frequencies max_freq_events]);
    
end


transient_params=extract_transient_levels(full_snips_name,1:length(auto.ctime),auto, ...
    round(head.Fs),head.bufferTime, ...
    filter_transition_band,filter_passband_frequencies,run_options);


hh	=	waitbar(0,sprintf('Importing %s...',sel_names));
for JJ = 1:length(auto.ctime)
    if rem(JJ,500) == 0
        waitbar(JJ/length(auto.ctime),hh);
    end
    Data.Events(JJ)			=	Defaults.Template;
    Data.Events(JJ).start_time=	datenum(1970,1,1,0,0,auto.ctime(JJ));
    Data.Events(JJ).author	=	'JAVA Energy Processor+snips conversion';
    Data.Events(JJ).duration	=	num2str(transient_params.level.t_Malme(JJ));
    Data.Events(JJ).min_freq	=	num2str(auto.features(Iminn,JJ));
    Data.Events(JJ).max_freq	=	num2str(auto.features(Imaxx,JJ));
    
    if (Data.Events(JJ).duration>0)
        Data.Events(JJ).noise_se_dB		=	num2str(10*log10(transient_params.noise.SE(JJ)));
        Data.Events(JJ).noise_rms_dB		=	num2str(20*log10(transient_params.noise.rms(JJ)));
        Data.Events(JJ).noise_peakpsd_dB	=	0;
        
        Data.Events(JJ).signal_se_dB		=	num2str(10*log10(transient_params.level.SE_Malme(JJ)));
        Data.Events(JJ).signal_rms_dB		=	num2str(20*log10(transient_params.level.rms_Malme(JJ)));
        Data.Events(JJ).signal_peakpsd_dB	=	0;
        
        Data.Events(JJ).SNR_rms			=	num2str(transient_params.level.SNR(JJ));
        Data.Events(JJ).SNR_rms_dB= num2str(10*log10(transient_params.level.SNR(JJ)));
    end
end
close(hh);


end %import_JAVA_Energy_into_notes

%	Adds new event to notes structure and re-sorts chronologically
function	[Event_set, ii]		=	add_event(Event_set, Event)

%	If new set, just return the single event
if	isempty(Event_set)
    Event_set	=	Event;
    ii			=	1;
else
    %	Append new event(s)
    N			=	length(Event_set) + 1;
    Event_set	=	[Event_set(:); Event(:)];
    Start_times	=	cell2mat({Event_set.start_time});
    %	Check if already sorted, i.e. sequential additions
    if	issorted(Start_times)
        ii	=	N;	%	index of 1st inserted item
    else
        [~,I]		=	sort(Start_times);
        Event_set	=	Event_set(I);
        %	Find new index of 1st inserted event
        ii			=	find(I == N);
    end
end

end  %add_event
% --------------------------------------------------------------------

%	Checks existing notes data structure for format compliance
function	Data	=	check_notes(Data)

if	isempty(Data)
    return;
end

req_fields	=	{'Description', 'Template', 'Events'};
if	~isstruct(Data) || ~all(isfield(Data, req_fields))
    errmsg	=	'Existing Data must be a structure containing the fields: \n ';
    for		ii	=	1:length(req_fields)
        errmsg	=	[errmsg req_fields{ii} ', ']; %#ok<AGROW>
    end
    errmsg(end)		=	[];		errmsg(end)		=	[];
    error(errmsg, []);
end

if	~iscell(Data.Description)
    error('Description must be a cell array of strings');
end

if	~isstruct(Data.Template)
    error('Template must be a struct with the desired data field names');
end

N		=	length(Data.Description);
names	=	fieldnames(Data.Template);
if	N ~= length(names)
    error('# of fields in Template do not match # of descriptions');
end

if	~strcmp(names{1}, 'start_time')
    error('First field in Template (and Events) must be start_time');
end

if	isempty(Data.Events)
    return;
else
    names2	=	fieldnames(Data.Events);
end

if	~strcmp(names, names2)
    error('Template and Events have differring structures');
end

%	At this point assume most everything is ok
%	return

%	JIT: this is quick, even with several hundred records
%	no need for waitbar
%	Code to convert all fields (except start_time to strings)
% hh = waitbar(0,sprintf('Checking Data'));

for ii	=	2:N
    
    %	!?!? Existing strings are automatically passed through, but do I
    %	need to specifically check for data types other than numeric?
    Data.Template.(names{ii})		=	num2str(Data.Template.(names{ii}));
    for	jj	=	1:length(Data.Events)
        %         if rem(jj,200)==0
        %             waitbar((ii*length(Data.Events)+jj)/(N*length(Data.Events)),hh);
        %         end
        
        %AARON: it is useful to have fields that have substructures, so I
        %will allow annotations to have cell matricies and subfields
        try
            Data.Events(jj).(names{ii})		=	num2str(Data.Events(jj).(names{ii}));
        end
    end
end
%close(hh)


end %check_notes
% --------------------------------------------------------------------

%	Merge Events, inserting missing fields as needed
function	[Data]		=	merge_data(Data1, Data2)

%	Combine field names, to account for duplicates and missing entries
names1	=	fieldnames(Data1.Events);
names2	=	fieldnames(Data2.Events);
[names, i1, i2]	=	union(names1, names2, 'stable');

%	Combine descriptions in the same order
Data.Description	=	[Data1.Description(i1), Data2.Description(i2)];

%	Use dummy entry to add new field names to each record set
Data1.Events(end+1)	=	Data1.Events(end);
Data2.Events(end+1)	=	Data2.Events(end);
for	ii	=	1:length(names)
    Data1.Events(end).(names{ii})	=	[];
    Data2.Events(end).(names{ii})	=	[];
    
    Data.Template.(names{ii})		=	[];
end
%	Delete dummy entries
Data1.Events(end)	=	[];
Data2.Events(end)	=	[];

%	Populate template with default values from each set
for ii	=	1:length(i1)
    Data.Template.(names1{i1(ii)})	=	Data1.Template.(names1{i1(ii)});
end
for ii	=	1:length(i2)
    Data.Template.(names2{i2(ii)})	=	Data2.Template.(names2{i2(ii)});
end

%	Merge records, now that they are field consistent
Data.Events	=	add_event(Data1.Events, Data2.Events);


%	Manually sort field names so that start_time is always first, and
%	comments last
perm	=	1:length(names);
is		=	find(strcmp('start_time',names));
ic		=	find(strcmp('comments',names));
perm([is ic])	=	[];
perm	=	[is, perm(:).', ic];

names	=	names(perm);
Data.Description	=	Data.Description(perm);
Data.Template		=	orderfields(Data.Template, names);
Data.Events			=	orderfields(Data.Events, names);


end  %merge_data
% --------------------------------------------------------------------

%	Remove obsolete fields from old events, but save data into comments
function	Data		=	clean_data(Data, obsolete_fields)
if	nargin < 2
    obsolete_fields		=	{'noise_db', 'peak_db'};	% default list
end
if	ischar(obsolete_fields)
    temp				=	{obsolete_fields};
    obsolete_fields		=	temp;
    clear temp;
end

%	Find which obsolete fields are present
names	=	fieldnames(Data.Events);
for	ii	=	1:length(obsolete_fields)
    test	=	strcmp(obsolete_fields{ii}, names);
    if	any(test)
        i_of(ii)	=	find(test);
    else
        i_of(ii)	=	0;
    end
end
%	reduce set to those present, and needing removal
obsolete_fields(~logical(i_of))	=	[];
i_of(~logical(i_of))			=	[];

%	loop through events, and append info to comments
for	ie	=	1:length(Data.Events)
    event		=	Data.Events(ie);
    comments	=	event.comments;
    for ii	=	1:length(obsolete_fields)
        name		=	obsolete_fields{ii};
        new_txt		=	num2str(event.(name));
        comments	=	sprintf([comments '\n' name '=' new_txt]);
    end
    event.comments	=	comments;
    Data.Events(ie)	=	event;
end

%	remove obsolete fields
Data.Template			=	rmfield(Data.Template, obsolete_fields);
Data.Events				=	rmfield(Data.Events, obsolete_fields);
Data.Description(i_of)	=	[];


end %clean_data
% --------------------------------------------------------------------

%	Overlays events within current window
function	handles	=	plot_events(handles)


%	Set note text accordingly
i_sel		=	handles.notes.i_sel;
N			=	length(handles.notes.Data.Events);
msg			=	{[num2str(i_sel) ' of']; num2str(N)};
set(handles.text_notenum, 'String', msg);

%%If no notes are opened, leave gracefully
if isempty(handles.notes.file_name)
    return
end

%	Try deleting existing events from window
try
    delete(handles.notes.h_show);
end


%	Disable edit, next/prev buttons if nothing is shown
if	isempty(handles.notes.Data) || isempty(handles.notes.Data.Events)...
        || ~handles.notes.show
    set(handles.pushbutton_notes_edit,'Enable','off');
    set(handles.checkbox_notes_delete,'Enable','off');
    %	Disable prev/next buttons
    disable_notes_nav(handles);
    return;
end
Events	=	handles.notes.Data.Events;
Description	=	handles.notes.Data.Description;

h_axes	=	handles.axes1;
%	Window limits
Times	=	xlim(h_axes);
Times	=	handles.tdate_start + datenum(0,0,0,0,0,Times);

%	Event times
Start_Times	=	cell2mat({Events.start_time});

%	Find events that lie within window
i_show	=	find((Times(1) <= Start_Times) & (Start_Times <= Times(2)));

%	Plot rectangles for visible events
h_show	=	gobjects(1,length(i_show));
x_show	=	zeros(1,length(i_show));

%x_show =h_show;
sel_vis	=	false;
%axes(h_axes);
for	ii	=	1:length(i_show)
    axes(h_axes); %axes must be called inside loop for all annotations to be visible.
    ie		=	i_show(ii);
    event	=	Events(ie);
    x		=	event.start_time - Times(1);	x	=	x*24*60*60;
    x_show(ii)=x;
    y		=	str2double(num2str(event.min_freq))/1000;
    width	=	str2double(num2str(event.duration));
    height	=	(str2double(num2str(event.max_freq))...
        - str2double(num2str(event.min_freq)))/1000;
    
    %AARON fix for dodgy detections
    height	=	max([1e-5 height]);
    width	=	max([1e-5 width]);
    
    h_show(ii)	=	rectangle('Position',[x,y,width,height],...
        'Curvature',[0.3],...
        'LineWidth',2,'LineStyle','-',...
        'EdgeColor','w'); %#ok<AGROW>
    
    %	Highlight selected event with a different color
    if	ie == handles.notes.i_sel
        set(h_show(ii),'EdgeColor', 'm');
        sel_vis		=	true;
        handles.notes.h_info	=	show_event_info(event, Description, handles.notes.h_info);
    end
end

handles.notes.i_show	=	i_show;
handles.notes.h_show	=	h_show;
handles.notes.x_show	=	x_show;
%guidata(h_axes, handles);

%	Enable prev/next buttons
if	length(Events) >= 1
    disable_notes_nav(handles,'on');
else	%	Disable prev/next buttons
    disable_notes_nav(handles,'off');
end

%	Enable edit/delete, if selected event is visible
if sel_vis
    set(handles.pushbutton_notes_edit,'Enable','on');
    set(handles.checkbox_notes_delete,'Enable','on');
else
    set(handles.pushbutton_notes_edit,'Enable','off');
    set(handles.checkbox_notes_delete,'Enable','off');
end

end %plot_events
% --------------------------------------------------------------------

%	Pops up new window with event details
function	hdlg	=	show_event_info(Event, Description, hdlg)

if nargin < 3 || isempty(hdlg) || ~all(ishandle(hdlg))
    hdlg(1)	=	dialog('Windowstyle', 'normal');
    pos		=	get(hdlg(1), 'Position');
    pos(3)	=	pos(3)*.75;
    pos(4)=pos(4)*1.5;
    set(hdlg(1),'Position',pos);
else
    figure(hdlg(1));
    clf(hdlg(1));
end

%	Fetch details
Title	=	'Annotation details';
names	=	fieldnames(Event);
Message	=	[Description(:).'; names(:).'];

%%AARON: some annotations will not be easily presentable, such as
%%'automated' or 'localization'.
Ibad=[];
for	ii	=	1:length(names)
    if isstruct(Event.(names{ii}))
        Ibad=[Ibad ii];
        continue
    end
    if iscell(Event.(names{ii}))
        temp			=	cell2mat(Event.(names{ii})');
        Message{2,ii}	=	[char(9) char(9) temp(:).'];
        
        continue
    end
    
    temp			=	num2str(Event.(names{ii}));
    if size(temp,1)>1
        Ibad=[Ibad ii];
        continue
    end
    Message{2,ii}	=	[char(9) char(9) temp(:).'];
    
end
Igood=setdiff(1:length(names),Ibad);
%	comments don't need to be indented
Message{2,ii}	=	temp;
%	start_time should be formated
Message{2,1}	=	datestr(Event.(names{1}), 'yyyy-mm-dd HH:MM:SS.FFF');

Message=Message(:,Igood);

%	Put details in window
set(hdlg(1), 'Name', Title);
htxt	=	uicontrol(hdlg(1), 'Style', 'edit',...
    'Enable', 'on',...
    'Units', 'normalized',...;
    'Position', [0 0 1 1],...
    'Max', 2, 'Min', 0,...;
    'HorizontalAlignment','left',...
    'String', Message(:));

% If localization info exists, plot it!
if isfield(Event,'localization')&&~isempty(getfield(Event,'localization'))
    
    if nargin < 3 || isempty(hdlg) || ~all(ishandle(hdlg))||length(hdlg)<2
        hdlg(2)	=	dialog('Windowstyle', 'normal');
        pos		=	get(hdlg(2), 'Position');
        pos		=	get(hdlg(2), 'Position');
        pos(1)=pos(1)/2;
        set(hdlg(2),'Position',pos);
        
    else
        figure(hdlg(2));
        clf(hdlg(2));
    end
    
    %Identify which station we are...
    
    Istation=str2num(Event.Istation);
    
    DASAR_coords=[Event.localization.station_position.easting Event.localization.station_position.northing];
    bearings=Event.localization.bearings_all;
    Igood=find(~isnan(bearings));
    
    VM=Event.localization.location;
    A=Event.localization.major;
    B=Event.localization.minor;
    ANG=Event.localization.ellipse_ang;
    [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,35,Istation);
    
    %Compare stored with recomputed results
    %     org.VM=Event.localization.location;
    %     org.A=Event.localization.major;
    %     org.B=Event.localization.minor;
    %     org.ANG=Event.localization.ellipse_ang;
    %
    %     [VM,Qhat,~,outcome] = vmmle_r(bearings(Igood),DASAR_coords(Igood,:),'h');
    %     mean_coords=mean(DASAR_coords(Igood,:));
    %     CRITVAL=4.60517; %chi2inv(0.90,2);
    %
    %     [~,A,B,ANG,Baxis] = ellipsparms(Qhat,CRITVAL,mean_coords,VM);
    %
    %     [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,35,Istation);
    %
    %%used to demonstrate original stored data.
    %disp('Original stored parameters:');
    %org
    %
    %     figure
    %     [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,org.VM,org.A,org.B,org.ANG,35,Istation);
    %     title('original stored parameters');
    %pause
    %close
end

end %show_event_info
% --------------------------------------------------------------------

%	enable/disable any buttons that require an intial update first
function	toggle_initial_buttons(handles, opt)

%	enable audio controls, now that x exists
set(handles.pushbutton_playsound,'Enable',opt);
set(handles.pushbutton_pausesound,'Enable','off'); %always off at first

%	enable new note button
if	~isempty(handles.notes.file_path)
    set(handles.pushbutton_notes_new,'Enable',opt);
end

%	enable next/prev window increment
set(handles.pushbutton_next,'Enable',opt);
set(handles.pushbutton_prev,'Enable',opt);
end  %toggle_initial_buttons
% --------------------------------------------------------------------

%	updates whole figure window according to selected event
function	handles		=	update_events(handles)

i_sel		=	handles.notes.i_sel;
start_time	=	handles.notes.Data.Events(i_sel).start_time;
tlen		=	datenum(0,0,0,0,0,handles.tlen);

%	New start date for window
new_date	=	start_time - tlen/2;
set(handles.edit_datestr, 'String', datestr(new_date,'dd-mmm-yyyy HH:MM:SS.FFF'));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

%	update figure window
pushbutton_update_Callback(handles.pushbutton_update, [], handles);
handles		=	guidata(handles.pushbutton_update);

end %update_events
% --------------------------------------------------------------------

%	Helper function to get current users id/name
function	user_name	=	getusername()
if	isunix || ismac
    user_name = getenv('USER');
elseif ispc
    user_name = getenv('USERNAME');
else
    error('Unregognized system');
end
end %getusername
% --------------------------------------------------------------------

% --- Executes on button press in pushbutton_notes_stats.
function pushbutton_notes_stats_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notes_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


X=[];Y.name=[];
%Provide option to reload a notes file that does not match current load
ReloadButton = questdlg('Would you like to load annotations associated with multiple files?', ...
    'Reload Annotations?');

switch ReloadButton
    case 'Yes'
        
        dialog_title	=	'Select location for "Multiple" annotation files:';
        start_path		=	handles.mydir;
        folder_name		=	uigetdir(start_path, dialog_title);
        
        if isnumeric(folder_name)
            return;
        end
        mydir=handles.mydir;
        
        handles.mydir=folder_name;
        cd(handles.mydir);  %Change location to be inside this folder, since we assume it is a data analysis folder.
        
        if isfield(handles,'myfile')
            myfile_org=handles.myfile;
        else
            myfile_org=[];
        end
        handles.myfile=['Multiple*'];
        handles		=	load_notes_file(handles, folder_name);  %Bulk Load Event detector
        handles.myfile=myfile_org;
        if ~isfield(handles,'tdate_start')
            handles.tdate_start=[];
        end
end


ButtonName = questdlg('What kind of plot?', ...
    'Statistical Analysis of Annotations', ...
    '1-D histogram', '2-D histogram', '2-D boxplot','1-D histogram');


X=get_histogram_vars(ButtonName, handles.notes.Data.Events,handles.notes.Data.Description);

if isempty(X.start)
    return
end
hprint=gobjects(0);
switch ButtonName
    case '1-D histogram'
        for Ix=1:(length(X.start)-1)
            hprint(Ix)=figure;
            edges=X.start(Ix):X.dx:X.start(Ix+1);
            Igood=find(X.var>=X.start(Ix)&X.var<=X.start(Ix+1));
            
            [Nclick,tbin]=histc(X.var(Igood),edges);
            
            
            %%Plot and format the axes
            bar(edges,Nclick,'histc');grid on;
            set(gca,'fontweight','bold','fontsize',14)
            set(gca,'xtick',edges(1):X.label_inc:edges(end));
            if strcmp(X.name,'start_time')
                datetick('x',X.style,'keeplimits','keepticks');
                xlabel('Time');
                ylabel('Counts');
                title(sprintf('%s of %s, %s to %s',ButtonName,X.name,datestr(edges(1)),datestr(edges(end))),'interp','none');
            else
                xlabel(X.label);
                ylabel('Counts')
                title(sprintf('%s of %s, %6.2f to %6.2f',ButtonName,X.name,(edges(1)),(edges(end))),'interp','none');
            end
            xlim([edges(1) edges(end)])
            set(gca,'fontweight','bold','fontsize',14);
            
        end
    case '2-D histogram'
        Y=get_histogram_vars(ButtonName, handles.notes.Data.Events,handles.notes.Data.Description);
        
        if isempty(Y.start)
            return
        end
        %[2 Npoint]
        hprint=zeros(1,length(X.start)-1);
        for Ix=1:(length(X.start)-1)
            
            edges=X.start(Ix):X.dx:X.start(Ix+1);
            edges2=Y.start:Y.dx:Y.start(2);
            Igood=find(X.var>=X.start(Ix)&X.var<=X.start(Ix+1));
            
            XX=[X.var(Igood);Y.var(Igood)];
            
            %[Nclick,tbin]=histc(X.var(Igood),edges);
            [N,printname,Ibin,hh,hprint(Ix)]=hist2D(XX,edges,edges2,{X.label,Y.label},1);
            if strcmp(X.name,'start_time')
                axes(hh(1));
                set(gca,'fontweight','bold','fontsize',14)
                set(gca,'ytick',edges(1):X.label_inc:edges(end));
                datetick('y',X.style,'keeplimits','keepticks');
                
                axes(hh(2));
                set(gca,'fontweight','bold','fontsize',14)
                set(gca,'ytick',edges(1):X.label_inc:edges(end));
                datetick('y',X.style,'keeplimits','keepticks');
            elseif strcmp(Y.name,'start_time')
                axes(hh(1));
                set(gca,'fontweight','bold','fontsize',14)
                set(gca,'xtick',edges2(1):Y.label_inc:edges2(end));
                
                datetick('x',Y.style,'keeplimits','keepticks');
                
                axes(hh(2));
                set(gca,'fontweight','bold','fontsize',14)
                set(gca,'xtick',edges2(1):Y.label_inc:edges2(end));
                
                datetick('x',Y.style,'keeplimits','keepticks');
                
            end
        end
    case '2-D boxplot'
        Y=get_histogram_vars(ButtonName, handles.notes.Data.Events,handles.notes.Data.Description);
        
        if isempty(Y.start)
            return
        end
        hprint=zeros(1,length(X.start)-1);
        
        for Ix=1:(length(X.start)-1)
            
            edges=X.start(Ix):X.dx:X.start(Ix+1);
            edges2=Y.start:Y.dx:Y.start(2);
            Igood=find(X.var>=X.start(Ix)&X.var<=X.start(Ix+1));
            
            XX=[X.var(Igood);Y.var(Igood)];
            
            %[Nclick,tbin]=histc(X.var(Igood),edges);
            %[N,printname,Ibin,hh]=hist2D(XX,edges,edges2,{X.label,Y.label},1);
            hprint(Ix)=plot_data_boxplot_percentile(X.var(Igood),X.dx,Y.var(Igood),[Y.start(1) Y.start(end)],X.label_inc,[],X.style);
            %hprint=plot_data_boxplot_percentile(Tnew, time_inc,sumPSD,ylimits,xlabel_inc,percentile);
            
            if strcmp(X.name,'start_time')
                %datetick('x',X.style,'keeplimits','keepticks');
                xlabel('Time');
                ylabel(Y.label);
                title(sprintf('%s of %s, %s to %s',ButtonName,X.name,datestr(edges(1)),datestr(edges(end))),'interp','none');
            else
                xlabel(X.label);
                ylabel(Y.label)
                title(sprintf('%s of %s, %6.2f to %6.2f',ButtonName,X.name,(edges(1)),(edges(end))),'interp','none');
            end
            
        end
        
    otherwise
end

%%Plot results

ButtonName2 = questdlg('Print and save Data?', ...
    'Save stats...');
switch ButtonName2
    case 'No'
        return
end
%Nfigs=sort(get(0,'Child'));
%Nfigs=Nfigs(Nfigs<150);

%for II=1:length(Nfigs)
try
    for Ix=1:length(hprint)
        edges=X.start(Ix):X.dx:X.start(Ix+1);
        
        if strcmp(X.name,'start_time')&&isempty(Y.name)
            save_str=sprintf('Stats_%s_%s_%s_%s', ButtonName,X.name,datestr(edges(1),30),datestr(edges(end),30));
        elseif isempty(Y.name)
            save_str=sprintf('Stats_%s_%s_%s_%s', ButtonName,X.name,(edges(1)),(edges(end)));
        elseif strcmp(X.name,'start_time')&&~isempty(Y.name)
            save_str=sprintf('Stats_%s_%s_%s_%s_%s', ButtonName,X.name,Y.name,datestr(edges(1),30),datestr(edges(end),30));
        elseif ~isempty(Y.name)
            save_str=sprintf('Stats_%s_%s_%s_%s_%s', ButtonName,X.name,Y.name,(edges(1)),(edges(end)));
            
        end
        
        orient landscape
        print(hprint(Ix),'-djpeg','-r300',save_str);
        saveas(hprint(Ix), save_str, 'fig');
        %end
        close(hprint(Ix));
    end
    
    
    if strcmp(X.name,'start_time')&&isempty(Y.name)
        save_str=sprintf('Stats_%s_%s_%s_%s', ButtonName,X.name,datestr(X.start(1),30),datestr(X.start(end),30));
    elseif isempty(Y.name)
        save_str=sprintf('Stats_%s_%s_%s_%s', ButtonName,X.name,X.start(1),X.start(end));
    elseif strcmp(X.name,'start_time')&&~isempty(Y.name)
        save_str=sprintf('Stats_%s_%s_%s_%s_%s', ButtonName,X.name,Y.name,datestr(X.start(1),30),datestr(X.start(end),30));
    elseif ~isempty(Y.name)
        save_str=sprintf('Stats_%s_%s_%s_%s_%s', ButtonName,X.name,Y.name,X.start(1),X.start(end));
        
    end
    
    save(save_str,'X','Y');
    
    uiwait(msgbox([save_str ' mat, fig and jpg file written to ' pwd],'replace'));
catch
    uiwait(errdlg('Figures have been deleted','Can''t save figures'));
    
end
end %pushbutton_notes_stats_Callback


%%%%%%%%%%%load_and_display_spectrogram%%%%
function	handles	=	load_and_display_spectrogram(handles)

cla;

%tdate_start		=	handles.tdate_start;
tlen	=	handles.tlen;
%if tlen<1e-2 %somehow window length is way too short
%    tlen=1;
%   set(handles.edit_winlen,'String',num2str(tlen));
%end
if isempty(tlen)
    tlen    =   str2double(get(handles.edit_winlen,'String'));
    
    handles.tlen=tlen;
end
mydir	=	pwd;
Ichan	=	eval(get(handles.edit_chan,'String'));

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
    
    contents=get(handles.popupmenu_ovlap,'String');
    ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
    ovlap=min([1-1/Nfft ovlap]);
    
    
end  %if PSD

handles.old_display_view=[];
if isfield(handles,'display_view')
    handles.old_display_view=handles.display_view;
end

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if strcmp(handles.display_view,'Spectrogram')||strcmp(handles.display_view,'New Fig')
    
    if strcmp(handles.display_view,'Spectrogram')
        axes(handles.axes1);
    else
        figure;
    end
    
    if length(x(:,1))<Nfft/2
        return
    end
    if ~(strcmp(handles.filetype,'PSD'))
        [S,FF,TT,B] = spectrogram(x(:,1),hanning(Nfft),round(ovlap*Nfft),Nfft,Fs);
        %B=(2*abs(B).^2)/(Nfft*Fs); %Power spectral density...
        %     For real signals, variable 'B'
        %     returns the one-sided modified periodogram estimate of the PSD of each
        %     segment; for complex signals and in the case when a vector of
        %     frequencies is specified, it returns the two-sided PSD.
        handles.sgram.T		=	TT;
        handles.sgram.F		=	FF;
        handles.sgram.B		=	B;
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
            
            imagesc(TT,FF/1000,10*log10(B));%
        end
    else
        
        ppsd=10*log10(x);
        imagesc(t,FF/1000,ppsd);
        handles.sgram.T		=	t;
        handles.sgram.F		=	FF;
        handles.sgram.B		=	ppsd;
        handles.sgram.Nfft	=	hdr.Nfft;
        handles.sgram.ovlap	=	ovlap;
        handles.sgram.Fs	=	Fs;
    end
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
    if get(handles.checkbox_grayscale,'Value')==1,
        colormap(flipud(gray));
    else
        colormap(jet);
    end
    colorbar;
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
    if ~strcmp(handles.display_view,'Spectrogram')
        title(get(handles.text_filename,'String'));
    end
elseif strcmp(handles.display_view,'Time Series') %%Time series
    
    
    %%Check that we are looking at acoustic data and are in spectrogram
    %%mode
    if isempty(strfind(handles.myfile,'Press'))
        % msgbox('Enter min and max frequency, or hit return to skip filtering:','modal');
        if strcmp(handles.old_display_view,'Spectrogram')
            ButtonName = questdlg('Filter? (If yes, click on two frequency values in spectrogram)');
        else
            ButtonName='No';
        end
        
        %if ~isempty(tmp)&&size(tmp,1)==2
        %%Check that are on current spectrogram view
        if strcmp(ButtonName,'Yes')
            tmp=ginput(2);
            
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
        y=x(:,1)/1000; %Pressure conversion
    end
    
    t=(1:length(x(:,1)))/Fs;
    xlabel('Time (sec)');
    if isfield(hdr,'calunits')&&strfind(hdr.calunits,'mPa')
        plot(handles.axes1,t,1000*y);grid on;
        ylabel('uPa');
    else
        plot(handles.axes1,t,y);grid on;
        ylabel('Amplitude');
    end
elseif strcmp(handles.display_view,'Correlogram') %%Correlogram
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

if handles.file_flags.multichannel
    status='on';
else
    status='off';
end
for II=1:length(handles.buttongroup.array)
    set(handles.buttongroup.array(II),'vis',status);
    set(handles.buttongroup.array(II),'enable',status);
end


if strcmpi(lower(handles.filetype),'gsi')
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
    
end
end %function load_and_display_spectrogram

function [handles,errorflag]	=	set_slider_controls(handles,filetype)
%%Set min, max and other times associated with slider controls
mydir			=	pwd;
errorflag		=	0;
handles.tdate_min=	-1;
handles.tdate_max=	-1;
%
% if isempty(handles.mydir)
%     handles.mydir=pwd;
% end
%cd(handles.mydir);

try
    [x,t,Fs,tmin,tmax]=load_data(filetype,-1,10,1,handles);
catch %no file selected
    %errordlg(sprintf('No %s file selected',filetype));
    errorflag=1;
    cd(mydir);
    return
end


set(handles.slider_datestr,'Min',0);
set(handles.slider_datestr,'Max',1);

%	5sec and 10% increments
T_len		=	(tmax-tmin)*24*60*60;
big_step	=	0.1;					% 10%
small_step	=	0.01;				% 5s

set(handles.slider_datestr,'sliderstep',sort([small_step big_step]));

set(handles.slider_datestr,'Value',0.0);
%handles.tdate_start		=	0.5*(tmin+tmax);
handles.tdate_start		=	tmin;
set(handles.edit_datestr,'String',datestr(handles.tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'));

set(handles.text_mintime,'String',datestr(tmin,0));
set(handles.text_maxtime,'String',datestr(tmax,0));
handles.tdate_min	=	tmin;
handles.tdate_max	=	tmax;

currentFs=str2num(get(handles.edit_fmax,'String'));
if currentFs==0
    set(handles.edit_fmax,'String',Fs/2000);
    set(handles.edit_fmin,'String',0);
end
%slider_step(1) = datenum(0,0,0,0,0,handles.tlen)/(maxx-minn);
%slider_step(2) = min([datenum(0,0,0,0,5,0) 0.1*(maxx-minn)])/(maxx-minn);
%set(handles.slider_datestr,'sliderstep',slider_step)
% keyboard
cd(mydir);

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
elseif strfind(local_machine,'Jan-Straleys')
    
    startup_info.base_directory='/Users/janstraley/SEASWAP_2011';
    
    
    startup_info.default_directory='/Users/janstraley/SEASWAP_2011';
    startup_info.default_inputfiledir='/Users/janstraley/SEASWAP_2011';
    
    %%Default annotation file
    startup_info.annotation_file='annotated.txt';
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
    startup_info.calibration_DASAR2007_dir='';
elseif strfind(local_machine,'Janice')
    
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
    
    
elseif ~isempty(strfind('dgrebner',local_machine))
    
    startup_info.base_directory='/Users/dgrebner/Desktop/';
    startup_info.default_directory=startup_info.base_directory;
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    startup_info.annotation_file='annotated.txt';
    startup_info.calibration_DASAR2007_dir='';
    
    startup_info.function_handles_filename='mt_specgram_handles.mat';
    
    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
elseif ~isempty(strfind('thode',local_machine))
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function load_data.m%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,t,Fs,tmin,tmax,head]	=	...
    load_data(filetype,tdate_start,tlen,Ichan,handles)
%%% tmin,tmax, t are datenumbers
%   x rows are samples, columns are channels
%   head.geom.rd:  If multiple channels are present, this gives spacing
persistent  Fs_keep keyword space

mydir=handles.mydir;
myfile=handles.myfile;
teager=get(handles.checkbox_teager,'Value');
set(handles.edit_normal_rotation,'Vis','Off');  %Set off for the moment
set(handles.text_normal_rotation,'Vis','Off');  %Set off for the moment

Fs=[];
x=[];
t=[];
tmin=[];
tmax=[];
head=[];
head.multichannel=false;
head.linked=false;

filetype	=	upper(filetype);

switch filetype
    case 'PSD'
        [x,F,t,Tabs,params]=read_Java_PSD(fullfile(mydir,myfile),tdate_start,tlen);
        if ~isempty(t)&&length(t)>1
            Fs=1./(t(2)-t(1));
            t=t-t(1);
        else
            Fs=1;
        end
        head=params;
        
        tmin=head.tstart_file;
        tmax=head.tend_file;
    case 'MAT'
        
        simulated=load(fullfile(mydir,myfile));
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
        head.geom.rd=simulated.rd;
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
        head.multichannel=false;
        
        if isempty(keyword)
            prompt = {'Enter a keyword for GSI calibration [DASARC]:'};
            dlg_title = 'DASAR calibration';
            num_lines = 1;
            def = {'DASARC'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            keyword=answer{1};
        end
        % keyword=input('Enter a keyword for GSI calibration [DASARC]:','s');
        %             if isempty(keyword)
        %                 keyword='DASARC';
        %             end
        %         end
        x=calibrate_GSI_signal(x, keyword);
        
        Fs=head.Fs;
        
        tmin=datenum(1970,1,1,0,0,head.ctbc);
        tmax=tmin+datenum(0,0,1,0,0,0);
        head.Nchan=length(Ichan);
        head.linked=true;
        
    case 'MT'
        %[x,t,Fs]=load_mt_mult(handles.mydir,tdate_start,tlen);
        %head=read_mt_header([mydir filesep myfile]);
        head=read_mt_header(fullfile(mydir, myfile));
        
        tmin=head.tstart;
        tmax=head.tend;
        Fs=head.Fs;
        
        if tdate_start>0
            tdate_vec=datevec(tdate_start-tmin);
            nsec=tdate_vec(6)+60*tdate_vec(5)+3600*tdate_vec(4);
        else
            %uiwait(msgbox('No start time requested'));
            head.multichannel=false;
            head.linked=false;
            
            return
        end
        %Load accelerometer data as desired
        %% Template: CB_Accel_X_2014_02100958.mt
        if strfind(head.abbrev,'Accel')
            mystr='XYZ';
            Ispace=strfind(myfile,'_');
            Ispace=Ispace(2)+1;
            x=[];
            for Iacc=1:3
                myfile_temp=myfile;
                myfile_temp(Ispace)=mystr(Iacc);
                try
                    [x0,t]=load_mt([mydir filesep myfile_temp],nsec,tlen);
                    if isempty(x)
                        x=zeros(length(x0),3);
                    end
                    x(:,Iacc)=x0;
                catch
                    fprintf('Channel %s of %s not available...\n',mystr(Iacc),myfile_temp);
                    
                end
                
            end
            head.Nchan=3;
            
            %%Select channel requested...or rotate channels...
            % Axis conventions assume 2011 definitons: X is long axis, Y
            % and Z are normal.
            
            set(handles.edit_normal_rotation,'Vis','On');
            set(handles.text_normal_rotation,'Vis','On');
            normal_angle=(pi/180)*str2num(get(handles.edit_normal_rotation,'String'));
            
            if normal_angle~=0
                xrot=x;
                xrot(:,2)=x(:,2)*cos(normal_angle)+x(:,3)*sin(normal_angle);
                xrot(:,3)=-x(:,2)*sin(normal_angle)+x(:,3)*cos(normal_angle);
                x=xrot;
                head.normal_rotation=(180/pi)*normal_angle;
                head.transformation.description{1}=' Normal Rotation: %6.2f deg';
                head.transformation.value{1}=head.normal_rotation;
            else
                head.transformation.description{1}=' No rotation';
                head.transformation.value{1}=[];
            end
            
            x=x(:,Ichan);
            %myfile(Ispace)=mystr(Ichan);
            
            head.myfile=myfile;
            
            % guidata(handles.myfile);  %Saves new file name...
            %guidata(hObject, handles);  %Saves new file name...
            %guidata(hObject, handles);
            %
        else
            if Ichan>1
                disp('WARNING: load_data: MT can only have one channel');
                Ichan=1;
            end
            if tdate_start>0
                [x,t]=load_mt(fullfile(mydir , myfile),nsec,tlen);
            end
            head.Nchan=1;
        end
    case 'SIO'
        %x,t,Fs,tmin,tmax,head
        %load_data(filetype,tstart_min,tdate_start,tlen,Ichan,handles)gff
        %set(handles.text_channel,'String','Channel [-angle (deg)]');
        sio_chc=get(handles.togglebutton_ChannelBeam,'String');
        [~, head] = sioread(fullfile(mydir,myfile));
        head.multichannel=true;
        
        %Extract start date and time
        %if ~exist('tstart_min') || isempty(tstart_min) || tstart_min < 0
        Idot=strfind(myfile,'.');
        
        if isfield(head,'date')
            tmin=head.date;
        else
            success=0;
            for II=1:(length(Idot)-1)
                if success || Idot(II+1)-Idot(II)~= 12
                    continue
                end
                template=myfile((Idot(II)+1):(Idot(II+1)-1));
                try
                    tmin=datenum(2000+str2num(template(1:2)),0,str2num(template(3:5)), ...
                        str2num(template(6:7)),str2num(template(8:9)),str2num(template(10:11)));
                    %tstart_min=tmin;
                    success=1;
                catch
                    tmin= now;
                end
            end
        end
        %end
        
        if isfield(head,'Fs')
            Fs=head.Fs;
        else
            Fs_default=25000;
            if isempty(Fs_keep)
                Fs=input('Enter a sampling frequency in Hz per channel [25000]:','s');
                if isempty(Fs),Fs=Fs_default;end
                head.Fs=Fs;
                Fs_keep=Fs;
            else
                Fs=Fs_keep;
            end
        end
        
        %%	Data parameters needed
        head.np		=	head.PperChan;
        tmax =tmin+datenum(0,0,0,0,0,head.np/Fs);
        head.Nchan	=	head.N_Chan;
        
        if ~exist('tdate_start') || isempty(tdate_start) || tdate_start < 0
            x=[];
            return
        end
        
        if tdate_start>0
            np_start=1+round((tdate_start-tmin)*24*3600*Fs);
        else
            
            np_start=1;
        end
        if np_start>head.np
            errordlg('load_data:  SIO start time exceeds points available!');
            return
        end
        
        
        %%%If beamforming is desired for a look at a given direction...
        beamform_data=0;
        get_geometry=0;
        if strcmpi(sio_chc,'angle')&&~strcmpi(Ichan,'all')
            beamform_data=1;
            get_geometry=1;
        elseif strcmpi(Ichan,'all')
            Ichan=1:head.Nchan;
            beamform_data=0;
            get_geometry=1;
        end
        
        if beamform_data==1
            thta=-Ichan;
            Ichan=1:head.Nchan;
        end
        
        %If loading single channel, check that request is reasonable...
        if beamform_data==0
            if max(Ichan)>head.Nchan
                errordlg(sprintf('load_data:  SIO channel request: %i, max channels: %i',max(Ichan),head.Nchan));
                return
            end
            
            if max(Ichan)<1
                errordlg(sprintf('load_data:  SIO channel request: %i is less than 1',max(Ichan)));
                return
            end
            
        end
        
        npi=round(tlen*Fs);
        if np_start+npi>head.np
            errordlg('load_data:  SIO end time exceeds points available!');
            return
        end
        
        
        
        if get_geometry==1
            if isempty(space)
                prompt = {'Enter spacing[m] between elements for SIO file:'};
                dlg_title = 'SIO file spacing';
                num_lines = 1;
                def = {'0.1'};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
                space=eval(answer{1});
                fprintf('Half-wavelength frequency: %6.2f Hz\n',1500/(2*space));
            end
            head.geom.rd=(0:(head.Nchan-1))*space;
            
            
        end
        
        
        [x,~]=sioread(fullfile(mydir,myfile),np_start,npi,Ichan);
        %Data arranged so that time are rows, columns are channels
        
        %Flip data ....
        
        if size(x,2)>1
            x=fliplr(x);
        end
        
        if beamform_data==1
            try
                space=head.geom.rd(2)-head.geom.rd(1);
                xtot=delaynsum(x,thta,space,Fs,Ichan);
                x=xtot;
                head.thta=thta;
            catch
                disp('load_data: sioread beamform failure');
                keyboard
            end
        end
        
        %Set time
        t=(1:size(x,1))/Fs;
        
        
    case 'DAT'
        [x,tmin,tmax,fs]=read_dat_file(fullfile(mydir,myfile),[],-1,tlen,0);
        
        if isempty(fs)
            fs=input('Enter sampling rate in Hz:');
        end
        Fs=fs;
        [x,tmin,tmax]=read_dat_file([mydir '/' myfile],Fs,tdate_start,tlen,0); %output in uPa
        %[x,tmin,tmax]=read_dat_file([mydir '/' myfile],Fs,tdate_start,tlen,1);  %Voltage output
        t=(1:length(x))/fs;
        
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
        t=(1:length(x))/fs;
        head.Nchan=1;
        
    case 'MDAT'
        
        
        [x,head]=read_synchronized_mdat_files(fullfile(mydir,myfile),tdate_start,tlen);
        Fs=head.fs;
        sio_chc=get(handles.togglebutton_ChannelBeam,'String');
        head.multichannel=true;
        t=(1:length(x))/Fs;
        
        tmin=head.tfs;
        tmax=head.tfe;
        head.Nchan=size(x,2);
        beamform_data=0;
        get_geometry=1;
        if strcmpi(sio_chc,'angle')&&~strcmpi(Ichan,'all')
            beamform_data=1;
            get_geometry=1;
        elseif strcmp(Ichan,'all')
            Ichan=1:head.Nchan;
            beamform_data=0;
            get_geometry=1;
        end
        
        if beamform_data==1
            thta=-Ichan;
            Ichan=1:head.Nchan;
        end
        
        if max(Ichan)>head.Nchan
            disp(sprintf('Channel too high, max channel %i',head.Nchan));
            Ichan=head.Nchan;
        end
        
        try
            x=x(:,Ichan);
            for II=Ichan
                fprintf('Logged depth of channel %i is %6.2f m\n',II,head.geom.rd(II));
            end
        end
        
        x=x';
        
        if beamform_data==1
            try
                space=head.geom.rd(2)-head.geom.rd(1);
                xtot=delaynsum(x',thta,space,Fs,Ichan);
                x=xtot;
                head.thta=thta;
            catch
                disp('load_data: MDAT beamform failure');
                keyboard
            end
        end
        
        %Set time
        t=(1:length(x))/Fs;
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
        
    case 'M4V'
        info	=	audioinfo(fullfile(mydir,myfile));
        info_vid = mmfileinfo(fullfile(mydir,myfile));
        %[~,Fs]		=	audioread(fullfile(mydir,myfile),1,'native');
        Nsamples	=	info.TotalSamples;
        handles.Fs	=	info.SampleRate;
        Fs=info.SampleRate;
        tmin=datenum([1970 1 1 0 0 0]);
        tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
        %tlen=Nsamples/handles.Fs;
        sens=1;
        
        tdate_vec	=	datevec(tdate_start - tmin);
        nsec		=	tdate_vec(6) + 60*tdate_vec(5) + 3600*tdate_vec(4);  %Ignores differences in days
        N1			=	1 + round(nsec*handles.Fs);
        N2			=	N1 + round(tlen*handles.Fs);
        
        try
            [x,Fs]		=	audioread(fullfile(mydir,myfile),[N1 N2],'native');
            
            
        catch
            x=[];
            
            t=[];
            head.Nchan=0;
            return
        end
        
        head.Nchan	=	size(x,2);
        if head.Nchan>1
            head.multichannel=true;
        end
        
        %%Option for "adaptive processing" to remove interference...
        % Set channel equal to -1 to give difference between signal
        
        if ~strcmp(Ichan,'all')&Ichan>0
            x		=	x(:,Ichan);
        elseif Ichan<0
            x=x(:,2)-x(:,1);
        end
        
        t	=	(1:length(x))/Fs;
        
        x			=	double(x)*sens;
    case 'WAV' %% Includes SUDAR data
        
        
        info	=	audioinfo(fullfile(mydir,myfile));
        %[~,Fs]		=	audioread(fullfile(mydir,myfile),1,'native');
        Nsamples	=	info.TotalSamples;
        handles.Fs	=	info.SampleRate;
        Fs=info.SampleRate;
        
        [head.cable_factor,sens]=get_ADAT24_cable_factor;
        
          
        try
            [SUDAR_true,tmin,tmax,Fs]=get_SUDAR_time(mydir,myfile); %Check whether a sUDAR file exists
            if SUDAR_true
                sens=(10^(186/20))/(2^15);
                
            end
            if ~SUDAR_true
                tmin	=	convert_date(myfile,'_');
                if isempty(tmin)
                    tmin=datenum([1970 1 1 0 0 0]);
                end
                tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
            end
        catch
            disp([myfile ': convert_date failure']);
            try
                tmin=datenum(get(handles.text_mintime,'String'));
            catch
                minn	=	input('Enter start date in format [yr mo day hr min sec] or hit return: ');
                if isempty(minn)
                    minn=[1970 1 1 0 0 0];
                end
                tmin	=	datenum(minn);
            end
            tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
        end
        
        
        tdate_vec	=	datevec(tdate_start - tmin);
        nsec		=	tdate_vec(6) + 60*tdate_vec(5) + 3600*tdate_vec(4);  %Ignores differences in days
        N1			=	1 + round(nsec*handles.Fs);
        N2			=	N1 + round(tlen*handles.Fs);
        
        try
            [x,Fs]		=	audioread(fullfile(mydir,myfile),[N1 N2],'native');
            
            
        catch
            x=[];
            
            t=[];
            head.Nchan=0;
            return
        end
        
        head.Nchan	=	size(x,2);
        if head.Nchan>1
            head.multichannel=true;
        end
        if ~strcmp(Ichan,'all')
            x		=	x(:,Ichan);
        end
        
        t	=	(1:length(x))/Fs;
        
        x			=	double(x)*sens;
        
end

if isempty(tmin) || isempty(tmax)
    disp('load_data: Warning, tmin and tmax should never be empty when exiting..');
end


%%Store whether multichannel data, regardless of original file format.
if min(size(x))>1
    head.multichannel=true;
end

if ~isfield(head,'linked')
    head.linked=false;
end

if ~isfield(head,'multichannel')
    head.multichannel=false;
end

%%%Optional Teager-Kaiser filtering...
if teager
    %%Assume that x is in form [ channel time]
    x=x(:,2:end-1).^2-x(:,1:end-2).*x(:,3:end);
end


end  %function load_data

function [cable_factor,sens]=get_ADAT24_cable_factor
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

cable_factor=2*pi*(110e-9)*140.0;  %Unit resistance 140 ohm, capacitance 110 nF
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

if strfind(rawfile,'.sio')
    
    [y,t,head]=readsiof(rawfile,cbegin,tlen,formatt);
    y=y-0.5*(2^16);  %Remove DC bias in A/D converter
    
elseif strfind(rawfile,'.gsi')
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

end  %function readGSIfile

function [thet,kappa,tsec]=get_GSI_bearing(hObject,eventdata,handles,tmp)
% Freq or temp in kiloHz
%tsec: seconds into time series that data are selected...
% tmp is previous ginput
thet=-1;  %Start with failed result
kappa=-1;
tsec=-1;
tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
%Ichan=str2double(get(handles.edit_chan,'String'));
[x,t,Fs,tstart,tend,head]=load_data(handles.filetype,tdate_start,tlen,'all',handles);

disp('Click on two extreme corners, click below axis twice to reject:');
if ~exist('tmp','var')
    tmp=ginput(2);
end
if (isempty(tmp)||any(tmp(:,2)<0))
    return
end
tsec=min(tmp(:,1));
n=round(Fs*sort(tmp(:,1)));

freq=1000*sort(tmp(:,2));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

[thet0,kappa,sd]=extract_bearings(x(n(1):n(2),:),0.25,Nfft,Fs,freq(1),freq(2),200);

if ~isempty(strfind('T2007',handles.myfile))
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
%Iend=strfind(titlestr,'.gsi')-1;
set(handles.text_filename,'String',sprintf('%s/%s %6.2f degrees... ',handles.mydir,handles.myfile,thet));
%keyboard;
end %function get_GSI_bearing

function nhout=read_num(fid,N)
nhout=str2double(char(fread(fid,N,'char')'));
end

function chout=read_char(fid,N)
chout=char(fread(fid,N,'char')');
end

function tabs=parse_date(str)

Istart=strfind(str,'=');
year=str2double(str((end-4):end));
tm=datenum(str((end-13):(end-5)),14)-datenum('00:00:00',14);
day=str2double(str((end-15):(end-14)));
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

function [thet,kappa,sd]=extract_bearings(y,bufferTime,Nfft,Fs,fmin,fmax,Nsamples)
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




function bearings=bnorm(bearings)
Ibig=find(bearings>=360);
bearings(Ibig)=bearings(Ibig)-360;

Ilow=find(bearings<0);
bearings(Ilow)=bearings(Ilow)+360;

end

function [Ksout,Ns,EE_sort,VV]=extractKsbest_contour(x,ovlap,Nfft,chann, frange,fr,fbad,Fs,M,keep_zeros,nowin)
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
EE_sort=[];VV=[];
figure;
if ~exist('nowin', 'var'),
    nowin=0;
elseif  nowin==1
    nowin=1;
else
    nowin=0;
end

Nel=length(chann);
if ~exist('M', 'var')||M<0, M=Nfft;end

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
    if ~isempty(fbad)
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

if exist('keep_zeros', 'var') %%Keep all frequency bins, even if no power...
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

function B=MV_beamforming(Ks,angles,freq,Lz,c)
B=zeros(length(freq),length(angles));
if size(Lz,2)>1,
    Lz=Lz.';
end
for If=1:length(freq),
    for Iang=1:length(angles)
        % lambda=1500/freq(If);
        w=exp((1i*2*pi*Lz*freq(If)/c)*sin(angles(Iang)*pi/180));
        
        %B(If,Iang)=real(w'*Ks{If}*w)/(norm(Ks{If})*norm(w).^2);
        w=w/norm(w);
        B(If,Iang)=1./real(w'*(squeeze(Ks(:,:,If))\w));
        
    end
    
    
end

%plot(angles,10*log10(B(If,:)));
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


%% MFP
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
    disp(sprintf('Frequencies to process, hard wired: %s',mat2str(freq,4)));
    if length(unique(Igood))~=length(Igood),disp('WARNING! Frequencies requested via SNRmin not available...');end
end

Kstot=Kstot(:,:,Igood);
disp(sprintf('Loading model: %s',model_name));
model=load(model_name);
prompt={'Number of modes to model:'};
name='Parameters for Matched Field Processing';
numlines=1;
defaultanswer={'all'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if strfind(answer{1},'all')
    max_mode= size(model.kr{1},2);
    
else
    max_mode=min([eval(answer{1}) size(model.kr{1},2)]);
end
model.U=model.U(:,:,1:max_mode);
model.kr{1}=model.kr{1}(:,1:max_mode);
%keyboard

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
            
            Ibad = isnan(real(phase));
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
    colorbar('fontweight','bold','fontsize',14);
    %title(sprintf('Frequencies: %s',mat2str(round(freq),4)));
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
        disp(savename)
        save([savename '.mat'], 'ranges','depths','amb','amb_tot','SNRmin','freq','data_name','model_name');
        
        figure(gcf)
        orient landscape
        
        print(gcf,'-djpeg','-r300',[savename '.jpg']);
    end
end %tilt
end %matched_field_processor

%% Correlograms
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
end %create_incoherent_correlogram

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

%Trim away portions of autocorrelation that will not have creaks (ici_range)
Iwant=find(tt>ici_range(1)&tt<ici_range(2));

%Remove times not of interest, normalize zero lag to one.
XC=XC(Iwant,:);
XC_eq=XC_eq(Iwant,:);
%XC0=XC0(Iwant,:);
%signn=signn(Iwant,:);
tt=tt(Iwant);
%Quick click removal...





end  %create_coherent_correlogram

%% Ellipse
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

if ischar(C),C=C(:);end;

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
        an=ang(k);
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
if isnan(S(1)) || (detS<=0) || min(diag(S))<0,
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
    if ((theta>-pi/2) && (theta<pi/2)  && (ang>(pi/2+theta))), % Adjust angle for
        ang = ang-pi;                                          %   estimated position
    elseif ((theta>=pi/2) && (ang<(theta-pi/2))),             %   relative to center
        ang = ang+pi;                                          %   of DASAR array.
    elseif ((theta<=-pi/2) && (ang<(3*pi/2+theta))),
        ang = ang+pi;
    end
end
end

%% write_covmat

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
end %write_covmat

%% Audio files
function audio_stop(player, ~, handles)
%	Set controls
ispaused	=	(player.CurrentSample > 1) &&...
    (player.CurrentSample < player.TotalSamples);
if	~ispaused
    set(handles.pushbutton_playsound, 'String', 'Play');
    set(handles.pushbutton_pausesound, 'String', 'Pause');
    set(handles.pushbutton_pausesound, 'Enable', 'off');
end

end %audio_stop


function audio_timer(player, ~, handles)

%	Get handle to marker line
hline	=	handles.hline;

% check if sound is playing, then only plot new marker
if strcmp(player.Running, 'on')
    % get the currently playing sample #
    x	=	player.CurrentSample;
    Fs	=	handles.audioFs;
    
    % change position of marker line
    set(hline,'XData',x/Fs*[1 1]);
    drawnow expose;
end
end %audio_timer

% --------------------------------------------------------------------

%%
%[X.var,X.start,X.dx,X.label_inc,X.label,X.name,X.style]=get_histogram_vars(ButtonName,handles.notes.Data.Events);

function X=get_histogram_vars(ButtonName,Events,Description)
%Input:  ButtonName: horizontal array of data
%         Events: annotation date..
%Output:  [X.var,X.start,X.dx,X.label_inc,X.label,X.name,X.style]

X.start=[];
X.dx=[];
X.label_inc=10;
X.style=[];
Nx=length(Events);

var_names=fieldnames(Events(1));
titstr=sprintf('Select an X variable for %s',ButtonName);
I_name=menu(titstr,var_names);
X.name=var_names{I_name};
X.label=Description{I_name};


switch X.name
    case {'author','comments','sig_type','call_type'}
        uiwait(msgbox('Inappropriate variable selected','Bad variable!','modal'));
        return
    case 'start_time'
        X.var=([(Events.(X.name))]);
        
        tspan=datevec(X.var(end)-X.var(1));
        %Make initial guess based on span
        if tspan(3)+tspan(4)/24>0.9 %If greater than a day
            interval=60; % one hour increments
        elseif tspan(4)+tspan(5)/60>0.75  %If greater than 45 minutes
            interval=10; %10 minute increments
        elseif tspan(4)+tspan(5)/60>0.5  %If greater than 45 minutes
            interval=2; %2 minute interval
        else
            interval=1; %1 minute interval
        end
        
        dinterval=datenum(0,0,0,0,interval,0);
        xmin=floor(min(X.var)/dinterval)*dinterval;
        xmax=ceil(max(X.var)/dinterval)*dinterval;
        
        prompt1={'Start Time','End Time', ...
            'Time per plot, minutes("all" puts all data in single window)','bin width (minutes)','tick label increment (minutes)'};
        def1={datestr(xmin), datestr(xmax),'all',num2str(interval,2),num2str(2*interval,4)};
        
        dlgTitle1=sprintf('Parameters for %s',X.name);
        answer=inputdlg(prompt1,dlgTitle1,1,def1);
        try
            tabs(1)=datenum(answer{1});
            tabs(2)=datenum(answer{2});
            X.dx=datenum(0,0,0,0,str2num(answer{4}),0);
            
        catch
            uiwait(msgbox('Entry variables not in datenumber format','Bad entry!','modal'));
            return
        end
        
        
        if strcmp(answer{3},'all')
            X.start=[tabs(1) tabs(2)];
        else
            x_window=datenum(0,0,0,0,eval(answer{3}),0);
            X.start=unique([tabs(1):x_window:tabs(2) tabs(2)]);
        end
        
        %%Create tick labels
        try
            if str2num(answer{5})>60
                X.style='HH';
            else
                X.style=15;
            end
            X.label_inc=datenum(0,0,0,0,str2num(answer{5}),0);
        catch
            X.label_inc=datenum(0,0,0,0,str2num(answer{4}),0);
        end
        
        %Find approximate spacing of tick labels
        
        
        
    otherwise
        X.var=zeros(1,Nx);
        
        for Ix=1:Nx
            X.var(Ix)=str2num(Events(Ix).(X.name));
        end
        Nbins=20;
        %X.var=sort((X.var));
        X.dx=round(100*((max(X.var)-min(X.var))/Nbins))/100;
        
        xmin=floor(min(X.var)/X.dx)*X.dx;
        xmax=ceil(max(X.var)/X.dx)*X.dx;
        prompt1={'Start Value','End Value', ...
            'bin width ','tick label increment '};
        
        def1={num2str(xmin), num2str(xmax),num2str(X.dx),num2str(2*X.dx)};
        dlgTitle1=sprintf('Parameters for %s ',X.name);
        answer=inputdlg(prompt1,dlgTitle1,1,def1);
        try
            xx(1)=str2num(answer{1});
            xx(2)=str2num(answer{2});
            X.dx=str2num(answer{3});
            
        catch
            uiwait(msgbox('Entry variables not numeric!','Bad entry!','modal'));
            return
        end
        X.start=[xx(1) xx(2)];
        %%Create tick labels
        try
            X.label_inc=str2num(answer{4});
        catch
            X.label_inc=str2num(answer{3});
        end
        
        
end



end %get_histogram_vars
% --------------------------------------------------------------------

%% Accelerometer rotation files
function edit_normal_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_normal_rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_normal_rotation as text
%        str2double(get(hObject,'String')) returns contents of edit_normal_rotation as a double
end %edit_normal_rotation_Callback
% --------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function edit_normal_rotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_normal_rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end %edit_normal_rotation_CreateFcn

% --------------------------------------------------------------------
function image_processor_Callback(hObject, eventdata, handles)
% hObject    handle to image_processor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

param.morph=[];
cbegin=3600*24*(handles.tdate_start-datenum(1970,1,1,0,0,0));

chc_list=    load_image_parameters(1);

param.Fs=handles.Fs;
contents=get(handles.popupmenu_Nfft,'String');
param.Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});
contents=get(handles.popupmenu_ovlap,'String');
param.ovlap=(str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100);

%param.morph=[];
param.morph.background=[];
param.filter=[];
param.merge=[];
param.median=[];

%Load everything
%param=load_image_parameters(param,chc_list{end},get(hObject, 'Label'));

yess=1;
while yess==1
    Ichc=menu('Select a group of image processing parameters to define ...', chc_list);
    param=load_image_parameters(param,chc_list{Ichc},get(hObject, 'Label'));
    
    [features,final_image]=extract_image_features(handles.x, cbegin,param,2);
    yess=menu('Redo?','Yes','No');
    
end

end %image_processor_Callback

%% Image processing files
function [out_param]=load_image_parameters(input_param,batch_chc,label_str)

chc_list={'Defaults','Equalization','Filtering','Ridge Extraction','Morphological Processing','All'};

%Check if input parm fields present, otherwise assign
out_param=input_param;

if nargin==1
    out_param=chc_list;
    return
end

do_all=strcmp(batch_chc,chc_list{end})||strcmp(batch_chc,'Defaults');
load_defaults=strcmp(batch_chc,'Defaults');

%%%%%%%%%%%%%%%%%%%%
%%%%Equalization%%%%
%%%%%%%%%%%%%%%%%%%%

if strcmp(chc_list{2},batch_chc)||do_all
    
    Ip=0;
    
    Batch_type	=	[label_str ' Equalization parameters'];
    Batch_vars.MinFreq='0'; Ip=Ip+1;Batch_desc{Ip}	=	'Minimum Frequency (Hz) to crop initial image';
    Batch_vars.MaxFreq='500'; Ip=Ip+1;Batch_desc{Ip}	=	'Maximum Frequency (Hz) to crop initial image';
    Batch_vars.eq_time	=	'1';            Ip=Ip+1;Batch_desc{Ip}	=	'Time (sec) used to create background  noise estimate. Used if ''equalization'' is Inf';
    Batch_vars.equalization = '[]';        Ip=Ip+1;Batch_desc{Ip}	=	'Precomputed noise equalization spectrum; first column: Frequecies in Hz; Secon column: PSD (dB) value of noise';
    
    out_param.morph=merge_image_params(out_param.morph,Batch_vars,Batch_desc,Batch_type,load_defaults);
    if isempty(out_param.morph.equalization)
        out_param.morph= rmfield(out_param.morph,'equalization');
    end
    Batch_type=[];Batch_vars=[];Batch_desc=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Image Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(chc_list{3},batch_chc)||do_all
    
    Ip=0;
    
    Batch_type	=	[label_str ' Median image filtering parameters'];
    
    Batch_vars.on='0';Ip=Ip+1;Batch_desc{Ip}	=	'0: none, 1: Conduct median filtering';
    Batch_vars.size='[0.2 20]';  Ip=Ip+1;Batch_desc{Ip}	=	'Median filter size in units of [sec Hz] for spectrogram';%
    out_param.median=merge_image_params(out_param.median,Batch_vars,Batch_desc,Batch_type,load_defaults);
    Batch_type=[];Batch_vars=[];Batch_desc=[];
    
    %parameters for gaussian filtering..
    Ip=0;
    
    Batch_type	=	[label_str ' Gaussian image filtering parameters'];
    
    Batch_vars.on='0';  Ip=Ip+1;Batch_desc{Ip}	=	'0: none, 1:asymetric gaussian, 2: symmetric gaussian';%
    Batch_vars.size='[0.2 0.2]';  Ip=Ip+1;Batch_desc{Ip}	=	'Gaussian filter size in units of [sec Hz] for spectrogram';%
    Batch_vars.sigma='0.5';  Ip=Ip+1;Batch_desc{Ip}	=	'Sigma in units of pixel size';%
    out_param.filter=merge_image_params(out_param.filter,Batch_vars,Batch_desc,Batch_type,load_defaults);
    Batch_type=[];Batch_vars=[];Batch_desc=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ridge Extraction
%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(chc_list{4},batch_chc)||do_all
    
    Ip=0;
    Batch_type	=	[label_str ' Ridge extraction parameters'];
    
    Batch_vars.threshold_chc='local_peaks'; Ip=Ip+1;Batch_desc{Ip}	=	'Algorithm for picking local ridge maximum (local_peaks[default],otsu,reconstruction)';
    Batch_vars.SNRmin	=	'10';           Ip=Ip+1;Batch_desc{Ip}	=	'Min dB level, for thresholding background noise';
    
    Batch_vars.local_bandwidth='[0 500]'; Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum value of local ridge bandwidth (Hz)';
    Batch_vars.time_band_product='[1 Inf]'; Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum value of time bandwidth product (s-Hz)';
    
    Batch_vars.dynamic_range	=	'3.11';	Ip=Ip+1;Batch_desc{Ip}	=	'%dB below ridge maximum a pixel is allowed to have';
    Batch_vars.dynamic_range2	=	'7';	Ip=Ip+1;Batch_desc{Ip}	=	'dB below maximum a pixel is allowed to have, horizontal (regional) maximum';
    
    out_param.morph=merge_image_params(out_param.morph,Batch_vars,Batch_desc,Batch_type,load_defaults);
    
    out_param.morph.local_bandwidth=parse_min_max(out_param.morph.local_bandwidth);
    out_param.morph.time_band_product=parse_min_max(out_param.morph.time_band_product);
    
    
    Batch_type=[];Batch_vars=[];Batch_desc=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Morphological Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(chc_list{5},batch_chc)||do_all
    Ip=0;
    Batch_type	=	[label_str ' Morphological processing parameters'];
    
    Batch_vars.gap_f='11';  Ip=Ip+1;Batch_desc{Ip}	=	'Ridge trace: vertical dilation object size in Hz';
    Batch_vars.gap_t='0.04';  Ip=Ip+1;Batch_desc{Ip}	=	'Ridge trace: horizontal dilation object size in sec';
    
    
    Batch_vars.duration='[0.15 5]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum duration of segment allowed (sec)';
    Batch_vars.robust_fmin='[0 500]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum  of robust minimum frequency (Hz)';
    Batch_vars.robust_fmax='[0 500]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum  of robust maximum frequency (Hz)';
    Batch_vars.robust_bandwidth='[0 500]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum  of robust bandwidth (Hz)';
    Batch_vars.Orientation='[-90 90]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum  of segment orientation (deg)';
    Batch_vars.Centroid_freq='[0 500]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum  of centroid frequency (Hz)';
    Batch_vars.Eccentricity='[0 1]';Ip=Ip+1;Batch_desc{Ip}	=	'minimum and maximum  of segment eccentricity ';
    
    out_param.morph=merge_image_params(out_param.morph,Batch_vars,Batch_desc,Batch_type,load_defaults);
    Batch_type=[];Batch_vars=[];Batch_desc=[];
    
    out_param.morph.duration=parse_min_max(out_param.morph.duration);
    out_param.morph.robust_fmin=parse_min_max(out_param.morph.robust_fmin);
    out_param.morph.robust_fmax=parse_min_max(out_param.morph.robust_fmax);
    out_param.morph.robust_bandwidth=parse_min_max(out_param.morph.robust_bandwidth);
    out_param.morph.Orientation=parse_min_max(out_param.morph.Orientation);
    out_param.morph.Centroid_freq=parse_min_max(out_param.morph.Centroid_freq);
    out_param.morph.Eccentricity=parse_min_max(out_param.morph.Eccentricity);
    
    Ip=0;
    Batch_type	=	[label_str ' Background Morphological processing parameters'];
    
    Batch_vars.on='1';  Ip=Ip+1;Batch_desc{Ip}	=	'1: Permit background contour processing; 0: ridge tracing only';
    Batch_vars.gap_f='20';Ip=Ip+1;Batch_desc{Ip}	=	'Background trace: vertical dilation object size in Hz';
    Batch_vars.gap_t='0.1';Ip=Ip+1;Batch_desc{Ip}	=	'Background trace: horizontal dilation object size in sec';
    
    out_param.morph.background=merge_image_params(out_param.morph.background,Batch_vars,Batch_desc,Batch_type,load_defaults);
    Batch_type=[];Batch_vars=[];Batch_desc=[];
    
    Ip=0;
    Batch_type	=	[label_str ' Segment merging processing parameters'];
    
    Batch_vars.ovlap='0.25';Ip=Ip+1;Batch_desc{Ip}	=	'minimum overlap between segmentes for merger (0-1.0)';
    Batch_vars.max_frequency_separation='50';Ip=Ip+1;Batch_desc{Ip}	=	'maximum frequency sepearation allowed for parallel merger (Hz)';
    Batch_vars.gap_t='0.1';Ip=Ip+1;Batch_desc{Ip}	=	'maximum time separation between segments allowed for sequential merger (sec)';
    Batch_vars.gap_f='20'; Ip=Ip+1;Batch_desc{Ip}	=	'maximum frequency separation between segments allowed for sequential merger (Hz)';
    
    out_param.merge=merge_image_params(out_param.merge,Batch_vars,Batch_desc,Batch_type,load_defaults);
    Batch_type=[];Batch_vars=[];Batch_desc=[];
    
    
end
%%%%%%%%%

out_param.morph.want_contour=[];

    function outt_param=merge_image_params(in_param,Batch_vars,Batch_desc,Batch_type,load_default)
        outt_param=in_param;
        
        if nargin<5
            load_default=0;
        end
        
        if load_default==0
            [Batch_vars, Batch_names]	=	input_batchparams(Batch_vars, Batch_desc, Batch_type,1);
            for I=1:length(Batch_desc)
                outt_param.(Batch_names{I})=Batch_vars.(Batch_names{I});
            end
        else
            Batch_names=fieldnames(Batch_vars);
            
            for I=1:length(Batch_names)
                try
                    outt_param.(Batch_names{I})=eval(Batch_vars.(Batch_names{I}));
                catch
                    outt_param.(Batch_names{I})=(Batch_vars.(Batch_names{I}));
                    
                end
            end
            
        end
        
    end

    function x=parse_min_max(x)
        warning('off');
        if ~isnumeric(x)
            temp=eval(x);
        else
            temp=x;
        end
        %temp=eval(x);
        x.min=temp(1);
        x.max=temp(2);
        warning('on');
    end
end %load_image_parameters

% --- Executes on button press in togglebutton_ChannelBeam.
function togglebutton_ChannelBeam_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_ChannelBeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_ChannelBeam

%First, check that multichannel data exist
if ~isfield(handles,'file_flags')||~handles.file_flags.multichannel
    set(hObject,'Value',0)
    set(hObject,'String','Channel')
    return
end

onn=get(hObject,'Value');
if onn
    
    set(hObject,'Value',1)
    set(hObject,'String','Angle')
    
else
    set(hObject,'Value',0)
    set(hObject,'String','Channel')
end
guidata(hObject, handles);
end


% --------------------------------------------------------------------
function batch_spectrogram_Callback(hObject, eventdata, handles)
% hObject    handle to batch_spectrogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Generate a bunch of spectrograms using current parameters.
mydir=pwd;
batch_dir = uigetdir(handles.mydir, 'Select a batch folder');
cd(batch_dir);
fnames=dir(['*.' handles.myext]);
for I=1:length(fnames)
    files{I}=fnames(I).name;
end
[SELECTION,OK] = listdlg('ListString',files);
files=files(SELECTION);

file_org=handles.myfile;

%	Stop audio playback if it's running
if	~isempty(handles.audioplayer) && isplaying(handles.audioplayer)
    stop(handles.audioplayer);
end

for I=1:length(files)
    handles.myfile=files{I};
    
    %yes_wav=get(handles.togglebutton_getwav,'value');
    tlenn=datenum(0,0,0,0,0,handles.tlen);
    Nshots=(handles.tdate_max-handles.tdate_min)/tlenn;
    handles.tdate_start=handles.tdate_min;
    [handles,~]		=	set_slider_controls(handles, handles.filetype);
    
    while ((handles.tdate_max-handles.tdate_start)>tlenn)
        set(handles.edit_datestr,'String',datestr(handles.tdate_start));
        
        edit_datestr_Callback(handles.edit_datestr, [], handles);
        
        set(handles.text_filename,'String',fullfile(handles.mydir, handles.myfile));
        
        %	Call actual function to produce figure
        handles		=	load_and_display_spectrogram(handles);
        pushbutton_print_Callback([], [], handles);
        %pushbutton_update_Callback(hObject, eventdata, handles);
        disp('done')
        handles.tdate_start=handles.tdate_start+tlenn;
        
        pause(1)
    end
    
end


end

%% Annotation Utilties
% --------------------------------------------------------------------
function Annotation_Utilities_Callback(hObject, eventdata, handles)
% hObject    handle to Annotation_Utilities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end



% --------------------------------------------------------------------
function airgun_detector_Callback(hObject, eventdata, handles)
% hObject    handle to airgun_detector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Convert bowhead whale files into annotation files....

%Load list of file names to be converted, along with criteria for filtering
% the results into annotation (e.g. geographic restrictions...)
[list_names,filter_params]=load_airgun_detector_params(handles.outputdir);
if isempty(list_names)
    return
end

for I=1:length(list_names)
    success_flag=convert_automated_airgun_into_annotations(list_names{I},filter_params);
    if success_flag==0
        uiwait(msgbox(sprintf('%s failed to process',list_names{I}),'Failed!','modal'));
        return
    else
        disp(sprintf('%s processed',list_names{I}));
        
    end
    
end

uiwait(msgbox(sprintf('Conversion finished')));
end

% --------------------------------------------------------------------
function bowhead_detector_Callback(hObject, eventdata, handles)
% hObject    handle to bowhead_detector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Convert bowhead whale files into annotation files....

%Load list of file names to be converted, along with criteria for filtering
% the results into annotation (e.g. geographic restrictions...)
[list_names,filter_params,station_position]=load_bowhead_detector_params(handles.outputdir);
if isempty(list_names)
    return
end

for I=1:length(list_names)
    sprintf('%s starting ...',list_names{I});
    success_flag=convert_automated_bowhead_into_annotations(list_names{I},filter_params,station_position);
    if success_flag==0
        uiwait(msgbox(sprintf('%s failed to process',list_names{I}),'Failed!','modal'));
        
    else
        disp(sprintf('%s procesed',list_names{I}));
        
    end
    
end

uiwait(msgbox(sprintf('Conversion finished')));


end



%% Time warping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%TIMEWARP functions below %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_timewarp.
function pushbutton_timewarp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_timewarp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;
Ichan	=	str2double(get(handles.edit_chan,'String'));  %Hardwire first channel

if ~strcmp(handles.filetype,'MDAT')
    [x0,~,Fs,~,~,hdr]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,Ichan,handles);
    Ichan=1;Nchan=1;
else
    Nchan=1:16;
    %Ichan=16-Ichan;  %Fix this; why upside down?
    [x0,~,Fs,~,~,hdr]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,Nchan,handles);
    Nchan=1:size(x0,1);
end

RR=10;  %Decimation factor
for I=Nchan
    x(I,:)=decimate(x0(I,:),RR,'FIR');
end
Fs=Fs/RR;

contents=get(handles.popupmenu_Nfft,'String');
Nfft=(contents{get(handles.popupmenu_Nfft,'Value')});


prompt={'Range guess(m)','water speed (m/s)','Nfft','N_window','N_window_w', ...
    'Time to clip off front of deconvolved signal (to ensure time 0 is first mode arrival)', 'Number of FM contour points:'};
def={'15000','1490',Nfft,'2*floor(0.5*31*Fs/200)+1','2*floor(0.5*301*Fs/200)+1','0.8','3'};
dlgTitle	=	sprintf('Warping parameters');
lineNo		=	ones(size(prompt));
answer		=	inputdlg(prompt,dlgTitle,lineNo,def);

if	isempty(answer)
    return;
end

r_guess=eval(answer{1});
c1=eval(answer{2});
Nfft=eval(answer{3});
N_window=round(eval(answer{4}));
N_window_w=round(eval(answer{5}));
n_limit=eval(answer{6});
Npts=eval(answer{7});

%Select FM-sweep guess...
fprintf('Select %i points from lowest-order mode (best guess of original FM structure:\n)',Npts);
tmp=ginput(Npts);
fstart=(tmp(1:(end-1),2)*1000);
fend=(tmp(2:end,2)*1000);
tsweep=abs(tmp(2:end,1)-tmp(1:(end-1),1));

minfreq=min([fstart(1) fend(end)]); maxfreq=max([fstart(1) fend(end)]);
frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
[N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.1],Fs);
B = firpm(N,Fo,Ao,W);
y=filter(B,1,x(Ichan,:)-mean(x(Ichan,:)));
%[s_w,mode_t,spectro_w,time_w, freq_w, filter_mask, Fe_w]=
modes=warp_transform(1,hilbert(y.'),r_guess,c1,Fs, Nfft, N_window,N_window_w,n_limit,tsweep,fstart,fend,[],hdr.geom.rd(Ichan));
%close all
MM=size(modes.s_filt);
if MM(1)>0
    
    %%mode_stack is unwarped time series that should be mode-dominated
    mode_stack=zeros(MM(1),MM(2),length(Nchan));
    %sw_stack=zeros(length(Nchan),length(s_w));
    for I=Nchan
        fprintf('%i of %i channels\n',I,length(Nchan));
        y=filter(B,1,x(I,:)-mean(x(I,:)));
        tmp=warp_transform(0,hilbert(y.'),r_guess,c1,Fs, Nfft, N_window,N_window_w,n_limit,tsweep,fstart,fend,modes.spectro_mask,hdr.geom.rd(I));
        mode_stack(:,:,I)=tmp.s_filt;
        
    end
end
close(10:11)
clear tmp
save temp


%depth_shift=-2;  %How much to shift my estimated depth by
%freq_want=[84 104];
%freq_want=[120:5:145];
prompt={'Frequencies at which to extract modes (Hz):'};
def={'[84 104]'};
dlgTitle	=	sprintf('Frequencies to extract modes');
lineNo		=	ones(size(prompt));
answer		=	inputdlg(prompt,dlgTitle,lineNo,def);

if	isempty(answer)
    return;
end

freq_want=eval(answer{1});

Nfft_filt=2^nextpow2(MM(1));
%Nfft_filt=512;
Uwarp=fft(mode_stack,Nfft_filt);  %Spectrum of each mode arrival
Fwarp=linspace(0,Fs,Nfft_filt);
Igood=find(Fwarp>=(minfreq)&(Fwarp<=(maxfreq)));
Uwarp=Uwarp(Igood,:,:);
%Sum energy in each modal arrival. Assign sign later.
Umode_broadband=squeeze(sum(mode_stack.^2,1))';

Nc=length(Nchan);
Nf=length(Igood);
%Ifreq_plot=[1 Nf/4 Nf/2 3*Nf/4 Nf];

%Option: pick frequency with most power.
%Fpower=squeeze(sum(sum(shiftdim(abs(Uwarp).^2,1))));
%Ifreq_plot=[1 Nf/2  Nf];
%Ifreq_plot=round(([ 84 98 105]-minfreq)*Nfft_filt/Fs);
Ifreq_plot=round((freq_want-minfreq)*Nfft_filt/Fs);

freq_plot=Fwarp(Igood(Ifreq_plot));
Nplot=length(freq_plot);

%Cyle through each mode and extract mode shape from various frequencies,
% as well as broadband time series
for Im=1:MM(2)
    
    %normalize both individual-frequency based and broadband modes.
    Umode=squeeze(Uwarp(:,Im,:)).';
    Umode_broadband(:,Im)=Umode_broadband(:,Im)/norm(Umode_broadband(:,Im));
    Umode_broadband(:,Im)=sqrt(Umode_broadband(:,Im));
    
    for II=1:Nf
        Umode(:,II)=Umode(:,II)/norm(Umode(:,II));
    end
    
    %Extact phase along array
    Uang_uncorrected=angle(Umode./(ones(Nc,1)*Umode(end,:)));
    if Im==1
        Umode1=Umode;  %We assume that mode 1 phase arises from tilt
    end
    
    tmp=(Umode.*conj(Umode1));  %Remove tilt or phase offsets
    Uang=angle(tmp./(ones(Nc,1)*tmp(end,:)));  %Reference phase to bottom element
    
    %Determine sign of mode: observe when corrected phase
    %   (relative to bottom mode) exceeds pi/2
    sgn_chnge=-(abs(Uang(:,Ifreq_plot))>pi/2);
    sgn_chnge=2*sgn_chnge+1;  %Convert 0 and -1 to 1 and -1;
    if Nplot>1
        sgn_chnge=median(sgn_chnge')';
    end
    Umode_broadband(:,Im)=sgn_chnge.*Umode_broadband(:,Im);
    Umode_shape(:,:,Im)=(sgn_chnge*ones(1,Nplot)).*abs(Umode(:,Ifreq_plot));
    
    %end
    %tmp=sum(Umode.^2);
    %[junk,Imax]=max(tmp);
    
    figure(20);
    subplot(MM(2),3,3*Im-2);
    plot(abs(Umode(:,Ifreq_plot)),hdr.geom.rd,'o-');grid;axis('ij')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('%s: abs value of mode %i',datestr(tdate_start),Im));
    if Im==1
        ylabel('depth (m)');
    end
    
    subplot(MM(2),3,3*Im-1);
    plot((Uang_uncorrected(:,Ifreq_plot)),hdr.geom.rd,'o-');grid;axis('ij')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('uncorrected phase of mode %i',Im));
    if Im==1
        ylabel('depth (m)');
    end
    
    subplot(MM(2),3,3*Im);
    plot((Uang(:,Ifreq_plot)),hdr.geom.rd,'o-');grid;axis('ij')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('tilt-corrected phase of mode %i',Im));
    xlabel('Phase (rad)');
    %     plot(Umode_fin(:,Ifreq_plot),hdr.geom.rd,'o-');grid;axis('ij')
    %     title(sprintf('Recovered mode %i',Im));
    if Im==1
        ylabel('depth (m)');
    end
    
    figure(19)
    subplot(1,MM(2),Im)
    plot(Umode_broadband(:,Im),hdr.geom.rd,'k',squeeze(Umode_shape(:,:,Im)),hdr.geom.rd,'g');grid on
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('Normalized mode %i',Im));
    ylabel('depth(m)');
    title(sprintf('Inverted mode %i',Im));
    axis('ij')
    legend('broadband');
    
    
end

%%Plot recovered tilt
U1_phase=angle(Umode1./(ones(Nc,1)*Umode1(end,:)));
U1_tilt=unwrap(U1_phase)./(ones(Nc,1)*2*pi*Fwarp(Igood)/1420);

%Simple broadband calculations
xt1=squeeze(mode_stack(:,1,:));
RR=10;
maxlag=round(50*RR*Fs/1490);

for I=1:(Nc-1)
    y1=resample(xt1(:,I),RR,1);
    y2=resample(xt1(:,I+1),RR,1);
    [cc,lags]=xcov(y1,y2,maxlag);
    [junk,Imax]=max(abs(cc));
    tilt_broadband(I+1)=1490*lags(Imax)/(RR*Fs);
    
end
tilt_broadband=cumsum(tilt_broadband);
tilt_broadband=tilt_broadband-tilt_broadband(end);  %Make bottom phone the origin

figure(18)
subplot(1,2,1)
plot(U1_phase(:,Ifreq_plot),hdr.geom.rd);axis('ij');grid
xlabel('rad');ylabel('depth (m)');
title('Array tilt in terms of mode 1 phase');

subplot(1,2,2);
plot(U1_tilt(:,Ifreq_plot),hdr.geom.rd);axis('ij');grid
hold on
plot(tilt_broadband,hdr.geom.rd,'ko')
title('Array tilt in terms of mode 1 offset');
xlabel('m offset');ylabel('depth (m)');
grid on
legend(strvcat(num2str(freq_plot','%2.1f'),'broadband'));

figure(20)
legend(num2str(freq_plot','%2.1f'));


for II=18:20
    figure(II)
    orient landscape
    print('-djpeg','-r300',sprintf('ModeTiltWarpInversion_%i',II));
end
end


%warp_transform.m
function modes= warp_transform(plot_chc,s_ok,r_guess,c1,Fe, Nfft, N_window,N_window_w,t_trim, ...
    t_sweep,fstart,fend,spectro_mask_in,rd) % You need a long window to see the warped modes)
% The time-frequency representations that are used in this script are computed using a free
% toolbox that you can find here : http://tftb.nongnu.org/
% Be sure to add the tftb-01 folder to your matlab file if you want this to
% work !
% Inputs:
%  plot_chc:  If 1, make plots
%  s_ok: input signal, single channel
%  r_guess: initial guess of signal range, in meters
%  c1: initial estimate of water sound speed, m/s
%  Fe: sampling rate, Hz
%  Nfft: FFT size
%  N_window:  number of FFT window points for original time series
%  N_window_w: number of FFT window points for warped time series
%  t_trim: number of sec to remove from time series after source deconvolution estimate
%  fstart: vector of start frequencies of piecemeal linear FM sweeps, Hz
%  fend: vector of end frequencies of pieceal FM sweeps, Hz
%  t_sweep: duration of these sweeps, sec
%  spectro_mask_in:  masking window for filtering, externally provided from processing a different channel.
% pt: propagated signal in a Pekeris waveguide with parameters
%       c1, c2, rho1, rho2: sound speed / density
%       D: depth
%       r: range
%       zs, zr: source/receiver depth
% Fe: sampling frequency
%  Output:
%       modes:
%           s_w:  warped deconvolved input signal
%           s_filt: IFFT of mode-selections
%           spectro_w: spectrogram of s_w, evaluated for every time sample
%           spectro_mask: blocked out region of interest for inverse filtering.  3D arrays
%           time_w, freq_w: time and frequency scale for spectro_w
%           Fe_w:  Effective sampling rate in Hz for s_w


%N_window=31; % you need a short window to see the modes
%N_window_w=301; % You need a long window to see the warped modes
%t_trim=400;% Don't forget to change that if you consider other data

if size(s_ok,2)>1
    s_ok=s_ok.';
end

if ~exist('rd')
    rd=[];
end
N_ok=length(s_ok);
Nfft_ok=2^nextpow2(N_ok);
% Looking at the received signal, you can guess that the source is a piecewise linear FM
% downsweep. This is wrong but it shows you that you don't need to know
% the source exactly
%source_est=fmlin(160,0.5,0.01);

iflaw=[];
for I=1:length(fstart)
    iflaw_seg=linspace(fstart(I)/Fe,fend(I)/Fe,round(Fe*t_sweep(I)));
    iflaw=[iflaw; iflaw_seg(2:end)'];
end

%source_est=fmlin(round(t_sweep*Fe),fstart/Fe,fend/Fe);
source_est=fmodany(iflaw);

source_est_f=fft(source_est,Nfft_ok);
phi=angle(source_est_f);

% Phase correction for estimated source signal

%Trim front of signal to ensure first mode arrival is time zero, if necessary...
%s_ok=s_ok(round(Fe*t_trim):end);

sig_prop_f=fft(s_ok,Nfft_ok);
sig_rec_f=sig_prop_f.*exp(-1i*phi); % note that only source phase is deconvoluted
sig_rec_t=ifft(sig_rec_f,Nfft_ok); % Signal in time domain after source deconvolution


% STFT computation of original data
tfr=tfrstft(s_ok,1:N_ok,Nfft,hamming(N_window));
spectro=abs(tfr).^2;
% Time and frequency axis of the original signal
time=0:1/Fe:(N_ok-1)/Fe;
freq=0:Fe/Nfft:Fe-Fe/Nfft;

%Trim front of signal to ensure first mode arrival is time zero, if necessary...
sig_rec_t=sig_rec_t(round(Fe*t_trim):end);

tfr_sig_rec=tfrstft(sig_rec_t,1:N_ok,Nfft,hamming(N_window));

% Figure

if plot_chc==1
    figure(10)
    subplot(2,2,1);
    imagesc(time, freq, 10*log10(spectro));caxis([80 130])
    set(gca,'fontweight','bold','fontsize',14);
    axis xy
    ylim([0 Fe/2])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title(sprintf('Original signal at %2.0f m depth',rd))
    colorbar
    
    subplot(2,2,2)
    imagesc(time, freq, 20*log10(abs(tfr_sig_rec)))
    set(gca,'fontweight','bold','fontsize',14);
    
    axis xy
    ylim([0 Fe/2])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Signal after source deconvolution')
    colorbar
    orient landscape
    caxis([80 130])
    %print('-djpeg',sprintf('DeconvSpectrogram_depth%2.0fm',rd));
    
end

for Iguess=1:length(r_guess)
    time=r_guess(Iguess)/c1:1/Fe:r_guess(Iguess)/c1+(N_ok-1)/Fe;
    
    % The warped signal will be s_w
    [s_w, Fe_w]=warp_temp(sig_rec_t,Fe,r_guess(Iguess),c1);
    M=length(s_w);
    
    tfr_w=tfrstft(s_w,1:M,Nfft,hamming(N_window_w));
    spectro_w=abs(tfr_w).^2;
    % Time and frequency axis of the warped signal
    time_w=0:1/Fe_w:(M-1)/Fe_w;
    freq_w=0:Fe_w/Nfft:Fe_w-Fe_w/Nfft; % this is actually not exact, but precise enough for a rough example
    % Figure
    
    if plot_chc==1
        subplot(2,2,3)
        imagesc(time_w, freq_w, spectro_w)
        set(gca,'fontweight','bold','fontsize',14);
        
        axis xy
        ylim([0 50])
        xlabel('Warped time (sec)')
        ylabel('Corresponding warped frequency (Hz)')
        title(sprintf('Warped signal square amp. assuming range %6.2f',r_guess(Iguess)))
        caxis([1e10 5e11]);colorbar
        
        subplot(2,2,4)
        imagesc(time_w, freq_w, 10*log10(abs(spectro_w)))
        set(gca,'fontweight','bold','fontsize',14);
        
        axis xy
        ylim([0 50])
        xlabel('Warped time (sec)')
        ylabel('Corresponding warped frequency (Hz)')
        title('dB value');
        caxis([80 120]);colorbar;
        
        orient landscape
        print('-djpeg','-r300',sprintf('WarpedSpectrogram_range%2.0fkm_depth%2.0fm',r_guess(Iguess)/1000,rd));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%

% To make it easier, filtering will be done by hand using the roipoly tool.
% See Matlab help for more information

% We work on a smaller TFR to make it easier

if ~exist('spectro_mask_in')||isempty(spectro_mask_in)
    %figure(5)
    spectro_mask_in=[];
    %RTF=abs(tfr_w).^2;%RTF=RTF/max(max(RTF));
    figure(10);
    %imagesc(time_w, freq_w, 10*log10(RTF));axis xy;caxis([90 120]); colorbar; ylim([0 50])
    BW_flag=0;
    Nroi=input('Enter number of modes to select, or return to skip:');
    modes.s_filt=[];
    modes.spectro_mask=[];
    if isempty(Nroi)
        return
    end
else
    BW_flag=1;
    spectro_mask=spectro_mask_in;
    Nroi=size(spectro_mask,3);
    fprintf('Input spectro_mask detected with %i modes.\n',Nroi);
end

s_filt=zeros(length(s_ok),Nroi);
for J=1:Nroi
    
    if BW_flag==0
        figure(10);
        BW0 = roipoly ;
        spectro_mask(:,:,J)=BW0;
    else
        BW0=squeeze(spectro_mask(:,:,J));
    end
    % To create a time-frequency mask, click several times on the spectrogram
    % to define the region you want to filter. The last click must be a
    % double-click to validate the mask creation
    %masque=zeros(size(tfr_w));
    %masque(1:n_limit-1,:)=spectro_mask;
    
    % Note that the mask is applied on the STFT (not on the spectrogram)
    mode_rtf_warp=BW0.*tfr_w;
    
    %%The trick below works if the FFT spectrogram is computed for every sample...
    [mode_temp_warp]=real(sum(mode_rtf_warp,1));
    s_filt(:,J)=iwarp_temp(mode_temp_warp,Fe_w,r_guess,c1,Fe,N_ok)';
    
    if plot_chc==1
        figure(11);
        
        subplot(Nroi,2,2*J-1)
        %RTF=abs(mode_rtf_warp/max(max(abs(tfr_w))));RTF=RTF./max(max(RTF));
        RTF=20*log10(abs(mode_rtf_warp));
        imagesc(time_w, freq_w, RTF);axis xy;ylim([0 50]);
        caxis([90 120]);colorbar
        
        subplot(Nroi,2,2*J);
        [~,FF,TT,BB] = spectrogram(s_filt(:,J),hanning(32),round(0.5*32),32,Fe);
        imagesc(TT,FF,10*log10(BB));%
        axis('xy');
        caxis([90 120]);colorbar
        title(sprintf('Recovered mode %i at depth: %6.2f m',J,rd));
        if J==Nroi
            orient tall
            print('-djpeg','-r300',sprintf('warp_tranform_%2.0fm',round(rd)))
        end
    end
end

modes.s_w=s_w;
modes.s_filt=s_filt;
modes.spectro_w=spectro_w;
modes.time_w=time_w;
modes.freq_w=freq_w;
modes.spectro_mask=spectro_mask;
modes.Fe_w=Fe_w;

end  %warp_transform


function [ t_w ] = warp_t(t,r,c)

t_w=sqrt(t.^2+r^2/c^2);

end

function [ t ] = iwarp_t(t_w,r,c)

t=sqrt(t_w.^2-r^2/c^2);
end


function [tfr,t,f] = tfrstft(x,t,N,h,trace)
%TFRSTFT Short time Fourier transform.
%	[TFR,T,F]=TFRSTFT(X,T,N,H,TRACE) computes the short-time Fourier
%	transform of a discrete-time signal X.  The overlap is very large:
%        A FFT is computed for every time sample, with one sample change.
%
%	X     : signal (analytic; only has positive frequencies)
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	H     : frequency smoothing window, H being normalized so as to
%	        be  of unit energy.      (default : Hamming(N/4)).
%           Number of time points per FFT column
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR   : time-frequency decomposition (complex values). The
%	        frequency axis is graduated from -0.5 to 0.5.
%	F     : vector of normalized frequencies.
%
%	Example :
%	 sig=[fmconst(128,0.2);fmconst(128,0.4)]; tfr=tfrstft(sig);
%	 subplot(211); imagesc(abs(tfr));
%	 subplot(212); imagesc(angle(tfr));
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-August 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
%	------------------- CONFIDENTIAL PROGRAM --------------------
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to
%	f.auger@ieee.org

[xrow,xcol] = size(x);
if (nargin < 1),
    error('At least 1 parameter is required');
elseif (nargin <= 2),
    N=xrow;
end;

hlength=floor(N/4);
hlength=hlength+1-rem(hlength,2);

if (nargin == 1),
    t=1:xrow; h = tftb_window(hlength); trace=0;
elseif (nargin == 2) || (nargin == 3),
    h = tftb_window(hlength); trace = 0;
elseif (nargin == 4),
    trace = 0;
end;

if (N<0),
    error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol~=1),
    error('X must have one column');
elseif (trow~=1),
    error('T must only have one row');
elseif (2^nextpow2(N)~=N),
    fprintf('For a faster computation, N should be a power of two\n');
end;

[hrow,hcol]=size(h); Lh=(hrow-1)/2;
if (hcol~=1)||(rem(hrow,2)==0)
    error('H must be a smoothing window with odd length');
end;

h=h/norm(h);

tfr= zeros (N,tcol) ;
if trace, disp('Short-time Fourier transform'); end;
for icol=1:tcol
    ti= t(icol);
    tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
    %if trace, disprog(icol,tcol,10); end;
    tfr(indices,icol)=x(ti+tau,1).*conj(h(Lh+1+tau));
end;
tfr=fft(tfr);
if trace, fprintf('\n'); end;

if (nargout==0),
    %tfrqview(abs(tfr).^2,x,t,'tfrstft',h);
elseif (nargout==3),
    if rem(N,2)==0,
        f=[0:N/2-1 -N/2:-1]'/N;
    else
        f=[0:(N-1)/2 -(N-1)/2:-1]'/N;
    end;
end




end

function [s_w, Fe_w]=warp_temp(s,Fe,r,c)

% Time warping function using isovelocity waveguide as a warping model.
% Inputs:
% s : signal that will be warped
% Fe : sampling frequency of s
% r : warping parameter (source/receiver range)
% c : warping parameter (water sound speed)
% Outputs :
% s_w : warped signal
% Fe_w : sampling frequency of s_w (i.e. sampling frequency
%        of the warped-time domain)

% NB : warping actually has a single parameter : r/c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%% Important %%%%%%%%%%%%%%%%%%%%%%%%
% Signal s should be defined such that its first sample
% corresponds to the first modal arrival. In an isovelocity or Pekeris
% waveguide, if an impulsive source explodes at time ts, then first sample
% of s should be take at time ts+r/c. Warping is way more sensitive to
% that setting than to the actual parameters r and c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Step 1: preliminary computations
s=real(s);
N=length(s);
dt=1/Fe;

tmax=(N-1)/Fe+r/c;
tmin=r/c;

%% Step 2: new time step
dt_w=iwarp_t(tmax,r,c)-iwarp_t(tmax-dt,r,c);

%% Step 3: new sampling frequency
Fe_w=2/dt_w;

%% Step 4: new number of points
t_w_max=iwarp_t(tmax,r,c);
t_w_min=0;
M= ceil(t_w_max*Fe_w);


%% Step 5: warped signal computation

s_w=zeros(M,1);

for m=1:M
    
    t_w=m/Fe_w;
    % energy conservation
    toto=sqrt(Fe/Fe_w)*sqrt(t_w/sqrt(t_w^2+r^2/c^2));
    
    % linear interpolation of signal from original time domain
    n=(warp_t(t_w,r,c)-r/c)*Fe;
    n1=floor(n);
    n2=n1+1;
    dn=n-n1;
    
    if n2<=N && n1>0
        a=(s(n1)-s(n2))/(n1-n2);
        val=s(n1)+a*dn;
        s_w(m)=toto*val;
    end
    
end


end

function [s_r]=iwarp_temp(s_w,Fe_w,r,c,Fe,N)

% Inverse (time) warping function using isovelocity waveguide as a warping model.
% Inputs :
% s_w : (warped) signal that will be unwarped
% Fe_w : sampling frequency of s_w
% r : warping parameter (source/receiver distance)
% c : warping parameter (source/receiver distance)
% Fe : sampling frequency of the new unwarped signal s
% N : number of points of the new unwarped signal s
% Sorties :
% s : signal after inverse warping


%%

M=length(s_w);
s_r=zeros(N,1);

for n=1:N
    
    t=n/Fe+r/c;
    % Energy conservation
    toto=sqrt(Fe_w/Fe)*sqrt(t/sqrt(t^2-r^2/c^2));
    
    % Interpolation
    m=(iwarp_t(t,r,c))*Fe_w;
    m1=floor(m);
    m2=m1+1;
    dm=m-m1;
    
    if m2<=M && m1>0
        a=(s_w(m1)-s_w(m2))/(m1-m2);
        val=s_w(m1)+a*dm;
        s_r(n)=toto*val;
    end
    
end


end




