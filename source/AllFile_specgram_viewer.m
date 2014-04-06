

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

% Last Modified by GUIDE v2.5 28-Feb-2014 10:04:21

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

path('~/Desktop/Ulysses',path);
path('~/Desktop/Ulysses/deps',path);

startupinfo		=	gui_startup_information;

%Make a splash screen
try
    s = SplashScreen( 'Splashscreen', '~/Desktop/Ulysses/104703.JPG', ...
        'ProgressBar', 'on', ...
        'ProgressPosition', 5, ...
        'ProgressRatio', 0.4 );
    s.addText( 30, 50, 'Ulysses/Ribbit', 'FontSize', 30, 'Color', [0 0 0.6] )
    s.addText( 30, 80, 'v1.0, March 3, 2014', 'FontSize', 20, 'Color', [0.2 0.2 0.5] )
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
%Make GSIbearing button, CSDM, and accelerometer buttons invisible
set(handles.pushbutton_GSIbearing,'Vis','off');
set(handles.pushbutton_GSI_localization,'Vis','off');
set(handles.pushbutton_CSDM,'Vis','off');
set(handles.pushbutton_Mode,'Vis','off');
set(handles.pushbutton_tilt,'Vis','off');
set(handles.pushbutton_modalfiltering,'Vis','off');
set(handles.edit_normal_rotation,'Vis','off');
set(handles.text_normal_rotation,'Vis','off');

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
File_types	=	{'MT';'SIO';'WAV';'GSI';'ADI';'DAT';'MDAT';'PSD';'MAT'};
File_descs	=	...
 {	'MT files';...
    'WAV files';...
    'SIO files';...
    'GSI files';...
    'ADI files';...
    'DAT files';...
    'MDAT files';...
    'PSD binary files'; ...
    'Simulations'};

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
persistent Batch_vars Batch_desc

Batch_type	=	get(hObject, 'Label');
Batch_mode	=	input_batchmode(Batch_type);

if	isempty(Batch_mode)
    errordlg('Processing cancelled!');
    return;
end

if isempty(Batch_vars)
    Batch_vars.threshold	=	'5';	Batch_desc{1}	=	'dB threshold between peak and surrounding values at +/-df_search';
    Batch_vars.df_search	=	'10';	Batch_desc{2}	=	'+-/Hz to examine around each local maximum';
    Batch_vars.f_min	=	get(handles.edit_fmin,'string');	Batch_desc{3}	=	'minimum frequency to examine (kHz)';
    Batch_vars.f_max	=	get(handles.edit_fmax,'string');	Batch_desc{4}	=	'maximum frequency to examine (kHz)';
    Batch_vars.sec_avg	=	'2';	Batch_desc{5}	=	'Averaging time of spectrogram (sec)';
    
end
Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, Batch_type);


if	isempty(Batch_vars)
    uiwait(errordlg('Vessel Processing cancelled!'));
    return;
end

try
    sec_avg=str2double(Batch_vars.sec_avg);
    threshold=str2double(Batch_vars.threshold);
    df_search=str2double(Batch_vars.df_search);
    f_min=1000*str2double(Batch_vars.f_min);
    f_max=1000*str2double(Batch_vars.f_max);
    if f_min>f_max
        return
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
            PfdB_delta2=   diff(PSD_dB',2)./((F(2)-F(1)).^2);
            
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
                    df_search=str2double(Batch_vars.df_search);
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
        
        Batch_desc{1}	=	'Time to begin loading (e.g. "here" to use visible start time, "datenum(2011,1,2,0,0,0)", "file start" for current file , "select" to choose file from list, "folder start" for first file in folder)';
        Batch_vars_bulkload.end_time	=	'folder end';
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
                print(hprint,'-djpeg',save_str);
            case 'Plot percentiles'
                
                pms=get_PSD_percentile_params(fmin,fmax);
                pms.label_style=get_datetick_style([]);
                
                pms.y_label='dB re 1uPa';
                Tabs=datenum(get(handles.edit_datestr,'String'))+datenum(0,0,0,0,0,Tnew);
                if length(pms.fmin)>3
                    suppress_output=1;
                else
                    suppress_output=0;
                end
                
                %%Process PSD matrix..
                for Iff=1:length(pms.fmin)
                    Igood= (F>=pms.fmin(Iff)&F<=pms.fmax(Iff));
                    sumPSD=10*log10((F(2)-F(1))*sum(PSD(Igood,:),1));  %Now power spectral density converted to power
                    pms.title=sprintf('Spectral power between %i and %i Hz, beginning %s',pms.fmin(Iff),pms.fmax(Iff),datestr(Tabs(1)));
                    [temp,hprint]=create_percentile_distributions(Tabs, sumPSD,pms, suppress_output);
                    if Iff==1
                        p_matrix=zeros(length(pms.fmin),length(pms.percentiles),size(temp.data,2));
                        XX=temp.x;
                    end
                    p_matrix(Iff,:,:)=temp.data;
                end
                
                %%Plot contour plots if enough frequency bands
                if suppress_output
                   for Ipp=1:length(pms.percentiles)
                      figure
                      imagesc(XX,pms.fmin,squeeze(p_matrix(:,Ipp,:)));
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
        Batch_vars_bulkload.start_time	=	'folder start';
     
        Batch_desc{1}	=	'Time to begin loading (e.g. "here" to use visible start time, "datenum(2011,1,2,0,0,0)", "file start" for current file , "select" to choose file from list, "folder start" for first file in folder)';
        Batch_vars_bulkload.end_time	=	'folder end';
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
        PSD_all=[];Tabs_all=[];
        Iplot=1;
        for I=Icurrent_file:Ifinal_file
            fname=Other_FileNames(I).name;
            fprintf('Processing %s...\n',fname);
            [PSD,F,Tsec,Tabs,params]=read_Java_PSD(fname,tabs_loop_begin,Inf);  %Read an entire file into memory
            sec_avg=(params.Nfft+params.dn*params.Nsamps)/params.Fs;
        
            Istrip=find(Tabs<=tabs_end);
            switch PSDButtonName
                case 'Display and Print Figure'
                    PSD_all=[PSD_all PSD(:,Istrip)];
                    Tabs_all=[Tabs_all Tabs(Istrip)];
                    if split_windows
                        [hprint(Iplot),save_tag{Iplot}]=image_PSD(min(Tabs_all), max(Tabs_all));
                        PSD_all=[];Tabs_all=[];
                        Iplot=Iplot+1;
                    end
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
            
            %If we move into additional files, start at the beginning
            tabs_loop_begin=0;
            
        end  %Icurrent_file
        
                 
        switch PSDButtonName
            
            case 'Plot percentiles'
                %Remove excess storage
                
                Icount=Icount-length(Istrip);
                if Icount<length(Tabs_all)
                    Tabs_all=Tabs_all(1:Icount);
                    PSD_all=PSD_all(:,1:Icount);
                    
                end
                pms.y_label='dB re 1uPa';
                yes=1;
                while yes
                    close_all_figures;  %Close all open figures..
                    if Nf>3
                        %%Process PSD matrix..
                        for Iff=1:length(pms.fmin)
                            %Igood= (F>=pms.fmin(Iff)&F<=pms.fmax(Iff));
                            sumPSD=PSD_all(Iff,:);  %Now power spectral density converted to power
                            pms.title=sprintf('Spectral power between %i and %i Hz, beginning %s',pms.fmin(Iff),pms.fmax(Iff),datestr(Tabs_all(1)));
                            [temp,hprint]=create_percentile_distributions(Tabs_all, sumPSD,pms, 1);
                            if Iff==1
                                p_matrix=zeros(length(pms.fmin),length(pms.percentiles),size(temp.data,2));
                                XX=temp.x;
                            end
                            p_matrix(Iff,:,:)=temp.data;
                        end
                        
                        %%Plot contour plots if enough frequency bands
                        clear hprint
                        for Ipp=1:length(pms.percentiles)
                            hprint(Ipp)=figure;
                            imagesc(XX,pms.fmin,squeeze(p_matrix(:,Ipp,:)));
                            axis('xy');
                            %datetick('x',date_tick_chc);
                            Ibin=get(gca,'xtick');
                            Istep=floor(pms.xlabel_inc/pms.x_inc);
                            Ibin=1:Istep:length(XX);
                            
                            %Ibin=XX(Ibin);
                            set(gca,'xtick',XX(Ibin));
                            if isempty(date_tick_chc) %auto adjustment has failed..
                                date_tick_chc=get_datetick_style('auto',pms.xlabel_inc,'datenumber');
                            end
                            datetick('x',date_tick_chc,'keepticks','keeplimits')
                            %set(gca,'xticklabel',);
                            
                            
                            
                            colorbar
                            caxis(pms.y_limits);
                            title(sprintf('%ith percentile, %6.2f s average, starting %s',100*pms.percentiles(Ipp),sec_avg,datestr(Tabs(1))));
                            xlabel('Date/Time','fontweight','bold','fontsize',14);
                            ylabel('Frequency (Hz)','fontweight','bold','fontsize',14);
                            grid on;
                            
                        end
                    else  %only one or two frequeny bands
                        if isempty(date_tick_chc) %auto adjustment has failed..
                            date_tick_chc=get_datetick_style('auto',pms.xlabel_inc,'datenumber');
                        end
                        pms.label_style=date_tick_chc;
                        pms.title=sprintf('Spectral power between %i and %i Hz, averaged %6.2f sec, beginning %s', ...
                        pms.fmin(Iff),pms.fmax(Iff),sec_avg,datestr(Tabs_all(1)));
                            
                        [~,hprint]=create_percentile_distributions(Tabs_all, squeeze(PSD_all), pms,0);
                        %[output(Iff,:,:),hprint]=create_percentile_distributions(Tabs, sumPSD,pms, suppress_output);
                        %datetick('x',date_tick_chc,'keeplimits','keepticks');
                        xlabel(pms.x_label,'fontsize',14,'fontweight','bold');
                        ylabel(pms.y_label,'fontsize',14,'fontweight','bold');
                        uiwait(msgbox('Please click on screen to print out percentiles used'));
                        gtext(sprintf('Percentiles: %s',num2str(pms.percentiles)),'fontweight','bold','fontsize',18)

                        
                    end  %if Nf?3
                    yes=menu('Redo formatting? (Time formatting only)','Yes','No');
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
                    save_tag=sprintf('PSD_%s_fmin%iHz_fmax%iHz_interval%isec', ...
                        datestr(Tabs_all(1),30),(pms.fmin(Iprint)),(pms.fmax(Iprint)),tmp);
                    print(hprint(Iprint),'-djpeg',save_tag)
                    saveas(hprint(Iprint), save_tag, 'fig');
                    save(save_tag,'Tabs_all','PSD_all','pms','params','titlestr','Batch_vars_bulkload','sec_avg');
            
                end
                
            case 'Display and Print Figure'
                
                
                if ~split_windows
                    [hprint(1),save_tag{1}]=image_PSD(min(Tabs_all), max(Tabs_all));
                end
                
                yess=menu('Save PSD Data and figure?','Yes','No');
                if yess==1
                    
                    Nfigs=sort(get(0,'Child'));
                    Nfigs=Nfigs(Nfigs<150);
                    
                    for II=1:length(Nfigs)
                        
                        
                        save_str=sprintf('PSD_%s', save_tag{II});
                        save(save_str,'pms','PSD_all','Tabs_all','params','titlestr','Batch_vars_bulkload','sec_avg');
                        uiwait(msgbox([save_str ' mat and jpg file written to ' pwd],'replace'));
                        
                        orient landscape
                        print(hprint(II),'-djpeg',save_str);
                        
                        %Plot zooms if desired
                        figure(hprint(II));
                        Tabs_limm=get(gca,'xlim');
                        Tabs_frame=unique([Tabs_limm(1):plot_interval:Tabs_limm(2) Tabs_limm(2)]);
                        
                        
                        for JJ=1:(length(Tabs_frame)-1)
                            xlim([Tabs_frame(JJ) Tabs_frame(JJ+1)]);
                            datetick('x',date_tick_chc,'keeplimits');
                            save_str=sprintf('PSDzoom_%s_%s', datestr(Tabs_frame(JJ),30), datestr(Tabs_frame(JJ+1),30));
                            titlestr=sprintf('Start time: %s, End Time: %s, seconds averaged: %6.2f', ...
                                datestr(Tabs_frame(JJ),30),datestr(Tabs_frame(JJ+1),30),sec_avg);
                            title(titlestr);
                            print(hprint(II),'-djpeg',save_str);
                            
                        end
                        
                    end %II
                    
                end %%yes
                
        end  %switch
    otherwise
        error('Batch mode not recognized');
end %switch Batch Mode

    
    function [hprint,save_tag]=image_PSD(twin1,twin2)
        
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
        titlestr=(sprintf('Start time: %s, End Time: %s, seconds averaged: %6.2f',datestr(twin1,30),datestr(twin2,30),sec_avg));
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

if exist('def_old','var')
    def=def_old;
    def{3}=num2str(1000*fmin);
    def{4}=num2str(1000*fmax);
    
else
    def = {'Date/Time','2*3600',num2str(1000*fmin),num2str(1000*fmax),'[70 120]','4*3600','[0.01 0.1 .25 .5 .75 0.9 .99]'};
    
end
yes=1;
while yes
    prompt = {'Time unit','Time increment (sec):','Min Frequencies (Hz) [Examples: ''[100 200 300]" or "100:10:300"]:', ...
        'Max Frequencies(Hz) [Examples: ''[100 200 300]" or "100:10:300"]','dB limits [min max]','tick increment (sec)','percentile'};
    dlg_title = 'Input for PSD statistics';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        pms=[];
        return
    end
    
    pms.x_label=answer{1};
    pms.x_inc=datenum(0,0,0,0,0,str2num(answer{2}));
    pms.fmin=str2num(answer{3});
    pms.fmax=str2num(answer{4});
    pms.y_limits=str2num(answer{5});
    pms.xlabel_inc=datenum(0,0,0,0,0,str2num(answer{6}));
    pms.percentiles=eval(answer{7});
    
    if rem(length(pms.percentiles),2)==0
        uiwait(msgbox('Must be an odd number of percentiles'));
    else
        yes=0;
    end
    
    pms.def=answer;
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
            print(hprint,'-djpeg',save_str);
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
        
        write_Java_script('PSD',param);
        ! ./masterPSD.scr > outt.txt &
        return
    case 'Load Bulk Processing'
        
        
        Batch_vars_bulkload.start_time	=	'folder start';
        Batch_desc{1}	=	'Time to begin loading (e.g. "here" to use visible start time, "datenum(2011,1,2,0,0,0)", "file start" for current file , "select" to choose file from list, "folder start" for first file in folder)';
        Batch_vars_bulkload.end_time	=	'folder end';
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
            print(hprint,'-djpeg',save_str);
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

dialog_title	=	'Select an example Bulk file to load: Note that I am looking in analysis directory, not data directory';
if isfield(handles,'myfile')
    [~,token,extt] = fileparts(handles.myfile);
else
    token=[];
end


[FileName,DirectoryName] = uigetfile([token '*.psd'],dialog_title);
if isnumeric(FileName)
    disp('No file selected');
    return;
end

fprintf('Changing working directory to %s\n',DirectoryName);
cd(DirectoryName);
Iscore=min(strfind(FileName,'_chan'));

token2=FileName((Iscore+1):end);  %PSD files having the same input parameters..

Other_FileNames=dir(['*_' token2]);
if isempty(Other_FileNames)
    uiwait(errordlg('load_PSD_Bulk_Run cannot find matches to file you selected:','Files not found!'));
    return;
end


for II=1:length(Other_FileNames)
    [~,~,~,~,params]=read_Java_PSD(Other_FileNames(II).name,0,Inf,1); %Read header only
    tabs_file_start(II)=params.tstart_file;
    tabs_file_end(II)=params.tend_file;

end

tabs_folder_start=tabs_file_start(1);  %datenumber of start of data in folder (assuming all filenames have chronological order)
tabs_folder_end=tabs_file_end(end);  %datenumber of end of data in folder (assuming all filenames have chronological order)


[~,FF,~,~,~]=read_Java_PSD(Other_FileNames(1).name,0,10);  %Get FF vector...

% [~,FF,~,~,params]=read_Java_PSD(Other_FileNames(Icurrent_file).name,0,Inf,1);
% tabs_end=params.tend_file;  %datenumber of end of data in focal file (assuming all filenames have chronological order)
% tabs_start=params.tstart_file;  %datenumber of start of data in focal file (assuming all filenames have chronological order)



%%Translate start time.  Need to define tabs_start and Icurrent_file
switch lower(Batch_vars_bulkload.start_time)
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

% Asks user for specified batch processing parameters
function	Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, Batch_type)

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
dlgTitle	=	['Please enter the required parameters for processing ' Batch_type];


%	Create input dialogbox
answer		=	inputdlg(Batch_desc, dlgTitle, num_lines, defaults,'on');


%	Put variables into output structure
if	isempty(answer)
    Batch_vars	=	[];
    return;
end

for	ii	=	1:N_fields
    Batch_vars.(names{ii})	=	answer{ii};
end
end


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
    set(hObject,'String', datestr(new_date, 0));
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
set(handles.edit_datestr,'String',datestr(datenumm,0));

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
Ichan='all';
%[x,t,Fs,tstart,junk,hdr]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,Ichan,handles);
[x,~,Fs,~,~,hdr]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,Ichan,handles);

if ~isempty(strfind(lower(computer),'mac'))
    Islash=strfind(handles.mydir,'/');
else
    Islash=strfind(handles.mydir,'\');
end
chan=get(handles.edit_chan,'String');

save_name=sprintf('soundsamp%s_%s',handles.mydir((Islash(end)+1):end),datestr(tdate_start,30));
disp(['Saving ...' save_name]);

save_path	=	fullfile(pwd, save_name);  %AARON: save to local directory, not server

try
    if size(x,2)>size(x,1)
        x=x';
    end
    for Ichan=1:size(x,2)
        xfilt(:,Ichan)=filter(handles.b,1,x(:,Ichan));
    end
    wavwrite(xfilt/(1.1*max(max(abs(xfilt)))),Fs,save_path);
    
catch
    disp('No filtering desired... exists; saving raw acoustic data to WAV');
    xfilt=[];
    wavwrite(x/(1.1*max(max(abs(x)))),Fs,save_path);
    
end
save(save_path,'x','xfilt','Fs','tdate_start','hdr');

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
        %[x,t,Fs,tstart,junk1,head]=load_data(handles.filetype,handles.tdate_min , tdate_start,tlen,'all',handles);
        [x,~,Fs,~,~,head]=load_data(handles.filetype,handles.tdate_min , tdate_start,tlen,'all',handles);
        
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
            set(hh,'fontweight','bold','fontsize',14,'ycolor','y')
            
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
                print(I,'-djpeg',[save_path '.jpg']);
                save([save_path '.mat'],'x','Fs');
                close(I);
            end
            
        end
        
        return
    end
end  %MDAT

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if ~strcmp(handles.display_view,'New Fig')
    figchc=[];
    axes(handles.axes1);
else
    figure(gcf);
    chcc=get(0,'child');
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
save_path	=	fullfile(pwd, save_name);
disp(['Printing %s ...' save_name]);

orient landscape
print(figchc,'-djpeg',save_path);
%print(figchc,'-dtiff',save_path);

%print('-depsc',save_path);

end

% --- Executes on button press in pushbutton_selectpoints.
function pushbutton_selectpoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt={'Number of selections'};
def={'3'};
start_time=handles.tdate_start;
lineNo=1;
answer=inputdlg(prompt,'',lineNo,def);
tmp=ginput(str2num(answer{1}));
disp('Absolute times:')
disp(datestr(start_time+datenum(0,0,0,0,0,tmp(:,1)),0));
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
[x,~,Fs,tstart]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,Ichan,handles);

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

set(hObject,'SelectionChangeFcn',@uipanel_type_SelectionChangeFcn);


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

prompt1={'File of locations','Site','UTM location to compute range'};
dlgTitle1='Parameters for GSI localization...';
%def1={'/Volumes/ThodePortable2/2010_Beaufort_Shell_DASAR/DASAR_locations_2010.mat', '5','[4.174860699660919e+05 7.817274204098196e+06]'};
%'/Volumes/Data/Shell2010_GSI_Data/DASARlocations/DASAR_locations_2010.mat',
def1={[filesep fullfile('Volumes','Data','Shell2010_GSI_Data','DASARlocations','DASAR_locations_2010.mat')], '5','[4.174860699660919e+05 7.817274204098196e+06]'};
answer=inputdlg(prompt1,dlgTitle1,1,def1);
locs=load(answer{1});
Isite=str2double(answer{2});
Ikeep=find(~isnan(theta)&theta>0);  %Remove absent or unselected DASARs
Igood=find(~isnan(theta)&theta>-2);  %Remove absent DASARS (used to compute center of array)
%theta=theta(Ikeep);
%kappa=kappa(Ikeep);
%Note that the Site contains all DASAR locations,
%  whether they collected data or not
DASAR_coords=[locs.Site{Isite}.easting locs.Site{Isite}.northing];
[VM,Qhat,~,outcome] = vmmle_r(theta(Ikeep)',DASAR_coords(Ikeep,:),'h',kappa(Ikeep)');
mean_coords=mean(DASAR_coords(Igood,:));
CRITVAL=4.60517; %chi2inv(0.90,2);

[~,A,B,ANG,Baxis] = ellipsparms(Qhat,CRITVAL,mean_coords,VM);

figure
Ikeep=find(~isnan(theta(Igood))&theta(Igood)>0);
[~,xg,yg]=plot_location(DASAR_coords(Igood,:),theta(Igood),Ikeep,VM,A,B,ANG);
hold on
VA_cords=str2double(answer{3})/1000;
tmp=(VA_cords-[xg yg]);
plot(tmp(1),tmp(2),'o');
range=sqrt(sum((VA_cords-VM/1000).^2));
title(sprintf('Range of source from chosen location: %6.2f +/- %3.2f km, minor axis %3.2f km',range,A/1000,B/1000));
yes=menu('Print and Save?','Yes','No');
if yes==1
    orient landscape
    tstart=datestr(datenum(0,0,0,0,0,tsec)+handles.tdate_start,30);
    figure(1);
    print(1,'-djpeg',sprintf('Localization_G_%s.jpg',tstart));
    print('-djpeg',sprintf('Spectrogram_G_%s.jpg',tstart));
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
[x,~,Fs,~,~,head]=load_data(handles.filetype,handles.tdate_min , tdate_start,tlen,Ichan,handles);
x=x.';

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
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
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
            B=conventional_beamforming(Ksout.Kstot,angles,Ksout.freq,head.geom.rd,1495);
        case 3
            B=MV_beamforming(Ksout.Kstot,angles,Ksout.freq,head.geom.rd,1495);
        case 4
            B=conventional_beamforming(Ksout.Kstot,angles,Ksout.freq,head.geom.rd,1495);
            
            B2=MV_beamforming(Ksout.Kstot,angles,Ksout.freq,head.geom.rd,1495);
        case 5
            R=derive_reflection_coefficient2(Ksout.Kstot,angles,Ksout.freq,head.geom.rd,1495);
        case 6
            %Simple Pekeris waveguide determined from DASAR cutoff frequencies
            
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
    title(sprintf('%s, %i FFT, %i elements',datestr(tdate_start),Nfft,length(head.geom.rd)));
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
    write_covmat(Ksout.freq,Ksout.Kstot,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    
    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
    write_covmat(Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    
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
[x,~,Fs,tstart,junk,head]=load_data(handles.filetype,handles.tdate_min , tdate_start,tlen,'all',handles); %#ok<*ASGLU>
if size(x,2)>1
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
[x,t,Fs,tstart,junk,head]=load_data(handles.filetype,handles.tdate_min , tdate_start,tlen,'all',handles);

%x(1) is shallowest element...
if size(x,2)>1
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
[data.x,t,Fs,tstart,junk,head]=load_data(handles.filetype,handles.tdate_min , tdate_start,tlen,'all',handles);

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
            print(gcf,'-djpeg',sprintf('DepthEstimate_%s.jpg',savename_base));
            
            
            
            
            
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
set(handles.edit_datestr, 'String', datestr(new_date,0));
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
set(handles.edit_datestr, 'String', datestr(new_date,0));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

pushbutton_update_Callback(handles.pushbutton_update,eventdata,handles);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Beginning of annoation subroutines %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	new annotation stuff

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
function pushbutton_notes_new_Callback(hObject, eventdata, handles)
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

%	Initial signal type selection
sig_types	=	{'Pulsive','FM'};
choice		=	menu('Signal type?',sig_types);
if	choice == 0
    disp('Input cancelled');
    delete(hrec);
    return;
end
sig_type	=	sig_types{choice};

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
else
    Event	=	Data.Template;
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
    handles.notes.Data.Events(i_sel)	=	[];
    N	=	length(handles.notes.Data.Events);
    if	N == 0
        i_sel	=	[];
    elseif	i_sel > N
        i_sel	=	1;
    end
    handles.notes.i_sel	=	i_sel;
    
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
function pushbutton_notes_screen_Callback(hObject, eventdata, handles)
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


%%	Supporting functions, i.e. not auto-generated callbacks
function	status	=	dependency_check()

status	=	true;

end

function	handles	=	load_and_display_spectrogram(handles)

cla;

%tdate_start		=	handles.tdate_start;
tlen	=	handles.tlen;
if tlen<1e-2 %somehow window length is way too short
    tlen=1;
    set(handles.edit_winlen,'String',num2str(tlen));
end
if isempty(tlen)
    tlen    =   str2double(get(handles.edit_winlen,'String'));
    
    handles.tlen=tlen;
end
mydir	=	pwd;
Ichan	=	str2double(get(handles.edit_chan,'String'));  %Hardwire first channel

try
    [x,t,Fs,tstart,junk,hdr]=load_data(handles.filetype,handles.tdate_min,...
        handles.tdate_start,tlen,Ichan,handles);
    
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

handles.display_view=get(get(handles.uipanel_display,'SelectedObject'),'String');

if strcmp(handles.display_view,'Spectrogram')||strcmp(handles.display_view,'New Fig')
    
    if strcmp(handles.display_view,'Spectrogram')
        axes(handles.axes1);
    else
        figure;
    end
    
    if length(x(:,1))<Nfft
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
    
    %%Check that we are looking at acoustic data
    if isempty(strfind(handles.myfile,'Press'))
        disp('Enter min and max frequency, or return to see raw time series:');
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
        y=x(:,1)/1000; %Pressure conversion
    end
    
    t=(1:length(x(:,1)))/Fs;
    plot(t,y);grid on;
    xlabel('Time (sec)');ylabel('Amplitude');
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
end

handles.x	=	x;
handles.Fs	=	Fs;
%tmp=ginput(2)

if strcmpi(handles.filetype,'gsi')
    set(handles.pushbutton_GSIbearing,'vis','on');
    set(handles.pushbutton_GSI_localization,'vis','on');
else
    set(handles.pushbutton_GSIbearing,'vis','off');
    set(handles.pushbutton_GSI_localization,'vis','off');
end

if strcmpi(handles.filetype,'mdat')||strcmpi(handles.filetype,'wav')||strcmpi(handles.filetype,'mat')
    set(handles.pushbutton_CSDM,'vis','on');
    set(handles.pushbutton_Mode,'vis','on');
    set(handles.pushbutton_tilt,'vis','on');
    set(handles.pushbutton_modalfiltering,'vis','on');
else
    set(handles.pushbutton_CSDM,'vis','off');
    set(handles.pushbutton_Mode,'vis','off');
    set(handles.pushbutton_modalfiltering,'vis','off');
    set(handles.pushbutton_tilt,'vis','off');
    
end

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
end

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
    [x,t,Fs,tmin,tmax]=load_data(filetype,-1,-1,1,1,handles);
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
small_step	=	5/T_len;				% 5s
big_step	=	0.1;					% 10%
set(handles.slider_datestr,'sliderstep',[small_step big_step]);

set(handles.slider_datestr,'Value',0.5);
handles.tdate_start		=	0.5*(tmin+tmax);
set(handles.edit_datestr,'String',datestr(handles.tdate_start,0));

set(handles.text_mintime,'String',datestr(tmin,0));
set(handles.text_maxtime,'String',datestr(tmax,0));
handles.tdate_min	=	tmin;
handles.tdate_max	=	tmax;

set(handles.edit_fmax,'String',Fs/2000);
set(handles.edit_fmin,'String',0);
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


%%% load_data.m%%%%
function [x,t,Fs,tmin,tmax,head]	=	...
    load_data(filetype,tstart_min,tdate_start,tlen,Ichan,handles)
%%% tmin,tmax, t are datenumbers
%   x rows are samples, columns are channels

persistent  keyword
mydir=handles.mydir;
myfile=handles.myfile;
teager=get(handles.checkbox_teager,'Value');
set(handles.edit_normal_rotation,'Vis','Off');  %Set off for the moment
set(handles.text_normal_rotation,'Vis','Off');  %Set off for the moment
            

x=[];
t=[];
tmin=[];
tmax=[];
head=[];
filetype	=	upper(filetype);

switch filetype
    case 'PSD'
        [x,F,t,Tabs,params]=read_Java_PSD(fullfile(mydir,myfile),tdate_start,tlen);
        if ~isempty(t)
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
    case 'MT'
        %[x,t,Fs]=load_mt_mult(handles.mydir,tdate_start,tlen);
        head=read_mt_header([mydir filesep myfile]);
        
        tmin=head.tstart;
        tmax=head.tend;
        Fs=head.Fs;
        
        if tdate_start>0
            tdate_vec=datevec(tdate_start-tmin);
            nsec=tdate_vec(6)+60*tdate_vec(5)+3600*tdate_vec(4);
        else
            %uiwait(msgbox('No start time requested'));
            return
        end
        %Load accelerometer data as desired
        %% Template: CB_Accel_X_2014_02100958.mt
        if strfind(head.abbrev,'Accel')
            mystr='XYZ';
            Ispace=strfind(myfile,'_');
            Ispace=Ispace(2)+1;
            for Iacc=1:3
               myfile_temp=myfile;
               myfile_temp(Ispace)=mystr(Iacc);
               
               [x0,t]=load_mt([mydir filesep myfile_temp],nsec,tlen);
               if Iacc==1
                   x=zeros(length(x0),3);
               end
               x(:,Iacc)=x0;
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
             myfile(Ispace)=mystr(Ichan);
             
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
                [x,t]=load_mt([mydir filesep myfile],nsec,tlen);
            end
            head.Nchan=1;
        end
    case 'SIO'
        %x,t,Fs,tmin,tmax,head
        Fs=input('Enter a sampling frequency in Hz:','s');
        [~, head] = sioread(fullfile(mydir,myfile),1,[],[]);
        
        %%Assume start time is of form 14047225500
        
        %%	Data parameters needed
        Np		=	Header.PperChan;
        Nchan	=	Header.N_Chan;
        
        
        
        %-----------------------------------------------------------------------
        % sioread.m
        %
        % This program runs under windows, unix, and macs.
        %
        % function x=sioread(filename,p1,npi,channels);
        %
        % Inputs:
        % 	filename: Name of sio file to read
        % 	p1:	Point to start reading ( 0 < p1 < np)
        % 	npi: 	Number of points to read in
        % 	channels: Single number or vector containing the channels to read
        % 		(example-to read channels 4 thru 10 enter 4:10)
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
        t=(1:length(x))/Fs;
        
        tmin=head.tfs;
        tmax=head.tfe;
        head.Nchan=size(x,2);
        if strcmp(Ichan,'all')
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
        Nsamples	=	wavread([mydir '/' myfile],'size');
        [~,Fs]		=	wavread([mydir '/' myfile],1,'native');
        Nsamples	=	Nsamples(1);
        handles.Fs	=	Fs;
        
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
                tmin	=	convert_date(myfile,'_');
                tmax	=	convert_date(myfile,'_') + datenum(0,0,0,0,0,Nsamples/Fs);
            catch
                disp([myfile ': convert_date failure']);
                minn	=	input('Enter start date in format [yr mo day hr min sec]: ');
                if isempty(minn)
                    minn=zeros(1,6);
                end
                tmin	=	datenum(minn);
                tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
            end
        else
            tmin		=	tstart_min;
            tmax		=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
            tdate_vec	=	datevec(tdate_start - tmin);
            nsec		=	tdate_vec(6) + 60*tdate_vec(5) + 3600*tdate_vec(4);
            N1			=	1 + round(nsec*handles.Fs);
            N2			=	N1 + round(tlen*handles.Fs);
            [x,Fs]		=	wavread([mydir '/' myfile],[N1 N2],'native');
            
            if ~strcmp(Ichan,'all')
                x		=	x(:,Ichan);
            end
            
            t	=	(1:length(x))/Fs;
        end
        x			=	double(x)*sens;
        head.Nchan	=	size(x,2);
        
end


%%%Optional Teager-Kaiser filtering...
if teager
    %%Assume that x is in form [ channel time]
    x=x(:,2:end-1).^2-x(:,1:end-2).*x(:,3:end);
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

end

function [thet,kappa,tsec]=get_GSI_bearing(hObject,eventdata,handles)
%tsec: seconds into time series that data are selected...
thet=-1;  %Start with failed result
kappa=-1;
tsec=-1;
tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
%Ichan=str2double(get(handles.edit_chan,'String'));
[x,t,Fs,tstart,tend,head]=load_data(handles.filetype,handles.tdate_min,tdate_start,tlen,'all',handles);

disp('Click on two extreme corners, click below axis twice to reject:');
tmp=ginput(2);
if (isempty(tmp)||any(tmp(:,2)<0))
    return
end
tsec=min(tmp(:,1));
n=round(Fs*sort(tmp(:,1)));

freq=1000*sort(tmp(:,2));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

[thet0,kappa,sd]=extract_bearings(x(n(1):n(2),:),0.25,Nfft,Fs,freq(1),freq(2),50);

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
end

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

function [mux,Rbar,kappa] = vm_ests_uneq(x,options,flag) %#ok<STOUT>
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
idx = find(lx>0);                        % Get rid of 0-length vectors
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

function [Kstot,f,VV,EE_sort]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,Isnap,M,nowin,threshold,tiltdata)

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
if	~exist('Isnap', 'var')|| Isnap<0
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
        if exist('yesnorm', 'var')
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
        
        print(gcf,'-djpeg',[savename '.jpg']);
    end
end %tilt
end

function [MFP_replica_file_out,range_str_out,default_tilt_out,default_SNR_out]=load_MFP_scenario

%%MFP replicas can be divided into Pekeris models, and more complex models...

if strcmp(getenv('USER'),'thode')
    model_dir='/Users/thode/Projects/Arctic_2010/PropagationModeling.dir/';
elseif strcmp(getenv('USER'),'shabadi');
    model_dir='/Users/shabadi/Documents/Shima-SIO/PropagationModeling.dir/';
elseif strcmp(getenv('USERNAME'),'bonnelju');
    model_dir=fullfile('C:','Users','bonnelju','Desktop','SanDiego','PropagationModeling.dir');
    %model_dir='C:\Users\bonnelju\Desktop\SanDiego\PropagationModeling.dir\';
end

%%%%%%%%Short-Range Estimate%%%%%%%%%%%%%%%%
MFP_replica_file{2}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_10to500Hz.mat');
%MFP_replica_file{2}{1}=[model_dir 'Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_10to500Hz.mat'];
default_tilt{2}{1}='1.94';
range_str{2}{1}='100:5:5000';
default_SNR{2}{1}=[45.78 76.29 88.5 170.9 210.6 253.3];

MFP_replica_file{2}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_ShortRange_arctic_pekeris1707_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{2}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_ShortRange_arctic_pekeris1707_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{2}{2}='2.41';
range_str{2}{2}=range_str{2}{1};
default_SNR{2}{2}=default_SNR{2}{1};

%%%%%%%%%%%%%Medium range 7 km est%%%%%%%%%%%%%%%%
MFP_replica_file{3}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_10to500Hz.mat');
%MFP_replica_file{3}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_10to500Hz.mat'];
%/Users/thode/Projects/Arctic_2010/MidRangeMFPInversion.dir/Thermocline_11dB_eigenvector_inversion.dir/RunAA
%SSP EOF included, 40 pop, thermocline depth fixed, depth restricted to 55 m, gradient greater than -50 m/s
%receiver depth limited to +/1 m/s.
default_tilt{3}{1}='1.35';
range_str{3}{1}='5000:100:10000';
default_SNR{3}{1}=[112.9 116 119 122.1 164.8 167.8 170.9 174 177 180.1 235 238];  %Original 12 frequencies
%used for inversion
%default_SNR{3}{1}=[116 119 122.1 164.8 167.8 170.9 174 177 180.1 ]; %Plotting results

MFP_replica_file{3}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{3}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{3}{2}='1.585';
range_str{3}{2}=range_str{3}{1};
default_SNR{3}{2}=default_SNR{3}{1};

MFP_replica_file{3}{3}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{3}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{3}{3}='1.585';
range_str{3}{3}=range_str{3}{1};
default_SNR{3}{3}=default_SNR{3}{1};

%%%%%%%%%%%%%%%%%Long range estimate (17 km)%%%%%%%%%%%%%%%%%%%%%%%%
MFP_replica_file{4}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_10to500Hz.mat');
%MFP_replica_file{4}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_10to500Hz.mat'];
default_tilt{4}{1}='1.35';
range_str{4}{1}='5000:100:20000';
%default_SNR=12;
default_SNR{4}{1}=[74.01 78.58 85.45 88.5 90.03 92.32 93.08 95.37 99.95 103 107.6 109.9];

MFP_replica_file{4}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_FarRange_arctic_pekeris1692_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{4}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_FarRange_arctic_pekeris1692_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{4}{2}='1.585';
range_str{4}{2}=range_str{4}{1};
default_SNR{4}{2}=default_SNR{4}{1};

%%%%%%%%%%%%%%%%%35 km range estimate%%%%%%%%%%%%%%%%
MFP_replica_file{5}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105748_10to500Hz.mat');
%MFP_replica_file{5}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105748_10to500Hz.mat'];

default_tilt{5}{1}='3.93';
range_str{5}{1}='30000:100:40000';
range_str{5}{1}='100:100:40000';
default_SNR{5}{1}=[99.1800  102.2000  105.3000  108.3000  111.4000  114.4000  117.5000  120.5000 ...
    123.6000  126.6000  129.7000  132.8000  135.8000  138.9000  141.9000  145.0000 ...
    148.0000  151.1000  154.1000  157.2000  160.2000 ];
%default_SNR=default_SNR(1:(end-1));
% default_SNR=default_SNR((end-3):end)

MFP_replica_file{5}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_UltraRange_arctic_pekeris1553_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{5}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_UltraRange_arctic_pekeris1553_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{5}{2}='3.34';
range_str{5}{2}=range_str{5}{1};
default_SNR{5}{2}=default_SNR{5}{1};


%Averaged Pekeris model...
for I=2:5
    %%Averaged Pekeris model
    MFP_replica_file{I}{3}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_arctic_pekeris1673_KRAKEN_flat55m_10to500Hz.mat');
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

if ~exist('linel', 'var')
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
    kest = true;
else
    kest = false;
end
r = lower(r);
robust = strfind('mah',r);
if isempty(robust),
    error('Input r must be either "m", "a", or "h"');
end
failed = {'less than 2 bearings','negative variance estimates',...
    'solution behind DASARs','failed to converge'};
n = length(angle);

if n<=1  %If only one set of bearings present...
    outcome = failed{1};
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
        while (~converge)&&(iter<maxiter),
            iter = iter+1;
            xyold = xyhat;
            % d = dist(xyhat', [x y]');
            
            for JJ=1:length(x)
                d(JJ)=sqrt((xyhat(1)-x(JJ)).^2+(xyhat(2)-y(JJ)).^2);
            end
            if (robust>1) && (n>2), % Need 3 or more bearings to calculate weights for
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
            if ((n-sum(~w))>1) && (cond2(M1)<cond_num),
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
        if kest && (n>2),
            Cbar = (w*Cd'/sum(w))^(n/(n-2));   % Exponent is small sample size correction
            k = inv(2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar);
        end
        if Cd*w'>0,             % Weighted average of cosine differences (check
            VM = xyhat';          %   on bad solutions behind DASARs)
            if kest && (n==2),     % Cannot estimate Qhat with only 2 bearings
                Qhat = nan*ones(2);
                outcome = 'successful; 2 bearings; no Kappa';
            else
                k = k';
                cv = -(k.*sstar*c + k.*cstar*s)/2;
                M3 = [k.*sstar*s cv; cv k.*cstar*c];
                Qhat = inv(M3);
            end
            if ~kest || (n>2),
                if all(diag(Qhat)>0),
                    outcome = 'successful';       % Successful solution
                else
                    outcome = failed{2};        % Implausible variance estimate(s)
                end
            end
        else
            outcome = failed{3};          % Bad solution behind DASARs
        end % if all(Cd>0)
    else
        outcome = failed{4};            % No convergence
    end   % if converge
end     % if n<=1


if isempty(strfind(outcome,'successful'));
    %if ~isempty(strfind(failed,outcome))
    VM = [nan nan];
    Qhat = nan*ones(2);
    w = zeros(1,n);
end
end

function [brefa_table,Icol]=calibrate_bearing_Shell2007(cal_dir,fname,no_load_table)
%function [brefa_table,Icol]=calibrate_bearing_Shell2007(cal_dir,fname,no_load_table)
% Return nonlinear calibration data for 2007 Dasars
% Input:
%          cal_dir: base directory of calibration information
%          fname: string similar in structure to 'S108H0T20080921T000000'
%          no_load_table:  if exists, don't load table.
% Note:  the calibration is assigned based on name of file, not name of enclosing directory..
brefa_table=[];
Islash=1+max(strfind(fname,'/'));
if ~isempty(Islash)
    fname=fname(Islash:end);
end
site=str2double(fname(2));
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
    data=load(cal_dir);
    brefa_table=data.bgridt;
    
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

function audio_stop(player, ~, handles)
%	Set controls
ispaused	=	(player.CurrentSample > 1) &&...
    (player.CurrentSample < player.TotalSamples);
if	~ispaused
    set(handles.pushbutton_playsound, 'String', 'Play');
    set(handles.pushbutton_pausesound, 'String', 'Pause');
    set(handles.pushbutton_pausesound, 'Enable', 'off');
end

end

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
end

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

end  %disable

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


end

%	Checks Notes folder for existing files and loads them
%%%%%%Thode: load_notes_file
function	handles		=	load_notes_file(handles, new_folder)


%	First check that current notes are saved before proceeding
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


%	New notes file disables all notes controls, until otherwise activated
opt		=	'off';
disable_notes_nav(handles, opt);
set(handles.pushbutton_notes_new, 'Enable', opt);
set(handles.pushbutton_notes_save, 'Enable', opt);
set(handles.pushbutton_notes_edit, 'Enable', opt);
set(handles.checkbox_notes_show, 'Enable', opt);
set(handles.checkbox_notes_delete, 'Enable', opt);


%	Switch to new folder and check for existing files
if	exist('new_folder','var') && ~isempty(new_folder)
    folder_name	=	new_folder;
else
    folder_name	=	handles.notes.folder_name;
end
[~,fname,~]		=	fileparts(handles.myfile);

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
%%AARON:  Give user option of loading automated detector
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
    
end


%	New file with default data template
[Defaults.Description, Defaults.Template, edit_fields]	=	load_default_template();
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
        
    end  %ii loop through all names.
    
    %	Always merge with current template data, to update old notes with
    %	new fields
    Data	=	merge_data(Defaults, Data);
    Data.Events(1)	=	[];
    
    
    if	length(Sel) > 1 & manual_flag
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



%	enable relevant buttons
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
handles		=	load_gui_params(handles, GUI_params);

%	Set folder text box to selected directory + file_name
set(handles.edit_folder, 'String', handles.notes.file_path);

%	Set note index text accordingly
i_sel		=	handles.notes.i_sel;
N			=	length(handles.notes.Data.Events);
msg			=	{[num2str(i_sel) ' of']; num2str(N)};
set(handles.text_notenum, 'String', msg);

end

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
set(handles.edit_datestr, 'String', datestr(GUI_params.tdate_start));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

end

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

end

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
for	ii	=	1:N_fields
    value	=	Event.(names{ii});
    defaults{ii}	=	num2str(value);
end

%	Size of each prompt field
num_lines		=	ones(N_fields,1);
num_lines(end)	=	2;

%	title for window
dlgTitle	=	['Annotation for event at ' datestr(OldEvent.start_time)];

if	length(Description) ~= N_fields
    error('# of Descriptors differs from # of fields');
end

%	Create input dialogbox
options.Resize	=	'on';
answer		=	inputdlg(Description, dlgTitle, num_lines, defaults, options);


if	isempty(answer)
    NewEvent	=	[];
    return;
end

NewEvent	=	OldEvent;
for	ii	=	1:N_fields
    NewEvent.(names{ii})	=	answer{ii};
end


end

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
i_time_noise2	=	(Tmax < T) & (T <= Tmax+duration_noise);
%Shift one column earlier to avoid picking up first signal FFT bin
if (sum(i_time_noise1)*dT < 0.9*duration_noise) ||...
        (sum(i_time_noise2)*dT < 0.9*duration_noise)
    disp('Selected signal too close to start of window; SNR cannot be computed');
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


end

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

end

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
        Data.Events(jj).(names{ii})		=	num2str(Data.Events(jj).(names{ii}));
    end
end
%close(hh)


end

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


end

%	Remove obsolute fields from old events, but save data into comments
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


end

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
h_show	=	[];
sel_vis	=	false;
axes(h_axes);
for	ii	=	1:length(i_show)
    ie		=	i_show(ii);
    event	=	Events(ie);
    x		=	event.start_time - Times(1);	x	=	x*24*60*60;
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
guidata(h_axes, handles);

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

end

%	Pops up new window with event details
function	hdlg	=	show_event_info(Event, Description, hdlg)

if nargin < 3 || isempty(hdlg) || ~ishandle(hdlg)
    hdlg	=	dialog('Windowstyle', 'normal');
    pos		=	get(hdlg, 'Position');
    pos(3)	=	pos(3)/2;
    set(hdlg,'Position',pos);
else
    figure(hdlg);
    clf(hdlg);
end

%	Fetch details
Title	=	'Annotation details';
names	=	fieldnames(Event);
Message	=	[Description(:).'; names(:).'];
for	ii	=	1:length(names)
    temp			=	num2str(Event.(names{ii}));
    Message{2,ii}	=	[char(9) char(9) temp(:).'];
end
%	comments don't need to be indented
Message{2,ii}	=	temp;
%	start_time should be formated
Message{2,1}	=	datestr(Event.(names{1}), 'yyyy-mm-dd HH:MM:SS.FFF');

%	Put details in window
set(hdlg, 'Name', Title);
htxt	=	uicontrol(hdlg, 'Style', 'edit',...
    'Enable', 'on',...
    'Units', 'normalized',...;
    'Position', [0 0 1 1],...
    'Max', 2, 'Min', 0,...;
    'HorizontalAlignment','left',...
    'String', Message(:));


end

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
end

%	updates whole figure window according to selected event
function	handles		=	update_events(handles)

i_sel		=	handles.notes.i_sel;
start_time	=	handles.notes.Data.Events(i_sel).start_time;
tlen		=	datenum(0,0,0,0,0,handles.tlen);

%	New start date for window
new_date	=	start_time - tlen/2;
set(handles.edit_datestr, 'String', datestr(new_date));
edit_datestr_Callback(handles.edit_datestr, [], handles)
handles		=	guidata(handles.edit_datestr);

%	update figure window
pushbutton_update_Callback(handles.pushbutton_update, [], handles);
handles		=	guidata(handles.pushbutton_update);

end

%	Helper function to get current users id/name
function	user_name	=	getusername()
if	isunix || ismac
    user_name = getenv('USER');
elseif ispc
    user_name = getenv('USERNAME');
else
    error('Unregognized system');
end
end

%	helper function for default template
function	[Description, Template, edit_fields]	=	load_default_template()
%edit fields are fields you can edit...
Description	=	{'Start Time',...
    'Author',...
    'Pulse or FM?',...
    'Call Type',...
    'Min Freq (Hz)',...
    'Max Freq (Hz)',...
    'Duration (s)',...
    'Noise SEL (dB re 1 uPa^2-s)',...
    'Noise rms (dB re 1 uPa)',...
    'Noise peak PSD (dB re 1uPa^2/Hz)',...
    'Signal SEL (dB re 1 uPa^2-s)',...
    'Signal rms (dB re 1 uPa)',...
    'Signal peak PSD (dB re 1uPa^2/Hz)',...
    'Signal SNR (rms)', ...
    'Signal SNR (dB rms)', ...
    '# of pulses',...
    '# of harmonics',...
    'Modulation (Hz)',...
    'Confidence (1-5)',...
    'Comments'};

%	Default values
Template.start_time			=	0;
Template.author				=	'Your name';
Template.sig_type			=	'NA';
Template.call_type			=	'S1';
Template.min_freq			=	0;
Template.max_freq			=	5000;
Template.duration			=	10;

Template.noise_se_dB		=	0;
Template.noise_rms_dB		=	0;
Template.noise_peakpsd_dB	=	0;

Template.signal_se_dB		=	0;
Template.signal_rms_dB		=	0;
Template.signal_peakpsd_dB	=	0;

Template.SNR_rms			=	0;
Template.SNR_rms_dB			=	0;

Template.num_pulses			=	2;
Template.num_harmonics		=	-1;
Template.modulation			=	0;
Template.confidence			=	3;
Template.comments			=	'';

edit_fields			=	fieldnames(Template);
edit_fields(5:14)	=	[];

% params_extract.noise_se_dB		=	noise_se_dB;
% params_extract.noise_rms_dB		=	noise_rms_dB; %not 20 log because actually mean sqaure
% params_extract.noise_peakpsd_dB	=	noise_peakpsd_dB;
%
% params_extract.signal_se_dB		=	signal_se_dB;
% params_extract.signal_rms_dB	=	signal_rms_dB;  %not 20 log because actually mean sqaure
% params_extract.signal_peakpsd_dB=	signal_peakpsd_dB;
%
% %params_extract.SNR_rms=SNR_rms;
% params_extract.SNR_rms_dB		=	SNR_rms_dB;


end


%%	Functions that should be removed and referenced externally instead

function x=calibrate_GSI_signal(xin, keyword,RawFileName)

%calibrate_GSI_signal.m
% Convert raw A/D value of presure sensor in DASAR to uPa
% using sensitivity of 150 dB re 1uPa/V
% and peak A/D voltage of 2.5 V
%keyboard

if isempty(xin)
    x=[];
    return
end
if strcmp(keyword,'short')||~isempty(strfind(keyword,'DASAR2007'))
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
elseif strcmp(keyword,'short')||~isempty(strfind(keyword,'DASARC'))
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
elseif ~isempty(strfind(keyword,'NorthStar08')),
    
    %%This has a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    if ~isempty(strfind(RawFileName,'NS08A0'))||~isempty(strfind(RawFileName,'NA08Cx'))
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
elseif ~isempty(strfind(keyword,'Liberty08')>0)
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
hprint=[];
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
        print(hprint(Ix),'-djpeg',save_str);
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
end

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



end


function edit_normal_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_normal_rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_normal_rotation as text
%        str2double(get(hObject,'String')) returns contents of edit_normal_rotation as a double
end

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%TIMEWARP functions below %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





