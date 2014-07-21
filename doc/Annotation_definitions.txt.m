%%Instructions and details about Ulyssess annotations and notes
%  Aaron Thode
%  July 1, 2014


A.  Basic annotation data fields

The basic annotation data fields are defined in "load_default_template",
    which are first loaded in "load_notes_file", which in turn is called by
    "pushbutton_notes_select_Callback"

B. Fields associated with notes.
    ~handles.notes.saved
    ~handles.notes.readonly
    handles.notes.folder_name;
    
    set(handles.pushbutton_notes_new, 'Enable', opt);
    set(handles.pushbutton_notes_save, 'Enable', opt);
    set(handles.pushbutton_notes_edit, 'Enable', opt);
    set(handles.checkbox_notes_show, 'Enable', opt);
    set(handles.checkbox_notes_delete, 'Enable', opt);
    
    
%%%  Annotation template, Annotations template
%	 Defines default fields
%%load_default_template
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
