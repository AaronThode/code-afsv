% load_default_annotation_template.m%%%
%  Aaron Thode, July 22, 2014
%%%  Annotation template, Annotations template
%	 Defines default fields
%%load_default_template
function	[Description, Template, edit_fields]	=	load_default_annotation_template()
%edit fields are fields you can edit...
Description	=	{'Start Time',...
    'Author',...
    'HashTag', ...
    'Pulse or FM?',...
    'Call Type',...
    'Min Freq (Hz)',...
    'Max Freq (Hz)',...
    'Duration (s)',...
    'FM Slope (Hz/s)', ...
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
    'Number of multipath',...
    'Confidence (1-5)',...
    'Comments', ...
    'Link Names (if you want to link this annotation to other annotation files)', ...
    'Link HashTags (hashtags of linked annotations)'};

%	Default values
Template.start_time			=	0;
Template.author				=	'Your name';
Template.hash_tag           =   now;
Template.Istation='1';
Template.sig_type			=	'NA';
Template.call_type			=	'S1';
Template.min_freq			=	0;
Template.max_freq			=	5000;
Template.duration			=	10;
Template.slope              =   0;


Template.noise_se_dB		=	0;
Template.noise_rms_dB		=	0;
Template.noise_peakpsd_dB	=	0;

Template.signal_se_dB		=	0;
Template.signal_rms_dB		=	0;
Template.signal_peakpsd_dB	=	0;

Template.SNR_rms			=	0;
Template.SNR_rms_dB			=	0;

Template.num_pulses			=	0;
Template.num_harmonics		=	-1;
Template.multipath			=	0;
Template.confidence			=	3;
Template.comments			=	'';
Template.link_names         =   '';
Template.link_hashtags      =   '';

edit_fields			=	fieldnames(Template);
edit_fields(11:18)	=	[];

edit_fields(3:4)=[];

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