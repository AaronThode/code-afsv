function [list_names,filter_params]=load_airgun_detector_params(outputdir)
%function [list_names,filter_params]=load_airgun_detector_params(outputdir)
% list_names is a cell array of full pathnames to automated detection files
% filter_params is a structure with two-element arrays for filtering
%   automated results.  'easting' and 'northing' are built in

list_names=[];
filter_params=[];

% Location folder with automated results to convert...
dirname = uigetdir(outputdir, 'Select a folder with automated airgun detector results...');

if dirname==0
    return
end

mydir=pwd;
cd(dirname);

% Define a template for selecting files
prompt={'Enter a template string for a airgun detector file:'};
name='File Template';
numlines=1;
options.Resize='on';
options.WindowStyle='normal';
defaultanswer={'S*_airgun*.mat'};

answer=inputdlg(prompt,name,numlines,defaultanswer,options);
fnames=dir(answer{1});

% Look for possible matches in folder, and select all relevent files
for I=1:length(fnames)
    list_names{I}=fnames(I).name;
end

if isempty(fnames)
    uiwait(msgbox('No files found!','modal'));
    return
end

[Sel, OK]	=	listdlg('ListString', list_names,...
    'Name', 'Available files',...
    'PromptString', 'Select file(s) to load:',...
    'OKString', 'Load',...
    'CancelString', 'New', ...
    'ListSize',[320 300]);

list_names=list_names(Sel);
for I=1:length(list_names)
    list_names{I}=fullfile(dirname,list_names{I});
end
cd(mydir);



