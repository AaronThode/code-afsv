function [list_names,filter_params,station_locations]=load_bowhead_detector_params(outputdir)
%function [list_names,filter_params]=load_bowhead_detector_params(outputdir)
% list_names is a cell array of full pathnames to automated detection files
% filter_params is a structure with two-element arrays for filtering
%   automated results.  'easting' and 'northing' are built in
% station_locations:  structure with 'northing' and 'easting' fields for
%       each DASAR, from 'A'(1)  to 'G'

list_names=[];
filter_params=[];
station_locations=[];

% Location folder with automated results to convert...
dirname = uigetdir(outputdir, 'Select a folder with automated bowhead detector results...');

if dirname==0
    return
end

mydir=pwd;
cd(dirname);

% Define a template for selecting files
prompt={'Enter a template string for a location output file:', ...
    'Enter geographical restriction (West,Center,East) or 2x2 UTC matrix:'};
name='File Template';
numlines=1;
options.Resize='on';
options.WindowStyle='normal';
defaultanswer={'S5*_BearingInterval_Huber_FilteredLocations.mat','Center'};

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

% Determine UTC boundaries
Site5_easting=[
    412659
    418898
    412949
    419139
    413253
    413492];

Site5_northing=[
    7795035
    7798311
    7802045
    7805333
    7809065
    7816110];

%Save data to station_locations output variable
station_locations.northing=Site5_northing;
station_locations.easting=Site5_easting;

IDASAR=logical([0 0 1 1 1 0]);  %DASARS CDE

Site5_northing=Site5_northing(IDASAR);
Site5_easting=Site5_easting(IDASAR);

xbuffer=50; %A little extra shrinking of the center site width, in meters (makes sure box of Center is inside of all DASARS
filter_params.northing=[min(Site5_northing)+xbuffer max(Site5_northing)-xbuffer];  %Boxes at same latitude
dl=abs(min(Site5_easting)-max(Site5_easting));  %How wide should a box be?  For now, make width the same as the width of the 'center' box..
%dl=10000; %10 km box width

filter_params.easting=[min(Site5_easting)+xbuffer max(Site5_easting)-xbuffer];

switch(answer{2})
    case 'Center'
        %Keep the same
    case 'West'
        filter_params.easting=filter_params.easting-(2*xbuffer+dl);
    case 'East'
        filter_params.easting=filter_params.easting+(2*xbuffer+dl);
    otherwise
        return
        
end
filter_params.keyword=answer{2};
