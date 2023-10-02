function startup_info = gui_startup_information(app)

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
