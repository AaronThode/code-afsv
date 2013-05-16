%%%%%%%%%gui_startup_information.m%%%%%%%%%%%%%%%%
% Aaron Thode
% November 22, 2004
% Contains directory locations and names of various files
%   used by the Schultz GUI system
function [startup_info]=gui_startup_information

%%%%Directory where GSI_specgram_viewer.m and function_handles files will be
%%%%stored....
[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);


if strcmp(local_machine,'macmussel.ucsd.edu')

    startup_info.base_directory='/Users/Shared/MATLAB/AllFile_specgram_viewer';


    startup_info.default_directory='/Volumes/Shared/Projects';
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';

    %%Default annotation file
    startup_info.annotation_file='annotated.txt';
     startup_info.function_handles_filename='mt_specgram_handles.mat';

    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';
elseif strcmpi((computer),'maci')  %KatMacGSI
   

elseif strcmpi((computer),'mac')

    startup_info.base_directory='//Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/Arctic_2007/CommonScripts.dir';
    startup_info.default_directory='/Volumes/GSI_08/';
    startup_info.default_inputfiledir='/Users/thode/Projects/Insta-array/Alaska_Sperm';
    startup_info.annotation_file='annotated.txt';

    startup_info.function_handles_filename='mt_specgram_handles.mat';

    %%Name of default multipath storage file
    startup_info.multipath_filename='gui_multipath_selections.mat';

end