
%%%%%%%%%%load_program_location.m%%%%%%%
%  Aaron Thode
%  Provide location of JAVA program on multiple platforms...

function program_dir=load_program_location

[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);

if strcmp(local_machine,'macmussel.ucsd.edu')
    program_dir ='/Users/Shared/DataPreProcessor/LatestVersion/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'KatsMacPro'))
    %program_dir ='/Users/thode/DataPreProcessor/PreProcessV1.jar';
    program_dir ='/Users/khkim/Thode/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'thode'))
    program_dir ='/Users/thode/Documents/workspace/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'dgrebner'))
    program_dir ='/Users/dgrebner/Projects/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'melania'))
    program_dir ='/Users/melania/Projects/Matlab_IO/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'delphine'))
    %program_dir ='/Users/delphine/Delphine/DataPreProcessor/PreProcessV1.jar'; 
    program_dir='/Users/delphine/Delphine/MATLAB_SVN/DataPreProcessor/PreProcessV1.jar'; 

elseif strcmp(getenv('USER'),'dponcemorado')
    program_dir='~/Desktop/Graywhale_Part2/ThodeLabSoftware.dir/JAVA_programs/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'Janice'))|~isempty(findstr(local_machine,'alaska.edu'))
    program_dir='~/Desktop/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'uas.ad'))|~isempty(findstr(local_machine,'raptor'))
    program_dir='~/Desktop/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(findstr(local_machine,'ucsd.edu'))
    program_dir ='/Users/dponcemorado/Desktop/ThodeLabSoftware.dir/JAVA_programs/DataPreProcessor/PreProcessV1.jar';
else
	program_dir='~/Desktop/DataPreProcessor/PreProcessV1.jar';
    
    
end

