
%%%%%%%%%%load_program_location.m%%%%%%%
%  Aaron Thode
%  Provide location of JAVA program on multiple platforms...

function program_dir=load_program_location

[~,local_machine]=unix('hostname');
local_machine=deblank(local_machine);

if strcmp(local_machine,'macmussel.ucsd.edu')
    program_dir ='/Users/Shared/DataPreProcessor/LatestVersion/PreProcessV1.jar';
elseif strcmp(getenv('USER'),'fielduser')
    %program_dir ='/Users/thode/DataPreProcessor/PreProcessV1.jar';
    program_dir ='~/Desktop/DataPreProcessor/PreProcessV1.jar';
elseif ~isempty(strfind(local_machine,'thode-lt'))
    program_dir ='/Users/thode/Documents/workspace/DataPreProcessor/PreProcessV1.jar';
    
else
    %program_dir='/Applications/AllFile_specgram_viewer_Mac/application/PreProcessV1.jar';
    
    %if exist(program_dir,'file')~=2
    program_dir='~/Desktop/Ulysses/PreProcessV1.jar';
    if exist(program_dir,'file')~=2
        uiwait(warndlg('Please copy a file named "PreProcessV1.jar" into a folder named "Ulysses" onto your Desktop','JAVA program not found!','modal'));
    end
    %end
    
    
end


