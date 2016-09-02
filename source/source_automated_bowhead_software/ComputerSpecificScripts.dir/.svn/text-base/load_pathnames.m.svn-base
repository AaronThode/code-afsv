%function [rawdatadir,Icase,outputdir,param,manualdir]=load_pathnames(Icase,param)
%  Given a descriptive string, check which computer is being used and
%  assign correct pathnames to raw data files, output data files, etc.
%  Also appends useful directories to the MATLAB path
%  Inputs:
%       Icase: A string in the form 'NAMEYY_Site%i_XXX.morph.TAG'.
%              If the string contains the phrase 'Core2' then an alternate
%               set of pathnames can be assigned, to permit multi-core
%               processing on separate external hard drives.
%       param: a structure of variables.  The field
%                param.energy.dir_out defines the output directory of
%               JAVA CFAR program, which can be adjusted here if desired.
%  Outputs:
%       rawdatadir:  Base directory of acoustic data files
%       Icase:  Same as Icase, but with the phrase 'Core2' removed.
%       outputdir: Base directory for processed output, typically starts
%           with 'Processed'.
%       param: Same as input param, but with param.energy.dir_out
%           optionally modified.
%       manualdir: Location of manually-processed TSV files for use in
%           algorithm training.
function [rawdatadir,Icase,outputdir,param,manualdir,locationdir]=load_pathnames(Icase,param)


%%%First check if external hard drive is needed for specific datasets

[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);
rawdatadir=[];outputdir=[];manualdir=[];
locationdir='../DASARlocations';

%%%Set path name for File_IO and Java_IO scripts
if strcmp(local_machine,'macmussel.ucsd.edu')
    path(path,'/Users/Shared/MATLAB/File_IO');
    path(path,'/Users/Shared/MATLAB/Java_IO');
    path(path,'/Users/Shared/MATLAB/AllFile_specgram_viewer');
    path(path,'/Users/Shared/MATLAB/UsefulScripts.dir');
elseif ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
    path(path,'/Users/khkim/Thode/MATLAB/File_IO');
    path(path,'/Users/khkim/Thode/MATLAB/Java_IO');
    path(path,'/Users/khkim/Thode/MATLAB/AllFile_specgram_viewer');
    path(path,'/Users/khkim/Thode/MATLAB/UsefulScripts.dir');
    
    
elseif ~isempty(findstr(local_machine,'thode')) %Thode laptop
    %path(path,'/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/MATLAB/File_IO');
    %path(path,'/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/MATLAB/Java_IO');
    path(path,'/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/MATLAB/AllFile_specgram_viewer');
    path(path,'/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/MATLAB/UsefulScripts.dir');
end

if findstr(Icase,'BP10'),
    if ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/BP_2010_Disk2/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
        else
            rawdatadir='/Volumes/BP_2010/';  %Location root of "sio" binary data files
        end
        outputdir='/Users/khkim/Thode/2010_Arctic_Analysis/Processed2010';  %Location of output detections
        loc.base='/Users/khkim/Thode/2010_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/khkim/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    elseif strcmp(local_machine,'macmussel.ucsd.edu')
        
        rawdatadir='/Volumes/macmussel2/Arctic_2009/';  %Location root of "sio" binary data files
        
        outputdir='/Volumes/macmussel1/Arctic_2009/Processed';  %Location of output detections
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end

if findstr(Icase,'BP11'),
    if ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/BP_2011_Disk2/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
        else
            rawdatadir='/Volumes/BP_2011/';  %Location root of "sio" binary data files
        end
        outputdir='/Users/khkim/Thode/2011_Arctic_Analysis/Processed2011';  %Location of output detections
        loc.base='/Users/khkim/Thode/2011_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/khkim/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    elseif strcmp(local_machine,'macmussel.ucsd.edu')
        
        rawdatadir='/Volumes/macmussel2/Arctic_2009/';  %Location root of "sio" binary data files
        
        outputdir='/Volumes/macmussel1/Arctic_2009/Processed';  %Location of output detections
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end

if findstr(Icase,'BP12'),
    if ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/2012DATAB/BP_2012/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
        else
            rawdatadir='/Volumes/2012DATAB/BP_2012/';  %Location root of "sio" binary data files
        end
        outputdir='/Users/khkim/Thode/2012_Arctic_Analysis/Processed.dir';  %Location of output detections
        loc.base='/Users/khkim/Thode/2012_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/khkim/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    elseif strcmp(local_machine,'macmussel.ucsd.edu')
        
        rawdatadir='/Volumes/macmussel2/Arctic_2009/';  %Location root of "sio" binary data files
        
        outputdir='/Volumes/macmussel1/Arctic_2009/Processed';  %Location of output detections
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end

if findstr(Icase,'BP13'),
    if ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/BP_2013/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
        else
            rawdatadir='/Volumes/BP_2013/';  %Location root of "sio" binary data files
        end
        outputdir='/Users/khkim/Thode/2013_Arctic_Analysis/Processed.dir';  %Location of output detections
        loc.base='/Users/khkim/Thode/2013_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/khkim/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    elseif strcmp(local_machine,'macmussel.ucsd.edu')
        
        rawdatadir='/Volumes/macmussel2/Arctic_2009/';  %Location root of "sio" binary data files
        
        outputdir='/Volumes/macmussel1/Arctic_2009/Processed';  %Location of output detections
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        manualdir=[];
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end

if findstr(Icase,'AURAL11')
%     if strcmp(local_machine,'macmussel.ucsd.edu')
%         %rawdatadir='/Volumes/Shell2011_GSI_Data/';  %Location root of "sio" binary data files
%         rawdatadir='/Volumes/Shell2011_GSI_Data_Copy2/';
%         outputdir=sprintf('/Volumes/macmussel1/Arctic_2011/Processed2011');  %Location of output detections
%         loc.base='/Users/Shared/Projects/Automated_Arctic_Software';
%         manualdir='/Users/Shared/Projects/2011_Arctic_Analysis/TSV_files_Shell11/Shell11_ManualTSVFiles';
%         param.energy.dir_out='.';
%         %rundir=[loc.base '/BulkProcessingScripts2009'];
%         
%     elseif ~isempty(findstr(local_machine,'thode'))
%         loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror';
%         rawdatadir=[loc.base '/RawData'];  %Location root of "sio" binary data files
%         rawdatadir='/Volumes/Shell2010_GSI_Data_Disk3/';
%         outputdir=['/Volumes/ThodePortable2/2010_Arctic_Analysis/Processed2010/'];  %Location of output detections
%         manualdir=[loc.base '/2010_Arctic_Analysis/TSV_files_Shell10/Shell10_ManualTSVFiles'];
%         param.energy.dir_out='.';
%         %rundir=[loc.base '/BulkProcessingScripts2009'];
%         
    if ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/Shell2010_GSI_Data_Disk2/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
        else
            rawdatadir='/Users/khkim/Thode/2012_AURAL_OW/Data_like_OW_11/';  %Location root of "sio" binary data files
        end
        fprintf('Raw data directory is %s\n',rawdatadir);
        outputdir='/Users/khkim/Thode/2012_AURAL_OW/Processed2012';  %Location of output detections
        loc.base='/Users/khkim/Thode/2012_AURAL_OW';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        manualdir=[loc.base '/BulkProcessing.dir/TSV_files.dir'];
        param.Fs=32768;
        
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end
if findstr(Icase,'Shell11')
    if strcmp(local_machine,'macmussel.ucsd.edu')
        %rawdatadir='/Volumes/Shell2011_GSI_Data/';  %Location root of "sio" binary data files
        rawdatadir='/Volumes/Shell2011_GSI_Data_Copy2/';
        outputdir=sprintf('/Volumes/macmussel1/Arctic_2011/Processed2011');  %Location of output detections
        loc.base='/Users/Shared/Projects/Automated_Arctic_Software';
        manualdir='/Users/Shared/Projects/2011_Arctic_Analysis/TSV_files_Shell11/Shell11_ManualTSVFiles';
        param.energy.dir_out='.';
        %rundir=[loc.base '/BulkProcessingScripts2009'];
        
    elseif ~isempty(findstr(local_machine,'thode'))
        loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror';
        rawdatadir=[loc.base '/RawData'];  %Location root of "sio" binary data files
        rawdatadir='/Volumes/Shell2010_GSI_Data_Disk3/';
        outputdir=['/Volumes/ThodePortable2/2010_Arctic_Analysis/Processed2010/'];  %Location of output detections
        manualdir=[loc.base '/2010_Arctic_Analysis/TSV_files_Shell10/Shell10_ManualTSVFiles'];
        param.energy.dir_out='.';
        %rundir=[loc.base '/BulkProcessingScripts2009'];
        
    elseif ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/Shell2010_GSI_Data_Disk2/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2
        else
            rawdatadir='/Volumes/Shell2011_GSI_Data/';  %Location root of "sio" binary data files
        end
        fprintf('Raw data directory is %s\n',rawdatadir);
        outputdir='/Users/khkim/Thode/2011_Arctic_Analysis/Processed2011';  %Location of output detections
        loc.base='/Users/khkim/Thode/2011_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        manualdir=[loc.base '/TSV_files_Shell10//Shell10_ManualTSVFiles'];
        
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end


%%Shell 2010 data
if findstr(Icase,'Shell10'),
    if strcmp(local_machine,'macmussel.ucsd.edu')
        rawdatadir='/Volumes/Shell2010_GSI_Data_Disk3';  %Location root of "sio" binary data files
        outputdir=sprintf('/Volumes/macmussel1/Arctic_2010/Processed');  %Location of output detections
        loc.base='/Users/Shared/Projects/Automated_Arctic_Software';
        manualdir='/Users/Shared/Projects/2010_Arctic_Analysis/TSV_files_Shell10/Shell10_ManualTSVFiles';
        param.energy.dir_out='.';
        %rundir=[loc.base '/BulkProcessingScripts2009'];
        
    elseif ~isempty(findstr(local_machine,'thode'))  %%If a laptop, use external hard drive
        loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror';
        rawdatadir='/Volumes/ThodePortable2/2010_Beaufort_Shell_DASAR/';
        outputdir=['/Volumes/ThodePortable2/2010_Beaufort_Shell_DASAR/Processed2010/'];  %Location of output detections
        manualdir=[loc.base '/2010_Arctic_Analysis/TSV_files_Shell10/Shell10_ManualTSVFiles'];
        param.energy.dir_out='.';
        %rundir=[loc.base '/BulkProcessingScripts2009'];
    elseif ~isempty(findstr(local_machine,'dgrebner'))
        loc.base='/Users/dgrebner/Projects';
        rawdatadir=[loc.base '/ArcticProcessing_AudioData'];  %Location root of "sio" binary data files
        outputdir=[loc.base '/BulkProcessing/Processed_2009'];  %Location of output detections
        manualdir=[loc.base '/2009_Arctic_Analysis/TSV_files_Shell08'];
        
    elseif ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/Shell2010_GSI_Data_Disk2/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2
        else
            rawdatadir='/Volumes/Shell2010_GSI_Data/';  %Location root of "sio" binary data files
        end
        fprintf('Raw data directory is %s\n',rawdatadir);
        outputdir='/Users/khkim/Thode/2010_Arctic_Analysis/Processed2010';  %Location of output detections
        loc.base='/Users/khkim/Thode/2010_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        manualdir=[loc.base '/TSV_files_Shell10//Shell10_ManualTSVFiles'];
        
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end

if findstr(Icase,'Shell09'),
    if strcmp(local_machine,'macmussel.ucsd.edu')
        rawdatadir='/Volumes/Shell2009_GSI_Data/';  %Location root of "sio" binary data files
        outputdir=sprintf('/Volumes/macmussel1/Arctic_2009/Processed');  %Location of output detections
        loc.base='/Users/Shared/Projects/Automated_Arctic_Software';
        manualdir=['/Users/Shared/Projects/2009_Arctic_Analysis/TSV_files_Shell09/Shell09_ManualTSVFiles_Final_First12Hours'];
        param.energy.dir_out='.';
        %rundir=[loc.base '/BulkProcessingScripts2009'];
        
    elseif ~isempty(findstr(local_machine,'thode'))
        loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror';
        rawdatadir=[loc.base '/RawData'];  %Location root of "sio" binary data files
        outputdir=['/Volumes/ThodePortable2/2009_Arctic_Analysis/Processed2009'];  %Location of output detections
        manualdir=[loc.base '/2009_Arctic_Analysis/TSV_files_Shell09/Shell09_ManualTSVFiles_Final_First12Hours'];
        param.energy.dir_out='.';
        
        %rundir=[loc.base '/BulkProcessingScripts2009'];
    elseif ~isempty(findstr(local_machine,'dgrebner'))
        loc.base='/Users/dgrebner/Projects';
        rawdatadir=[loc.base '/ArcticProcessing_AudioData'];  %Location root of "sio" binary data files
        outputdir=[loc.base '/BulkProcessing/Processed_2009'];  %Location of output detections
        manualdir=[loc.base '/2009_Arctic_Analysis/TSV_files_Shell08'];
        
    elseif ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            rawdatadir='/Volumes/Shell2009_GSI_Data_Disk2/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
        else
            rawdatadir='/Volumes/Shell2009_GSI_Data/';  %Location root of "sio" binary data files
        end
        outputdir='/Users/khkim/Thode/2009_Arctic_Analysis/Processed2009';  %Location of output detections
        loc.base='/Users/khkim/Thode/2009_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/khkim/Arctic_2008/CommonScripts.dir');
        manualdir=[loc.base '/TSV_files_Shell09//Shell09_ManualTSVFiles_Final_First12Hours'];
        
    end
    %path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
    
    return
end

if findstr(Icase,'Shell07'),
    if strcmp(local_machine,'macmussel.ucsd.edu')
         Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            %error('No alternate disk available for 2007 data.')
            rawdatadir='/Volumes/2007GSI_copy/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
            disp(sprintf('Icase is now %s',Icase));
        else
            rawdatadir='/Volumes/2007GSI_copy/';  %Location root of "sio" binary data files
        end
        outputdir='/Volumes/macmussel1/Arctic_2007/Processed2007.dir';  %Location of output detections
        loc.base='/Users/Shared/Projects/2007_Arctic_Analysis';
        param.energy.dir_out='.';
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        manualdir=[loc.base '/TSV_files_Shell07'];
        rundir=[loc.base '/BulkProcessingScripts'];
        param.calibration_dir='/Volumes/2007GSI_copy';
        
    elseif ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        Iseg=findstr(Icase,'Core2')-1;
        if ~isempty(Iseg)
            %error('No alternate disk available for 2007 data.')
            rawdatadir='/Volumes/2007GSI_copy/';  %Location root of "sio" binary data files
            Icase2=[Icase(1:Iseg) Icase((Iseg+length('Core2')+1):end)];
            Icase=Icase2;
            disp(sprintf('Icase is now %s',Icase));
        else
            rawdatadir='/Volumes/2007GSI/';  %Location root of "sio" binary data files
        end
        outputdir='/Users/khkim/Thode/2007_Arctic_Analysis/Processed2007';  %Location of output detections
        loc.base='/Users/khkim/Thode/2007_Arctic_Analysis';
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        manualdir=[loc.base '/TSV_files_Shell07'];
        rundir=[loc.base '/BulkProcessingScripts'];
        param.calibration_dir='/Volumes/2007GSI';
        
    elseif ~isempty(findstr(local_machine,'thode'))
        loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror';
        rawdatadir=[loc.base '/RawData'];  %Location root of "sio" binary data files
        %outputdir=[loc.base '/2007_Arctic_Analysis/Processed_2007'];  %Location of output detections
        rundir=[loc.base '/2007_Arctic_Analysis/BulkProcessingScripts'];
        param.energy.dir_out='.';
        param.calibration_dir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/RawData';
        
        outputdir=['/Volumes/ThodePortable2/2007_Arctic_Analysis/Processed2007'];  %Location of output detections
        outputdir='.';
        manualdir=[loc.base '/2007_Arctic_Analysis/TSV_files_Shell07'];
        param.energy.dir_out='.';
        
         param.energy.dir_out='.';
        
        
    end
    path(path,rundir);
    %path(path,[rundir '/Localization_mfiles']);
end

if findstr(Icase,'Shell08'),
    if strcmp(local_machine,'macmussel.ucsd.edu')
        rawdatadir='/Volumes/macmussel2/Arctic_2008';  %Location root of "sio" binary data files
        outputdir=sprintf('/Volumes/macmussel1/Arctic_2008/Processed_2008');  %Location of output detections
        loc.base='/Users/Shared/Projects/2008_Arctic_Analysis';
        manualdir=[loc.base '/TSV_files_Shell08'];
        param.energy.dir_out='.';
        
        
    elseif ~isempty(findstr(local_machine,'thode'))
        loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror';
        rawdatadir=[loc.base '/RawData'];  %Location root of "sio" binary data files
        outputdir=[loc.base '/Processed_2008'];  %Location of output detections
        manualdir=[loc.base '/2009_Arctic_Analysis/TSV_files_Shell08'];
        param.energy.dir_out='.';
        
    elseif ~isempty(findstr(local_machine,'dgrebner'))
        loc.base='/Users/dgrebner/Projects';
        rawdatadir=[loc.base '/ArcticProcessing_AudioData'];  %Location root of "sio" binary data files
        outputdir=[loc.base '/BulkProcessing/Processed'];  %Location of output detections
        manualdir=[loc.base '/2009_Arctic_Analysis/TSV_files_Shell08'];
        param.energy.dir_out='.';
    elseif ~isempty(findstr(local_machine,'KatsMacPro'))  %KatMacGSI
        
        loc.base='/Users/khkim/Thode/2008_Arctic_Analysis';
        if isempty(findstr(Icase,'Core2'))  %Processor one
            rawdatadir='/Volumes/Shell2008_GSI_Data';  %Location root of "sio" binary data files
            rundir=[loc.base '/BulkProcessingScripts'];
        else
            rawdatadir='/Volumes/Shell2008_Disk3';  %Location root of "sio" binary data files
            rundir=[loc.base '/BulkProcessingScriptsCore2'];
        end
        param.energy.dir_out=[rawdatadir '/EnergyFiles.dir'];
        %path(path,'/Users/thode/Arctic_2008/CommonScripts.dir');
        outputdir=[loc.base '/Processed'];  %Location of output detections
        manualdir=[loc.base '/TSV_files_Shell08'];
    end
    % rundir=[loc.base '/BulkProcessingScripts2009'];
    
    
    return
end
