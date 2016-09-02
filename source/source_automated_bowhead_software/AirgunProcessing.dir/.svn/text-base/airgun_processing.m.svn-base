%%%%%airgun_processing.m%%%%%
%  Conducts a bulk run of automated call detector to Greeneridge Data
%  to locate signals that pass ICI test...

%   September 18, 2008 working version...
%   December 27, 2009 revision
%   April 8, 2010 Revised and commented


clear ;fclose('all');
!rm *detsum *snips
path(path,pwd);
%path('../BulkProcessingScripts2009',path);
%path('../CommonScripts.dir',path);
%path('../ComputerSpecificScripts.dir',path);

%%Confirm that this script is being run from the correct directory...
if isempty(findstr(pwd,'Airgun'))
    error('Go to correct directory');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% CHANGE HERE...%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%Load scenario%%%
%  The following four variables define the scenario to be run:
%   Icase_str:  A string in the form 'NAMEYY_Site%i_XXX.morph.TAG', to be
%               used by TOC_runs.m.
%               The string defines the output folder structure and file
%               names.  The text before the period defines the output
%               folder name, and the text after the final period defines
%               the text string tag used in the '*mat' file output.
%               NAME is either 'Shell' or 'BP', i.e. the owner of the
%               deployment.
%               -YY is the two-digit year.
%               -XX is a descriptive label to be included in the folder
%                name.  A special case is when XX contains the string
%                'airgun'.  Then TOC_runs.m uses parameters associated with
%                airgun processing.  Another special case is the string
%                'Core2', which signals load_pathnames.m to use alternate
%                directories for raw data and output files (permitting
%                multi-core processing using separate disk drives).
%               -TAG is a text label that is embedded into the final MATLAB
%                   *mat file.  The phrase '_airguns.mat' is also included
%                   in the final output file
%               See'TOC_runs.m' for further details.
%               Example:  Icase_str='Shell08_Site%i_airgun.morph.Final';
%  date_str:  Two element cell array with strings of year/month/day; e.g. '20090902'
%               is 2009 Sept. 2.  A set of multiple entries, e.g.
%               date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};
%               defines a set of non-contiguous dates
%  Site_vec:   vector of indicies between 1 and 5, representing a site
%                number.  BP Northstar is always Site 1.
%  params_chc: A text string that is used by TOC_params or
%               TOC_params_airgun to download the run parameters of the
%               automated algorithms
%               Example: params_chc='Shell08_airgun';

%%load_local_runparams.m is a m-file stored in the local directory that
%%will specify the run parameters.  Storing this m-file locally permits
% different dates, labels, etc to be defined on different machines.
param=[];


[Icase_str,date_str_local,Site_vec,params_chc,param, ...
    run_options.debug.ctimes_debug,run_options.debug.tol, ...
    DASAR_str_local, ...
    rawdatadir_local,outputdir_local,locationdir_local]=load_local_runparams(mfilename);

if ~isempty(run_options.debug.ctimes_debug),
    disp('debug times loaded');
end
%Unpack detection parameters from keyword, store in 'param' structure
param=TOC_params_airgun(params_chc,param);
%param.energy.nsamples=1*60*60*1000;

%%Run parameters.  'run_options' stores information that influence the
%   control flow of the program, such as the level of debug information to
%   show, whether to run all stages of the program, etc.  Differs from
%   'param' in that the run_options values generally do not concern
%   algorithm parameters/adjustments, or values that do not change
%   frequently.

run_options.debug.raw=0;    %If one, show debug energy detector information
run_options.debug.interval=0;%If one, show debug interval selector information.  If two, plot summary ICI,bearing, and frequency information
run_options.debug.snips=0;  %If one, compare snips data with direct download of raw data...
run_options.Ncalls_to_sample=20;  %Number of snip samples to read into memory at once
run_options.force_energy_detector_run=1;  %If 1, always force JAVA cell CFAR to run
run_options.load_single_DASAR_results=0;  %If one, load stage 3 output (morph processing) from file
%run_options.max_clip_count=5;  %Number of clippings before signal rejected.  Commented out because now I just
%       output number of clipped samples
run_options.bearing_alg='sel_ratio';  %Which algorithm to use for computing DASAR bearings.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%NEVER TOUCH BELOW UNDER NORMAL CIRCUMSTANCES%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Isite=Site_vec  %For each site to process...
    
    Icase=sprintf(Icase_str,Isite);  %Create a specific Icase string for each Site
    homedir=pwd;
     
    %%Load all pathnames used, including raw data directories, output
    %%directories, output directory for JAVA files.
    [rawdatadir,Icase,outputdir,param,junk,locationdir]=load_pathnames(Icase,param);
    if ~isempty(rawdatadir_local)
        rawdatadir=rawdatadir_local;
    end
    if ~isempty(outputdir_local)
        outputdir=outputdir_local;
    end
    if ~isempty(locationdir_local)
        locationdir=locationdir_local;
    end
    
    %If load_local_runparam has defined variables, override default values...
    [~,date_str,DASAR_str,keyword]=TOC_runs(Icase);
    if ~isempty(date_str_local)
        date_str=date_str_local;
    end
    
    if ~isempty(DASAR_str_local)
        if ~iscell(DASAR_str_local)
            DASAR_str=DASAR_str_local;
        else
            DASAR_str=DASAR_str_local{Isite};
        end
    end
    
    
    
    %%Expand date parameters into a string array of contiguous dates.
    date_range=expand_date_range(date_str);
    
    cd(homedir)
    
    %%Construct output folder hierarchy
    finaldir2=sprintf('%s/Site_%02i/',outputdir,Isite);
    finaldir1=sprintf('%s/Site_%02i/%s',outputdir,Isite,keyword.scenario);
    finaldir=sprintf('%s/Site_%02i/%s/%s',outputdir,Isite,keyword.scenario,keyword.algorithm);
    disp(sprintf('Output to be written to %s',finaldir));
    param.calibration_keyword=Icase;
    
    
    param.energy.debug=0;
    param.energy.Fs=param.Fs;
    eval(sprintf('!mkdir %s',finaldir2));
    eval(sprintf('!mkdir %s',finaldir1));
    eval(sprintf('!mkdir %s',finaldir));
    
    %%Loop through desired dates
    for Idate=1:size(date_range,1)
         %%%Clean up workspace if we are forcing a first-stage CFAR JAVA run...
        
        cd(homedir);
        mydir=pwd;
        %%For each Site and date, identify all DASARS operational on that date,
        %%      potentially restricted by DASAR_str.
        %       goodDASAR: a cell string array containing the names of each
        %       DASAR to process...
        [goodDASAR,goodFile,goodName,goodDASARstr]=find_DASAR_dates(date_range(Idate,:),Isite,DASAR_str,rawdatadir,Icase);
        if isempty(goodDASAR)
            disp(sprintf('%s was not assigned any DASARS at Site %i',date_range(Idate,:),Isite));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        %Process a single DASAR day
        %%%%%%%%%%%%%%%%%%%%%
        
        for Istation=1:length(goodName)  %For each DASAR desired
	    !rm *detsum *snips
            !rm *.scr
            try
                disp(sprintf('Processing %s',goodName{Istation}));
                
                %%Option to load specific times to activiate debug
                %%  commands..
                
                
                [raw_ctimes,interval_ctimes,airgun_shots]=process_one_unit_airgun(goodFile{Istation},goodName{Istation},param,run_options);
             
                if run_options.force_energy_detector_run==1
                    cd(param.energy.dir_out)
                    !rm *.snips *.detsum
                end
                cd(finaldir);
                fname_out=sprintf('%s_%s_airguns.mat',goodName{Istation},keyword.stage);
                
                save(fname_out,'airgun_shots','raw_ctimes','interval_ctimes','param');
                %Write time series files in text format.  FIle name will be
                %the DASAR name...
                
                %write_airgun_textfile([finaldir '/' goodName{Istation}],airgun_shots);
                write_airgun_textfile([goodName{Istation}],airgun_shots);
                
                cd(mydir);
                
            catch  %%In case of crash, explain why...
                disp(sprintf('Crash! couldn''t process %s',goodName{Istation}));
                errormsg=lasterror;
                disp(errormsg.message);
                %disp(length(errormsg.stack))
                if (length(errormsg.stack)>0)
                    disp(errormsg.stack(1).file)
                    disp(errormsg.stack(1).name)
                    disp(errormsg.stack(1).line)
                end
                keyboard;
            end
        end %Istation
        
    end %IDate
end %Isite


