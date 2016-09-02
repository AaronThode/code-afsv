%%%%%all_processing.m%%%%%
%  Conducts a bulk run of automated call detector
%   September 18, 2008 working version...
%   April 8, 2010 Revised and commented

function all_processing
clear ;fclose('all');close all
!rm *scr S* out*txt
%%Confirm that this script is being run from the correct directory...

if isempty(findstr(pwd,'Bulk'))
    error('Go to correct directory');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%CHANGE HERE...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Run options..

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
%                   in the final output file. Prevents possible
%                   writeovers while allowing individual DASAR files to be
%                   resused.
%               See'TOC_runs.m' for further details.
%               Example:  Icase_str='Shell08_Site%i_airgun.morph.Final';

%               alternate description:
%               a keyword string used to select time and spatial subset of deployment...
%               First field determines output folder and is used by TOC_runs to select
%                   appropriate site and company, which it then uses to select dates and
%               locations.  For example, if first field has 'Shell08' and 'Site5'
%                   strings, then that site is selected.
%               Second field determins subfolder
%               Third field determines what string in final filename
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
Site_vec=5;  %Vector of site numbers to process; Meant to be overridden

param=[];
run_options.debug.ctimes_debug=[];run_options.debug.tol=[];
[Icase_str,date_str_local,Site_vec,params_chc,param,run_options.debug.ctimes_debug,run_options.debug.tol,DASAR_str_local]=load_local_runparams(mfilename);

%Unpack detection parameters from keyword, store in 'param' structure
param=TOC_params(params_chc,param);
%%Uncomment to permit processing of a portion of a file.
%% Default is to process entire file
%param.energy.nstart=10*60*60*1000+0*60*1000;  %Define sample number to
%       start processing.
%param.energy.nsamples=0.25*60*60*1000;  %Number of samples to process.




run_options.debug.raw=0;    %If one, show debug energy detector information
run_options.debug.interval=0;%If one, show debug interval selector information.  If two, plot summary ICI,bearing, and frequency information
run_options.debug.morph=0;  %If one, show debug output for morph processing, and load debug_ctimes to learn when "failure" occures
run_options.debug.snips=0;  %If one, compare snips data with direct download of raw data...
run_options.debug.cross_channel=0;  %If 1, plot intermediate cross-channel output, if 2 plot more detailed input...
run_options.min_stations=2;  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria


run_options.Ncalls_to_sample=100;  %Number of snip samples to read into memory at once
run_options.force_energy_detector_run=1;  %If 1, always force JAVA cell CFAR to run
run_options.plot_manual_data=0;  %If one, plot spectrograms of manual detections along with logged frequency, duration, and type data
%param.energy.threshold=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%NEVER TOUCH BELOW UNDER NORMAL CIRCUMSTANCES%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for Isite=Site_vec
    Icase=sprintf(Icase_str,Isite);
    homedir=pwd;
    
    %%Load all pathnames used, including raw data directories, output
    %%directories, output directory for JAVA files.
    [rawdatadir,Icase,outputdir,param,manualdir]=load_pathnames(Icase,param);
    
    %Load default dates and DASAR limits (bad DASAR list)
    [~,date_str,DASAR_str,keyword]=TOC_runs(Icase);
    if ~isempty(date_str_local)
        date_str=date_str_local;
    end
    if ~isempty(DASAR_str_local)
        DASAR_str=DASAR_str_local;
    end
    
    run_options.filtering_stage=keyword.algorithm;  %'contour','morph','both','none': Determines what third stage processing takes place..
    
    
    
    %%Expand date parameters into a string of contiguous dates.
    date_range=expand_date_range(date_str);
    
    cd(homedir)
    
    param.energy.debug=0;
    param.energy.Fs=param.Fs;
    
    best_call_count=zeros(size(date_range,1),10,8);
    if exist('best_calls','var') 
        clear best_calls
    end
    %%Loop through desired dates
    for Idate=1:size(date_range,1)
        %%%Clean up workspace if we are forcing a first-stage CFAR JAVA run...
        !rm *scr *snips *detsum
        cd(homedir);
        mydir=pwd;
        
        fclose('all');
          
        
        try
            [goodDASAR,goodFile,goodName,goodDASARstr]=find_DASAR_dates(date_range(Idate,:),Isite,DASAR_str,rawdatadir,Icase);
            
        catch
            disp('Loading goodFile from processed data output');
            [goodDASAR,goodFile,goodName,goodDASARstr]=find_DASAR_dates_Processed(date_range(Idate,:),Isite,DASAR_str,outputdir,Icase);
            
        end
        
        
        if isempty(goodDASAR)
            disp(sprintf('%s was not assigned any DASARS at Site %i',date_range(Idate,:),Isite));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  load manual detection data in question... %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        searchstr=dir([manualdir '/*' date_range(Idate,5:end) 'S' int2str(Isite) '*.tsv']);
        
        
        if length(searchstr)==1
            %manualdata=load([manualdir '/' searchstr(1).name]);
            [individual_stored]=read_tsv([manualdir '/' searchstr(1).name],0,Inf,goodDASAR);
            manual=individual_stored;
            manual.tabs=manual.ctime;
            Igood=find(manual.ctime>0);
            manual.tabs(Igood)=datenum(1970,1,1,0,0,manual.ctime(Igood));
            
        else
            disp(sprintf('More than one result for %s',searchstr));
        end
        %%%Process entire algorithm, CFAR, interval removal, image
        %%% processing...
        for Istation=1:length(goodName)  %For each DASAR desired
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Process a single DASAR day--all three stages%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                best_calls{Idate,Istation}=process_one_unit_CollectManualSnippets(goodFile{Istation},goodName{Istation},  param,run_options,manual,Istation);
                best_calls{Idate,Istation}.site=Isite;
                
                for Iwctype=1:8  %Call types
                    best_call_count(Idate,Istation,Iwctype)=sum(best_calls{Idate,Istation}.wctype==Iwctype);
                end
            catch
                disp(sprintf('%s failed',goodName{Istation}));
                keyboard
            end
            
        end %goodDASAR
        
        save temp best_calls best_call_count
        
    end %IDate
    %%Output name includes 'morph', historical residue...
    fname_out=sprintf('CallSnips_%s_%s_Site%i',date_range(1,:),date_range(Idate,:),Isite);
    save(fname_out,'best_calls','best_call_count','goodDASAR');
    
end %Isite



end
