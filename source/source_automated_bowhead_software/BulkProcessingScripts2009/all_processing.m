%%%%%all_processing.m%%%%%
%  Conducts a bulk run of automated call detector
%   September 18, 2008 working version...
%   April 8, 2010 Revised and commented

function all_processing
clear ;fclose('all');close all

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


param=[];
run_options.debug.ctimes_debug=[];run_options.debug.tol=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%load_local_runparams.m is a m-file stored in the local directory that
%%will specify the run parameters.  Storing this m-file locally permits
% different dates, DASAR clusters, data directories, labels, etc to be defined on different machines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Icase_str,date_str_local,Site_vec,params_chc,param, ...
    run_options.debug.ctimes_debug,run_options.debug.tol, ...
    DASAR_str_local,rawdatadir_local,outputdir_local,locationdir_local]=load_local_runparams(mfilename);

%Unpack detection parameters from keyword, store in 'param' structure
param=TOC_params(params_chc,param);
%%Uncomment to permit processing of a portion of a file.
%% Default is to process entire file
%param.energy.nstart=10*60*60*1000+0*60*1000;  %Define sample number to
%       start processing.
%param.energy.nsamples=0.25*60*60*1000;  %Number of samples to process.


%%Run parameters.  'run_options' stores information that influence the
%   control flow of the program, such as the level of debug information to
%   show, whether to run all stages of the program, etc.  Differs from
%   'param' in that the run_options values generally do not concern
%   algorithm parameters/adjustments, or values that do not change
%   frequently.

%Example settings for common situations:
%   To compute all stages using trained network:
%      load_neuralnet_output=0,apply_neural_net=1,load_single_DASAR_results=0;
%   To apply a different neural network to image-processed data:
%       load_neuralnet_output=0,apply_neural_net=1,load_single_DASAR_results=1;  Be sure to change 'TAG' label in Icasestr
%   To redo linking and localization, using a different output file name:
%       load_neuralnet_output=0,apply_neural_net=1,load_single_DASAR_results=1;
%   To redo_linking and localization, rewriting over previous results...note link_calls=apply_neural_net
%       load_neuralnet_output=1,apply_neural_net=1,load_single_DASAR_results=0/1(doesn't matter);

%If one, load ouptut of neural network (perform cross-linking and localization only, don't rerun
% earlier stages)
run_options.load_neuralnet_output=0;

%If one, apply a neural network stored in 'param.nnet'.  If zero, creates a
%training set.
run_options.apply_neural_net=1;

%If one, load stage 3 output (image processing) from file, don't rerun
%   earlier stages.  Note that load_neuralnet_output needs to be 0 to use
%   this option
run_options.load_single_DASAR_results=1;
run_options.find_DASAR_dates_rawdata=1;  %If one, use raw data directories to determine available DASARs.  Otherwise use processed directories

run_options.debug.raw=0;    %If one, show debug energy detector information
run_options.debug.interval=0;%If one, show debug interval selector information.  If two, plot summary ICI,bearing, and frequency information
run_options.debug.morph=0;  %If one, show debug output for morph processing, and load debug_ctimes to learn when "failure" occures
run_options.debug.snips=0;  %If one, compare snips data with direct download of raw data...
run_options.debug.cross_channel=0;  %If 1, plot intermediate cross-channel output, if 2 plot more detailed input...
run_options.min_stations=2;  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria


run_options.Ncalls_to_sample=100;  %Number of snip samples to read into memory at once
run_options.force_energy_detector_run=1;  %If 1, always force JAVA cell CFAR to run
run_options.max_time_delay=7;  %Maximum distance in seconds allowed between DASARS when deciding whether to link
run_options.dt_slop=5;  %How much time tolerance to give when matching detections between DASARs
run_options.bearing_alg='sel_ratio';
run_options.plot_locations=0;  %If one, plot locations..
run_options.localization_alg='Huber'; %HuberCalibratedKappa, repeated_median, MedianHuber
run_options.filter_chc='min_weight'; %dt_error, min_weight
run_options.kappa_Nsamples=100;
param.morph.threshold_chc='local_peaks';
run_options.auto_location_filter=1;  %Final desperation filter on location features to strip out seismics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%NEVER TOUCH BELOW UNDER NORMAL CIRCUMSTANCES%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_options.link_calls=run_options.apply_neural_net;  %If one, link calls together using cross_channel_match
    

for Isite=Site_vec
    Icase=sprintf(Icase_str,Isite);
    homedir=pwd;
    
    %%Load all pathnames used, including raw data directories, output
    %%directories, output directory for JAVA files.
    [rawdatadir,Icase,outputdir,param]=load_pathnames(Icase,param);
    if ~isempty(rawdatadir_local)
        rawdatadir=rawdatadir_local;
    end
    if ~isempty(outputdir_local)
        outputdir=outputdir_local;
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
            try
                DASAR_str=DASAR_str_local{Isite};
            catch exception
                if (strcmp(exception.identifier,'MATLAB:badsubscript'))
                   if Isite==0
                       disp('Isite==0 desired, will use DASAR_str_local{6}');
                       DASAR_str=DASAR_str_local{6};
                   end
                end

            end
        end
    end
    
    if ~isempty(locationdir_local)
        locationdir=locationdir_local;
    end
    
    run_options.filtering_stage=keyword.algorithm;  %'contour','morph','both','none': Determines what third stage processing takes place..
    
    
    
    %%Expand date parameters into a string of contiguous dates.
    date_range=expand_date_range(date_str);
    
    cd(homedir)
    
    param.energy.debug=0;
    param.energy.Fs=param.Fs;
    
    %%Construct output folder hierarchy
    create_output_directories;
    
    %%Loop through desired dates
    for Idate=1:size(date_range,1)
        %%%Clean up workspace if we are forcing a first-stage CFAR JAVA run...
        if run_options.force_energy_detector_run==1
            cd(param.energy.dir_out)
            !rm *.snips *.detsum
            cd(homedir)
        end
        cd(homedir);
        mydir=pwd;
        try
            fclose('all');clear station locations;
            %%For each Site and date, identify all DASARS operational on that date,
            %%      potentially restricted by DASAR_str.
            %       goodDASAR: a cell string array containing the names of each
            %       DASAR to process...
            %
            %  We run find_DASAR_dates twice; once to find all DASARs available, and then again, but
            %       restricted to DASARs we only want to localize..
            
            if run_options.find_DASAR_dates_rawdata==1
                [goodDASAR,goodFile,goodName]=find_DASAR_dates(date_range(Idate,:),Isite,'*',rawdatadir,Icase);
                [goodDASAR_loc,goodFile_loc,goodName_loc]=find_DASAR_dates(date_range(Idate,:),Isite,DASAR_str,rawdatadir,Icase);
                
            else
                disp('Loading goodFile from processed data output');
                [goodDASAR,goodFile,goodName]=find_DASAR_dates_Processed(date_range(Idate,:),Isite,'*',outputdir,Icase);
                [goodDASAR_loc,goodFile_loc,goodName_loc]=find_DASAR_dates_Processed(date_range(Idate,:),Isite,DASAR_str,outputdir,Icase);
                
            end
            
            
            if isempty(goodDASAR)
                disp(sprintf('%s was not assigned any DASARS at Site %i',date_range(Idate,:),Isite));
            end
            
            if exist('station'),clear station raw_station;end
            
            %%%Process entire algorithm, CFAR, interval removal, image
            %%% processing...
            if run_options.load_neuralnet_output==1
                %%%If branch for recomputing linking only...that is, use
                %%%results through neural network processing stage, don't recompute neural network or bearing
                %%%estimation.
                
                load_stations_for_linking;
            else
                for Istation=1:length(goodName)  %For each DASAR desired, including DASARS with bad particle velocity info
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Process a single DASAR day--all three stages%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !rm *scr *detsum *snips
                    create_one_station(run_options.load_single_DASAR_results);  %calls process_one_unit
                    %Neural net filtering...
                end %goodDASAR
                
                %%If a bulk run, filter with Nnet.
                %% If building a training set, this step should be skipped
                if run_options.apply_neural_net==1&&exist('station','var')
                    disp(sprintf('Running neural network...'));
                    if exist(param.net.dir,'dir')==7
                        station=filter_with_Nnet(station,param.net.name,param.net.dir,param.net.threshold);
                        %station=filter_with_nnet(station_in,param.net.name,param.net.dir,[thresholds1(Ithresh1) thresholds2(Ithresh2)] );
                        
                        
                        %%Add a bearing feature to surviving segments (or include in training set)
                        for K=1:length(station)
                            disp(sprintf('Processing bearings for surviving segments on station %i',K));
                            station_tmp(K)= compute_bearings_station(station(K),goodFile{K},param,run_options.bearing_alg,0,run_options.kappa_Nsamples);
                        end
                        station=station_tmp;
                        clear station_tmp
                        
                    else
                        fprintf('The network %s does not exist!\n',param.net.dir);
                        yes=input('Hit return to continue without nnet filtering...');
                    end
                end  %run_options.apply_neural_net
         
                
                
            end  %run_options.load_neuralnet_output
            
            %%If no data have been loaded, continue to next date...
            if ~exist('station','var')
                continue
            end
            
            if run_options.link_calls==0
                %%If no neural network has been applied to stations, create
                %%  a training file set.
                
                cd(finaldir);
                fname_out=sprintf('%s_%s',goodName{end},'NoNeuralNet');
                save(fname_out,'station','raw_station','param','goodName','goodFile','Icase','Isite');
                cd(mydir);
                
            else %If we have applied neural network output and wish to link and localize..
                %%cross-channel matching
                save temp
                
                %%This internal (nested) function creates the locations object--not passed explicitly because
                %%it can be a huge variable...
                linking_index=cross_link_stations;
                
                %%Do a brute-force feature filtering of linked locations...
                if run_options.auto_location_filter==1
                    
                    [Ipass]=crude_filter_feature_locations(locations,param.final_filter.names,param.final_filter.criteria,run_options.min_stations(1),1);
                    %Ifail=setdiff(1:length(auto.locations{1}),Ipass);
                    fprintf('Filter feature activated, %i out of %i locations pass\n',length(Ipass), length(locations));
                    
                    locations=locations(Ipass);
                    locations_ctime=locations_ctime(Ipass,:);
                end
                
                %%Extract bearing from raw data using detection times and durations
                %%for each location
                
                if ~isempty(locations)
                    %compute_bearings_for_locations;
                    compute_position_for_locations;
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%Write singleton detections (no localizations)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                try
                    for K=1:length(station)
                        disp(sprintf('Processing singletons on station %i',K));
                        if ~isempty(station(K).indicies)
                            
                            %%Remove detections that have been used in localizations
                            Isingleton=setdiff(1:length(station(K).ctime_min),linking_index.I(K,:));
                            
                            %%Only keep station detections that are NOT localized
                            tmp=trim_station(station(K),Isingleton);
                            
                            Ipass=crude_filter_feature_stations(tmp,param.final_filter.names,param.final_filter.criteria);
                            disp(sprintf('%i out of %i singletons in station %i pass',length(Ipass),length(Isingleton),K));
                            station_single(K)=trim_station(tmp,Ipass);
                           % station_single(K)=compute_bearings_station(tmp,goodFile{K},param,run_options.bearing_alg,0,run_options.kappa_Nsamples);
                        else
                            disp(sprintf('Empty station at %i',K));
                            station_single(K)=station(K);
                           % station_single(K)= compute_bearings_station(station(K),goodFile{K},param,run_options.bearing_alg,0,run_options.kappa_Nsamples);
                     
                        end
                        
                        
                    end
                    tsv_out_single=sprintf('%s_%s_singleton',keyword.stage,run_options.localization_alg);
                    [fname_tsv,Nwrite]=write_tsv_singleton(station,goodName,tsv_out_single);
                catch
                    fprintf('Crash! date_range(Idate,:) =%s singleton failure\n',date_range(Idate,:));
                    
                end
                
                %%Write final output localization file, but only if more
                %%than one station exists
                cd(finaldir);
                
                if length(station)>1
                    fname_out=sprintf('%s_%s_%s_FilteredLocations',goodName{end},keyword.stage,run_options.localization_alg);
                    save(fname_out,'locations','locations_ctime','station','station_single','raw_station', ...
                        'linking_index','param','goodName','goodFile','Icase','Isite','run_options');
                else
                    disp('Will not save a FilteredLocations file for a single station');
                end
                cd(mydir);
                
                
            end  %run_options.link_calls
            eval(sprintf('!rm %s/*snips %s/*detsum',param.energy.dir_out,param.energy.dir_out));
            
        catch  %%Processing this date failed...
            fprintf('Crash! date_range(Idate,:) =%s failed\n',date_range(Idate,:));
            errormsg=lasterror;
            disp(errormsg.message);
            %disp(length(errormsg.stack))
            if (length(errormsg.stack)>0)
                disp(errormsg.stack(1).file)
                disp(errormsg.stack(1).name)
                disp(errormsg.stack(1).line)
            end
            cd(mydir)
            
        end
        if exist('locations')
            clear locations
        end
    end %IDate
end %Isite


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Nested functions%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Construct output folder hierarchy
    function create_output_directories
        finaldir2=sprintf('%s/Site_0%i/',outputdir,Isite);
        finaldir1=sprintf('%s/Site_0%i/%s',outputdir,Isite,keyword.scenario);
        finaldir=sprintf('%s/Site_0%i/%s/%s',outputdir,Isite,keyword.scenario,keyword.algorithm);
        disp(sprintf('Output to be written to %s',finaldir));
        
        eval(sprintf('!mkdir %s',finaldir2));
        eval(sprintf('!mkdir %s',finaldir1));
        eval(sprintf('!mkdir %s',finaldir));
    end

%%% Load precomputed stations (data that have passed through energy
%%% detection, interval detection, and neural network)
    function load_stations_for_linking
        current_params=param;
        
        cd(finaldir);
        
        
        %Output file will have a keyword name, 'Huber', and
        %   'FilteredLocations'
        
        fname_out=sprintf('%s_%s_%s_FilteredLocations',goodName{end},keyword.stage,run_options.localization_alg);
        
        %%Sometimes when reprocessing data DASAR G (or highest letter
        %%DASAR) is excluded.
        %%  Therefore, Test what files are present and check the DASAR letter of the
        %%Filtered results...
        
        fname_test=dir(sprintf('*%s*FilteredLo*mat',date_range(Idate,:)));
        letter_check=fname_test(1).name(5);
        
        if ~strcmpi(DASAR_str(end),letter_check)
            Islash=min(findstr(fname_test(1).name,'_'))-1;
            fname_out=sprintf('%s_%s_%s_FilteredLocations',fname_test(1).name(1:Islash),keyword.stage,run_options.localization_alg);
        end
        
        
        disp(sprintf('loading %s',fname_out));
        past=load(fname_out);
        
        right=0;  %Comparison quality check to ensure we are using the right data
        %if strcmp(past.Icase,Icase)
        %    right=1;
        %end
        if Isite==past.Isite
            right=right+1;
        end
        
        try
            for I=1:length(goodName)
                for J=1:length(past.goodName)
                    if strcmp(goodName{I},past.goodName{J})
                        right=right+1;
                        if isfield(past,'station')&&~isempty(past.station)
                            station(I)=past.station(J);
                            raw_station(I)=past.raw_station(J);
                            continue
                        else
                            station(I)=[];
                            raw_station(I)=[];
                        end
                        
                        
                        
                    end
                end
            end
        catch
            disp(sprintf('Problem with matching desired stations with %s',fname_out));
        end
        
        if right~=length(goodName)+1
            disp('Not station match');
        end
        clear past
        
        
        param=current_params;
        cd(mydir);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Process a single DASAR day--all three stages%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function create_one_station(no_force_recompute)
        fprintf('Processing %s\n',goodName{Istation});
        
        %%Option to load specific times to activiate debug
        %%  commands...
        %[run_options.debug.ctimes_debug,run_options.debug.tol]= load_debug_ctimes(Istation);
        if ~isempty(run_options.debug.ctimes_debug),
            disp('debug times loaded');
        end
        %end
        
        %%Option to skip processing to load stored results--useful if
        %%only neural net or cross-matching calculations needed.
        if no_force_recompute==0
            [best_calls,debug_ctimes,debug_durations,raw_detections,interval_detections]=process_one_unit(goodFile{Istation},goodName{Istation},param,run_options);
            cd(finaldir);
            
            %%Output name includes 'morph', historical residue...
            fname_out=sprintf('%s_%s',goodName{Istation},run_options.filtering_stage);
            save(fname_out,'debug_ctimes','debug_durations','best_calls','raw_detections','interval_detections','param');
            
        else  %Load previous stages through image processing
            current_params=param;
            cd(finaldir);
            %%If names are obtained through processed results, don't need
            %%to add a _morph to name...
            if isempty(findstr(goodName{Istation},run_options.filtering_stage))
                fname_out=sprintf('%s_%s',goodName{Istation},run_options.filtering_stage);
            else
                fname_out=goodName{Istation};
            end
            stored_data=load(fname_out);
            param=current_params;
            best_calls=stored_data.best_calls;
            debug_ctimes=stored_data.debug_ctimes;
            debug_durations=stored_data.debug_durations;
            raw_detections=stored_data.raw_detections;
            interval_detections=stored_data.interval_detections;
            clear stored_data
            
        end
        cd(mydir);
        
        %if strcmp(run_options.filtering_stage,'none'),
        %    continue;
        %end
        
        %if ~isempty(best_calls.features)
        %%Extract features for future linking.
        try
            if ~isempty(best_calls.features)
                
                station(Istation)=create_station(best_calls.equalization_freq,best_calls.equalization,best_calls.features, ...
                    best_calls.labeled_image, ...
                    debug_ctimes,debug_durations,param.feature.names,param.feature.global_names, ...
                    param.feature.index ,param.feature.Nsegments);
                raw_station(Istation).raw_detections=raw_detections;
                raw_station(Istation).interval_detections=interval_detections;
                %raw_station(Istation).morph_detections.ctime=station(Istation).ctime_debug;
                %raw_station(Istation).morph_detections.duration=station(Istation).duration_debug;
                raw_station(Istation).morph_detections.ctime=station(Istation).ctime_min;
                raw_station(Istation).morph_detections.duration=station(Istation).Totalduration;
            else
                station(Istation).feature.ctime=[];
                station(Istation).indicies=[];
                station(Istation).feature.ctime=[];
                station(Istation).ctime_debug=[];
                station(Istation).duration_debug=[];
                station(Istation).SNR=[];
                station(Istation).SEL=[];
                
                raw_station(Istation).raw_detections=[];
                raw_station(Istation).interval_detections=[];
                raw_station(Istation).morph_detections.ctime=[];
                raw_station(Istation).morph_detections.duration=[];
            end
        catch
            disp('Station build failure, most likely no signals detected');
            %station(Istation)=[];
            
        end
        %         else
        %             station(Istation)=[];
        %             disp('Empty station being created ');
        %             raw_station(Istation)=[];
        %
        %         end
        clear best_calls raw_detections interval_detections
        
    end

%%Implicit function that creates locations...
    function linking_index=cross_link_stations
        if run_options.debug.cross_channel>0,
            run_options.debug.cross_channel=2;
            run_options.debug.J_anchor=1;
            run_options.debug.J_anchor_start=1;
            run_options.debug.Nfft=param.Nfft;
            run_options.debug.ovlap=param.ovlap;
            run_options.debug.Fs=param.Fs;
            run_options.debug.morph=param.morph;
            run_options.debug.merge=param.merge;
        end
        %Match calls across stations.  
        param.feature.Fs=param.Fs;
        param.feature.dT=param.Nfft*(1-param.ovlap)/param.Fs; %Time duration of image pixel
        %save(sprintf('before_cross_channel_%s',date_range(Idate,:)));
        try
            %run_options.debug.Ncalls=50;
            run_options.debug.Ncalls=Inf;
            %Note that we link all DASARs, even 'bad' stations, because we
            %assume that the hydrohpones on the 'bad' stations are OK, and
            %thus can be used to link non-adjacent good DASARS...
            [locations, locations_ctime,linking_index]=cross_channel_match_min(Isite,station,param.feature,goodFile,run_options,run_options.debug);
            
            
        catch
            disp(sprintf('Crash! cross_channel_match: %s',goodFile{1}));
            locations=[];
            locations_ctime=[];
            errormsg=lasterror;
            disp(errormsg.message);
            if (length(errormsg.stack)>0)
                disp(errormsg.stack(1).file)
                disp(errormsg.stack(1).name)
                disp(errormsg.stack(1).line)
            end
            
        end
    end

%     function compute_bearings_for_locations
%         try  %Compute bearings for every station, even if not to be used in final localization
%             
%             tic
%             disp('Computing bearings...');
%             [locations,tet]=compute_bearings(locations,station,goodFile,param,run_options.bearing_alg,0,run_options.kappa_Nsamples);
%             toc
%             %save(sprintf('after_bearings_%s',date_range(Idate,:)));
%             
%         catch
%             disp(sprintf('crash! compute_bearings: %s',goodFile{1}));
%             errormsg=lasterror;
%             disp(errormsg.message);
%             if (~isempty(errormsg.stack))
%                 disp(errormsg.stack(1).file)
%                 disp(errormsg.stack(1).name)
%                 disp(errormsg.stack(1).line)
%             end
%         end
%     end

    function compute_position_for_locations
        %%Compute final location...
        try
            locations=compute_position(locations,goodFile_loc,goodFile,param,Icase,Isite,run_options);
            
            tsv_out=sprintf('%s_%s',keyword.stage,run_options.localization_alg);
            [fname_tsv,Nwrite]=write_tsv(locations,goodName,tsv_out);
            fid=fopen(sprintf('%s_qualitycheck.txt',tsv_out),'a');
            %if Idate==1
            %for JJJ=1:length(goodName)
            %goodDASAR2{JJJ}=auto.goodName{JJJ}(1:5);
            %fprintf(fid,'DASAR file used: %s\n',goodDASAR2{JJJ});
            %end
            %end
            fprintf(fid,'File: %s, number of detections across stations: %s\n',fname_tsv,mat2str(Nwrite));
        catch
            disp(sprintf('crash! compute_position: %s',goodFile{1}));
            errormsg=lasterror;
            disp(errormsg.message);
            
            if (length(errormsg.stack)>0)
                disp(errormsg.stack(1).file)
                disp(errormsg.stack(1).name)
                disp(errormsg.stack(1).line)
            end
        end
    end

end
