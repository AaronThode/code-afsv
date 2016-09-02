%%%%%%%master_setup.m%%%%%%

%clear ;
fclose('all');
mydir=pwd;
[rawdatadir,Icase,detectiondir,param,manualdir]=load_pathnames(Icase,[]);

if exist('rawdatadir_local')
   rawdatadir=rawdatadir_local; 
end

if exist('detectiondir_local')
   detectiondir=detectiondir_local; 
end

%%%%%%%%%CHANGE THESE%%%%%%%%%%%%%%%%%%%%%
 %run_options.strict_eval_tolerance:  If one, manual and auto detections must
%               overlap 50% in time AND frequency.  If two, compare manual
%               and auto detections via time tolerances only, with no overlap
%               requirement.  If three, try to match
%               by location...
run_options.strict_eval_tolerance=2;
run_options.eval_tol=2;  %Time tolerance for comparing a manual vs. automated detection

%           what do do with multiple automated detections that match a
%           given manual detection?
run_options.redundant_detections='count_as_true';  %'count_as_true' or 'ignore', 

 %If zero, load whales only as true, pinnipeds as false. 
 %If one, load all whales,pinnipeds and unknown marine mammal as true,
 %  everything else as false.
 %If two, use 2008 option, load whales only as true, everything else as
 %  false...
run_options.all_biologics=2; 


%Option for limiting times to be considered
run_options.max_hrs=24;

%Option to permit feature filtering of automated locations
run_options.auto_location_filter=0;

%Option for restricting number of manual stations involved in localization

run_options.min_stations=[2 9];  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria
run_options.success_only=0;  %Load only manual  results that have been succesfully localized..
run_options.center_dist_limit=100000;  %distance in meters; only used is .success_only==1;

run_options.plot_manual_statistics=0;  %Histogram features from manual TSV files
run_options.limit_manual_ctimes=0;  %If one, limit manual ctimes load to maximum value of automated results
run_options.debug_miss_compare=0;  %If one, plot debug output of evaluate_miss_fraction

run_options.auto_success_only=1; %Load only auto  results that have been succesfully localized..
run_options.auto_min_stations=[2 9];


% morph options
param.morph.threshold_chc='dual_threshold';

%compute_position options..
run_options.localization_alg={'keep stored result','Huber','medianHuber','HuberCalibratedKappa','Andrews'};
%Iloc_chc=menu('Select localization method:',run_options.localization_alg);
Iloc_chc=1;
run_options.min_weight=0;  %Minimum Huber weight required to include call in match...
run_options.Ikeep_only=0;


%% Neural network options...
% run_options.train_chc='trainscg';  %'trainscg' or 'trainlm for small values'
% run_options.plot_manual_statistics=0;
% miss_fraction=0.1;
% nhidden=10;
% run_options.pca=1; %Principal component analysis on inputs?
% run_options.debug_plot=0;
% run_options.type='mse'; %cross entropy, mse
% run_options.save_memory=1;  %Don't keep history..


%%%%%%%%%%%%%DON'T TOUCH BELOW!!!! %%%%%%%%%%%%%%%%%%%%%%%%%
%Download parameters of the automated run...
%   Site_vector: Site being examined
%   date_str: [Nday x 1] cell matrix that contains string of appropriate
%       day
%   DASAR_str: string of single-letters representing desired stations
%   keyword: contains 'scenario', 'algorithm', 'stage' fields, used to load
%       precise automated results desired

 run_options.localization_alg=run_options.localization_alg{Iloc_chc};
   
if Iloc_chc>1
    disp(sprintf('Recomputing locations using %s algorithm', run_options.localization_alg));
    run_options.plot_locations=0;
end

[~,date_str,DASAR_str,keyword,detection_stage]=TOC_runs(Icase);
%Load default dates and DASAR limits (bad DASAR list)
% [~,date_str,DASAR_str,keyword]=TOC_runs(Icase);
if ~isempty(date_str_local)
    date_str=date_str_local;
end
if ~isempty(DASAR_str_local)
    DASAR_str=DASAR_str_local;
end


run_options.calculate_airgun_characteristics=0;  %Never change..

run_options.filtering_stage=keyword.algorithm;  %'contour','morph','both','none': Determines what third stage processing takes place..
