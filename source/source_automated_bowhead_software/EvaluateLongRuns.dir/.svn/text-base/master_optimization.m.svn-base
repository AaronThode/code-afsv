%%%%%%%%%master_optimizaton.m%%%%%%%%
% Program for optimizing detection parameters for energy detection and
% contour tracers
%% 

%%%DELPHINE, define large variables here so they can be accessed by 'optimization_function' 
global manual     run_options param opt_param_chc  final Iruns

!rm master_optimization_run.txt
diary master_optimization_run.txt
Icase='Shell08_Site1_allHuber.morph.crosschecked';  %keyword to select time and spatial subset of deployment...

master_setup;

run_options.debug.raw=0;    %If one, show debug energy detector information
run_options.debug.interval=0;
run_options.debug.morph=0;  %If one, show debug output for morph processing, and load debug_ctimes to learn when "failure" occures
run_options.debug.snips=0;  %If one, compare snips data with direct download of raw data...

run_options.Ncalls_to_sample=200;  %Number of snip samples to read into memory at once
run_options.force_energy_detector_run=1;  %If 1, always force JAVA cell CFAR to run
run_options.load_single_DASAR_results=0;  %If one, load stage 3 output (morph processing) from file

params_chc='Shell08_initial';
param=TOC_params(params_chc);  %Unpack detection parameters from keyword

%param.energy.nstart=4*60*60*1000;
run_options.max_hrs=12;
param.energy.nsamples=run_options.max_hrs*60*60*1000;


param.energy.Fs=param.Fs;
param.morph.threshold_chc='reconstruction';

alg_chc='direct'; %'direct' or 'GA'
run_options.opt_param_chc='Morph';

final.criteria.tol=3;
final.criteria.fit_weights=[0.5*10^6 1;1000 1;100 1; 10 1];
final.criteria.fit_goal=[0.05 0.25];
final.Idebug.val=0;
final.criteria.filtering_stage='morph';  %'contour','morph','both','none'


Idate=menu('Which date?',date_str{1});
 param_org=TOC_params(params_chc);
        
master_load_manual_data;

run_options.max_ctime=maxtime;

if length(goodFile)>1,
    error('Only one station can be processed at a time');
end
    
objfun= @(opt_vec) optimization_function(opt_vec,manual,date_str_new,goodFile,goodName,goodDASAR,run_options,param);

%%Begin criteria loop,
for Icrit=1:size(final.criteria.fit_weights,1),
    Iruns=0;
    disp(sprintf('Criteria %i',Icrit));
    final.criteria.fit_weight=final.criteria.fit_weights(Icrit,:);
    disp(sprintf('fit weight is %s',mat2str(final.criteria.fit_weight)));
    %Now select which parameters to optimize..
    %param=param_org;
    
    %%DELPHINE, this routine is useful to convert a set of variable you
    %%want to optimize.  The desired variables are fields within the
    %%structure 'param'.  When you feed the 'extract' string into the
    %% subroutine, the program extracts the values from each variable and
    %% creates a vector 'opt_vec' that can then be fed into the
    %% optimization program.  You'll see that inside the
    %% 'optimization_function.m'  I'll use the same program to 
    %% unpack the vector into distinct variables using an 'assign' string.
    
    [param, opt_vec,LB,UB]=extract_optimization_parameters(param_org,[],'extract',run_options.opt_param_chc);
    % keyboard;
    %%Test optimization function
    %opt_vec= [0.195 0.429 0.304 0.0284 0.0147 0.45 0.521 0.979 0.006 0.791 0.0686 0.101 0.55 0.375 0.22 0.91 0.493 0.0603];
    %opt_vec=[0.382 0.241 0.304 0.153 0.0147 0 0.443 0.167 0.333 0.641 0.006 0.478 0.319 0.101 0.691 1 0.0479 0.875 0.681 0.0603 0.981 0 1]

    dumm=objfun(opt_vec);
   
    if strcmp(alg_chc,'GA')
        %options=gaoptimset('PopulationSize',40,'popinitrange',[LB; UB]);

        %DELPHINE, this line is used to set parameters for running the
        %genetic algorithm.  The values are covered in detail in the
        %instructions for the Direct Search/Genetic Algorithm toolbox
        options=gaoptimset('PopulationSize',[20 20],'Display','diagnose','InitialPopulation',opt_vec,'StallTimeLimit',Inf);

        %%DELPHINE, command for running the genetic algorithm.  Note that
        %%it calls optimization_function.m and uses the output of
        %%extract_parameters.m 
        [fin_vec fvall reason output population scores]=ga(@optimization_function,length(opt_vec), ...
            [],[],[],[],LB,UB,[],options);
        
        %param_org=TOC_params(params_chc);
        %DELPHINE, I use this to translate the optimized output into
        %distinct variables.  Note I use 'assign' instead of 'extract' here
        
        [paramGA, opt_vec,LB,UB]=extract_optimization_parameters(param_org, fin_vec,'assign',opt_param_chc);
        Iruns
        save finGAresult.mat paramGA fin_vec LB UB fvall param_org final opt_param_chc load_false_reviews extract_high_SNR_only Iruns fvall reason output population scores
        disp('Raw GA result finished');
        %options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.02,'TolMesh',.01, ...
        %    'TolX',0.005, 'Display','Diagnose','InitialMeshSize',0.3,'PollMethod','MADSPositiveBasis2N');
        
        %%DELPHINE, here are the commands to run a direct search instead of
        %%GA.  I would recommend trying this optimization scheme first, it
        %%works pretty well using these parameters...
        
        options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.1,'TolMesh',.1, ...
            'TolFun',0.1, 'Display','Diagnose','InitialMeshSize',0.3,'PollMethod','MADSPositiveBasis2N');
        [out_vec fval]=patternsearch(@optimization_function, fin_vec, [],[],[],[],LB,UB,[],options);

        param_org=TOC_params(params_chc);
        [param, opt_vec,LB,UB]=extract_optimization_parameters(param_org, out_vec,'assign',opt_param_chc);

        %%all_processing.m%%%
        param.no_energy_run_needed=0;

        dumm=optimization_function(out_vec);
        save finGADSresult.mat param paramGA fin_vec out_vec LB UB fval param_org final opt_param_chc load_false_reviews extract_high_SNR_only fvall reason output population scores
        eval(sprintf('!mv finGADSresult.mat finGAresult_criteria%i.mat',Icrit));

    else
        disp('Starting direct search');
        %Run 1 options:
        %options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.025,'TolMesh',.001, ...
        %'Display','Iter','InitialMeshSize',0.25);
        %Run 2 options
       % options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.0015,'TolMesh',.001, ...
        %    'TolFun',0.01, 'Display','Diagnose','InitialMeshSize',0.3,'PollMethod','MADSPositiveBasis2N');

%	options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.0015,'TolMesh',.001, ...
%            'TolFun',0.05, 'Display','Diagnose','InitialMeshSize',0.3,'PollMethod','MADSPositiveBasis2N');

        options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.1,'TolMesh',.1, ...
            'TolFun',0.1, 'Display','Diagnose','InitialMeshSize',0.3,'PollMethod','MADSPositiveBasis2N');
        
        %options = psoptimset('CompletePoll','On','Cache','On','CacheTol',.3,'TolMesh',.3, ...
         %   'MaxIter',3,'TolX',0.05, 'Display','Diagnose','InitialMeshSize',0.3,'PollMethod','MADSPositiveBasis2N');


        %[out_vec fval]=patternsearch(objfun, opt_vec, [],[],[],[],LB,UB,[],options);
        [out_vec fval]=patternsearch(objfun, opt_vec, [],[],[],[],[],[],[],options);
        
        %[out_vec fval]=patternsearch(@optimization_function, out_vec1, [],[],[],[],LB,UB,[],options);

        [param, opt_vec]=extract_optimization_parameters(param_org, out_vec,'assign',run_options.opt_param_chc);
        

        %%all_processing.m%%%
        param.no_energy_run_needed=0;

        %dumm=optimization_function(out_vec);
        %keyboard;
        save finDSresult param out_vec fval param_org final opt_param_chc 
        eval(sprintf('!mv finDSresult.mat finDSresult_criteria%i.mat',Icrit));
       

    end
end
