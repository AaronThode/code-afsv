%%%master_neural_net_train.m%%
%% Compare output of automatic detector with Greeneridge manual results
%% Aaron Thode
%%  Sept. 18, 2008

function master_neuralnet_train
clear ;close all;fclose('all');

%diary Net_training
%Icase='Site_4_short.morph.crosschecked.Shell07';  %Icase
%params_chc='Shell08_nnet';


[base_name,date_str_local,Site_vector,param_chc,DASAR_str_local]=load_local_runparams(mfilename);

%[base_name,date_str,Site_vec,param_chc]=load_local_runparams(mfilename);

%param=auto_param;
for Itype=1:-1:1  %cycle through different types of biologics
    fprintf('Running network type %i\n',Itype);
    Nwhale=zeros(5,4);
    Nmiss=zeros(5,4);
    Nfalse=zeros(5,4);
    Nred_true=zeros(5,4);
    Nred_false=zeros(5,4);
    for Isite=Site_vector
        Icase=sprintf(base_name,Isite);  %keyword to select time and spatial subset of deployment...
        disp(Icase);
        
        %Note: neural network parameters stored here!
        master_setup;
        
        %%all_biologics
        %If zero, load whales only as true, pinnipeds as false.
        %If one, load all whales,pinnipeds and unknown marine mammal as true,
        %  everything else as false.
        %If two, use 2008 option, load whales only as true, everything else as
        %  false...
        run_options.all_biologics=Itype-1
        run_options.max_hrs=12;
        
        %run_options.strict_eval_tolerance:  If one, manual and auto detections must
        %               overlap 50% in time AND frequency.  If two, compare manual
        %               and auto detections via time tolerances only, with no overlap
        %               requirement.  If three, try to match
        %               by location...
        run_options.strict_eval_tolerance=1;
        run_options.eval_tol=2;  %Time buffer for comparing a manual vs. automated detection time window.  Used only if strict_eval_tolerance=2;
        %           what do do with multiple automated detections that match a
        %           given manual detection?
        run_options.redundant_detections='count_as_true';  %'count_as_true' or 'ignore',
        if run_options.all_biologics==2
            disp('To reproduce 2008 results need to ignore extra detections');
            run_options.redundant_detections='ignore';
        end
        
        run_options.min_stations=[1 9];  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria
        run_options.success_only=0;  %Load only manual  results that have been succesfully localized..
        run_options.center_dist_limit=Inf;  %distance in meters; only used if .success_only==1
        
        %% Neural network options...
        run_options.train_chc='trainscg';  %'trainscg' or 'trainlm for small values'
        run_options.plot_manual_statistics=0;
        miss_fraction=0.1;
        nhidden=10;
        run_options.pca=0; %Principal component analysis on inputs? Not used for 2008 results
        run_options.debug_plot=0;
        run_options.type='mse'; %cross entropy, mse
        run_options.save_memory=1;  %Don't keep history..
        
        param=TOC_params(param_chc);
        [base_name,date_str]=load_local_runparams(mfilename);
        
        
        for Idate=1:length(date_str{1})
            %master_load_automated_data.m
            
            date_str_new{1}{1}=date_str{1}{Idate};
            date_str_new{2}{1}=date_str{2}{Idate};
            
            [date_range,date_range_cell]=expand_date_range(date_str_new);
            
            DASAR_coords=load_DASAR_coords(Icase,Isite);
            DASAR_coords=DASAR_coords(1:7,:);  %Only main array..
            if all(DASAR_coords(6,:)==0)
                DASAR_coords(6,:)=DASAR_coords(5,:);
            end
            
            
            %Load specific DASAR file names, used to access data.
            for Idd=1:size(date_range,1),
                [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates(date_range(Idd,:), ...
                    Isite,DASAR_str,rawdatadir, Icase);
            end
            
            %if run_options.reload_results==0,
            %    load temporary_data
            %else
            head=readgsif_header(goodFile{1}{1});
            mintime=head.ctbc;
            maxtime=head.ctbc+run_options.max_hrs*3600;
            
            [auto.locations,auto.locations_ctime,auto.stations,auto.raw_stations,automated_dir,auto_param,auto.error_area]= ...
                load_automated_results(keyword,detectiondir, goodName,date_range,maxtime,run_options,DASAR_coords);
            
            %%%%%%%%%%%%%%Load manual data%%%%%%%%%%%%%%%
            date_str_new{1}{1}=date_str{1}{Idate};
            date_str_new{2}{1}=date_str{2}{Idate};
            
            [date_range,date_range_cell]=expand_date_range(date_str_new);
            
            
            %Load specific DASAR file names, used to access data.
            
            for Idd=1:size(date_range,1),
                [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates(date_range(Idd,:), ...
                    Isite,DASAR_str,rawdatadir, Icase);
                % DASAR_coords{Idd}=load_DASAR_coords(Icase,Site_vector,goodFile{Idd});
                
            end
            
            head=readgsif_header(goodFile{1}{1});
            mintime=head.ctbc;
            maxtime=head.ctbc+run_options.max_hrs*3600;
            
            %Version that uses distances from array centerline as a
            %filter...DASAR_coords only used to compute centerline, not filter columns
            %of data...
            
            DASAR_coords2=load_DASAR_coords(Icase,Isite);
            %DASAR_coords2=DASAR_coords2(1:7,:);  %Only main array..
            if all(DASAR_coords2(6,:)==0)
                DASAR_coords2(6,:)=DASAR_coords2(5,:);
            end
            
            if run_options.all_biologics==1
                disp('Loading all biological signals into manual');
                [manual.localized,manual.individual,Istrip]=load_manual_results_TSV(manualdir, ...
                    goodDASAR,date_range,maxtime,run_options,DASAR_coords2,Isite);
            elseif run_options.all_biologics==0
                disp('Loading pinniped signals into separate variable');
                [manual.localized,manual.individual,Istrip,pinniped.localized,pinniped.individual]=load_manual_results_TSV(manualdir, ...
                    goodDASAR,date_range,maxtime,run_options,DASAR_coords2,Isite);
                
            elseif run_options.all_biologics==2
                disp('Loading all whale signals into manual, 2008 default');
                [manual.localized,manual.individual,Istrip]=load_manual_results_TSV(manualdir, ...
                    goodDASAR,date_range,maxtime,run_options,DASAR_coords2,Isite);
                
            end

            if Idate==1
                for I=1:length(goodDASAR{1})
                    whale_out{I,Isite}=[];
                    other_out{I,Isite}=[];
                end
            end
            
            
            %Station features that are not allowed to be used in neural
            %network..
            try
                run_options.skip_feature=param.net.skip_feature;
            catch
                run_options.skip_feature={'duration_debug','ctime','ctime_min','ctime_debug','miss'};
            end
            
            %Compute maximum ctime to permit..
            [whale_temp,other_temp,Ntrue,first_time,maxx,stats_all{Idate}]=evaluate_individual_stations(param,manual.individual{1},auto.stations{1}, ...
                auto.raw_stations{1},goodDASAR{1},run_options,mintime,maxtime);
            
            %Display debug info
            display_count_number(whale_temp,'whale_temp');
            display_count_number(other_temp,'other_temp');
            
            
            
            if run_options.all_biologics==0
                fprintf('Pinnipeds will be false class, searching for matches with %i manual pinniped locs\n',length(pinniped.individual{1}.ctime(:,1)));
                
                [other_temp,junk1,junk2,junk3,junk4,stats_seal{Idate}]=evaluate_individual_stations(param,pinniped.individual{1},auto.stations{1}, ...
                    auto.raw_stations{1},goodDASAR{1},run_options,mintime,maxtime);
                display_count_number(other_temp,'other_temp, pinnipeds');
                
                
            end
            %Remove assymetric structure from whale and other
            if ~isempty(Ntrue)
                for K=1:length(whale_temp),
                    whale_temp{K}=rmfield(whale_temp{K},{'Image','param','equalization'});
                    whale_out{K,Isite}=merge_feature_data(whale_out{K,Isite},whale_temp{K});
                    
                    if ~isempty(other_temp)
                        other_temp{K}=rmfield(other_temp{K},{'Image','param','equalization'});
                        other_out{K,Isite}=merge_feature_data(other_out{K,Isite},other_temp{K});
                    end
                    
                end
                
                display_count_number(whale_out(:,Isite),'whale_out to date:');
                display_count_number(other_out(:,Isite),'other_out to date:');
                
                % pause(5);
                
                Nwhale(Isite,:)=Nwhale(Isite,:)+sum(stats_all{Idate}.Nwhale);
                Nred_true(Isite,:)=Nred_true(Isite,:)+sum(stats_all{Idate}.Nred);
                
                Nmiss(Isite,:)=Nmiss(Isite,:)+sum(stats_all{Idate}.Nmiss);
                if run_options.all_biologics==0
                    try
                        Nfalse(Isite,:)=Nfalse(Isite,:)+sum(stats_seal{Idate}.Nwhale);
                        Nred_false(Isite,:)=Nred_false(Isite,:)+sum(stats_seal{Idate}.Nred);
                       
                    end
                else
                    try
                        Nfalse(Isite,:)=Nfalse(Isite,:)+sum(stats_all{Idate}.Nfalse);
                        Nred_false(Isite,:)=Nred_false(Isite,:)+0;
                        
                    end
                    
                end
            end
            
        end %Idate
        display_count_number(whale_out(:,Isite),sprintf('Final whale_out to date, Site %i:',Isite));
        display_count_number(other_out(:,Isite),sprintf('Final other_out to date, Site %i:',Isite));
        
    end %Isite
    Nwhale
    Nmiss
    Nfalse
    Nred_true
    Nred_false
    
    
    clear auto
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Option to load specific analysis of seal calls at certain sites and
    %%%days %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if run_options.all_biologics==0  %Banish pinnipeds to 'false' category
        [other_temp,Isite,stats_extra]=load_additional_days(param.net.extra,detectiondir,rawdatadir,base_name,run_options);
        display_count_number(other_temp,'additional pinnipeds from extra files:');
        
        for K=1:length(other_temp)
            other_out{K,Isite}=merge_feature_data(other_out{K,Isite},other_temp{K});
            
        end
        Nfalse(Isite,:)=Nfalse(Isite,:)+sum(stats_extra.Nwhale);
        Nred_false(Isite,:)=Nred_false(Isite,:)+sum(stats_extra.Nred);
        display_count_number(other_out(:,Isite),sprintf('Total other_out at Site %i now:',Isite));
        
        disp('Nfalse after additional data collected.')
        Nfalse
        Nred_false
        
        Nred_true_total=sum(Nred_true(:,4));
        Nred_false_total=sum(Nred_false(:,4));
        
    end
    Nwhale_total=sum(Nwhale(:,4))
    Nother_total=sum(Nfalse(:,4))
   
    save Training_set_temp
    [pattern,target,historry,net,threshold,preprocess_info,preprocess_info_pca,feature_names]=build_nnetwork_memory(whale_out,Nwhale_total, ...
        other_out,Nother_total,miss_fraction,nhidden,run_options);
    network_save_name=sprintf('Nnetwork_all_Sites%ito%i_%s_%icells_%s_%s',min(Site_vector),max(Site_vector),param.label,size(net.IW{1},1),date_str{1}{1},date_str{2}{end});
    if run_options.all_biologics==1
        network_save_name=[network_save_name '_NonBioVsBio'];
    elseif run_options.all_biologics==0
        network_save_name=[network_save_name '_WhaleVsPinniped'];
    elseif run_options.all_biologics==2
        network_save_name=[network_save_name '_WhaleVsNotWhale'];
    end
    
    keyboard
    mydir=pwd;
    cd(automated_dir);
    cd ..
    !mkdir Neural_Networks.dir
    cd Neural_Networks.dir
    fprintf('Writing results to %s/%s',pwd,network_save_name);
    save(network_save_name,'net','param','threshold','preprocess_info','preprocess_info_pca','nhidden','miss_fraction','feature_names','Nwhale','Nmiss','Nfalse','pattern','target','historry');
    figure(gcf);title(network_save_name)
    print('-djpeg',[network_save_name '.jpg']);
    param.net.dir=pwd;
    
    cd(mydir);
    param.net.name{1}=network_save_name;
    
    Npass_all=0;
    Nall_total=0;
    for Isite=Site_vector
        
        
        for Idate=1:length(date_str{1})
            master_load_automated_data;
            fprintf('Applying Site %i, %s neural net to original data \n',Isite, date_range(1,:));
            [station,Npass,Nall]=filter_with_Nnet(auto.stations{1},param.net.name,param.net.dir);
            Npass_all=Npass_all+sum(sum(Npass));
            Nall_total=Nall_total+sum(sum(Nall));
            
        end
        
        
    end
    
    
end %Itype
end


function [other_out,Isite,stats] = load_additional_days(input_info,detectiondir,rawdatadir,base_name,run_options,mintime,maxtime)
%%Load odd TSV files that contain addition pinniped or non-whale biological
%%information.
%date_str,Idate,Site_vector,Isite,DASAR_str,rawdatadir,manualdir,Icase
% run_options.max_hrs,
mydir=pwd;

%%Ensure that pinnipeds will be assigned to 'false' category
run_options.all_biologics=0;
        
%input_info{I}.dir='/Users/Shared/Projects/2009_Arctic_Analysis/TSV_files_Shell09/Shell09_PinnipedAnalysis';
for I=1:length(input_info);
    cd(input_info{I}.dir)
    fname=dir('*.tsv');
    cd(mydir);
    for J=1:length(fname)
        date_str{1}{1}=['20' fname(J).name(9:10) fname(J).name(3:6)];
        date_str{2}{1}=date_str{1}{1};
        Idate=1;
        Site_vector=str2num(fname(J).name(8));
        Isite=Site_vector;
        Icase=sprintf(base_name,Isite);  %keyword to select time and spatial subset of deployment...
        Idots=findstr(Icase,'.');
        if isempty(Idots),
            keyword.scenario=Icase;
        elseif length(Idots)==1,
            keyword.scenario=Icase(1:(Idots(1)-1));
            keyword.algorithm=Icase((Idots(1)+1):end);
        elseif length(Idots)==2,
            keyword.scenario=Icase(1:(Idots(1)-1));
            keyword.algorithm=Icase((Idots(1)+1):(Idots(2)-1));
            keyword.stage=Icase((Idots(2)+1):end);
        elseif length(Idots)==3,
            keyword.scenario=Icase(1:(Idots(1)-1));
            keyword.algorithm=Icase((Idots(1)+1):(Idots(2)-1));
            keyword.stage=Icase((Idots(2)+1):(Idots(3)-1));
        end
        
        DASAR_str='*';
        manualdir=input_info{I}.dir;
        
        
        date_str_new{1}{1}=date_str{1}{Idate};
        date_str_new{2}{1}=date_str{2}{Idate};
        [date_range,date_range_cell]=expand_date_range(date_str_new);
        
        %Load specific DASAR file names, used to access data.
        
        for Idd=1:size(date_range,1),
            [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates(date_range(Idd,:), ...
                Site_vector,DASAR_str,rawdatadir, Icase);
        end
        
        head=readgsif_header(goodFile{1}{1});
        mintime=head.ctbc;
        maxtime=head.ctbc+run_options.max_hrs*3600;
        
        DASAR_coords=load_DASAR_coords(Icase,Site_vector);
        if J==1
            for I=1:length(goodDASAR{1})
                other_out{I}=[];
            end
        end
        
        %Since all_biologics==0 pinnipeds are assigned to 'false' category
        [manual.localized,manual.individual,Istrip,pinniped.localized,pinniped.individual]=load_manual_results_TSV(manualdir, ...
            goodDASAR,date_range,maxtime,run_options,DASAR_coords,Isite);
        
        [auto.locations,auto.locations_ctime,auto.stations,auto.raw_stations,automated_dir,auto_param,auto.error_area]= ...
            load_automated_results(keyword,detectiondir, goodName,date_range,maxtime,run_options,DASAR_coords);
        %[whale_temp,other_temp,Ntrue,first_time,max_ctime,stats_all{Idate}
        [other_temp,junk1,junk2,junk3,junk4,stats]=evaluate_individual_stations(auto_param,pinniped.individual{1},auto.stations{1}, ...
            auto.raw_stations{1},goodDASAR{1},run_options,mintime,maxtime);
        
        for K=1:length(other_temp)
            other_temp{K}=rmfield(other_temp{K},{'Image','param','equalization'});
            other_out{K}=merge_feature_data(other_out{K},other_temp{K});
        end
        
    end
end
cd(mydir);
end


function display_count_number(data,strr)
fprintf('Length of %s ',strr);
Ntotal=0;
for I=1:length(data)
    if isempty(data{I})
        continue
    end
    count=length(data{I}.ctime_min);
    fprintf(' station %i: %i ',I,count);
    Ntotal=Ntotal+count;
end
fprintf('\nTotal %s count: %i\n\n',strr,Ntotal);
end

%%Nall_total and Nwhale+Nfalse should be close (close duplicates in time removed)

% for Idate=1:length(date_str{1})
%     %Idate=menu('Which date?',date_str{1});
%     date_str_new{1}{1}=date_str{1}{Idate};
%     date_str_new{2}{1}=date_str{2}{Idate};
%
%     [date_range,date_range_cell]=expand_date_range(date_str_new);
%
%     [auto.locations,auto.locations_ctime,auto.stations,auto.raw_stations,automated_dir,auto_param]= ...
%         load_automated_results(keyword,detectiondir, goodName,date_range,maxtime,run_options);
%
%     for Istation=1:length(whale),
%         disp(sprintf('Applying %s neural net to original data',date_range(1,:)));
%         station_tmp=auto.stations{1}(Istation);
%         station(Istation)=filter_with_Nnet(station_tmp,param,automated_dir);
%     end
% end
%
