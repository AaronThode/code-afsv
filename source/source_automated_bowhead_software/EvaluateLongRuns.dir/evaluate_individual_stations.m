%%%%evaluate_individual_stations.m%%%%
% function [whale,other,Ntrue]=evaluate_individual_stations(param,individual,stations,raw_stations,goodDASAR,run_options)
%    Divides station detections into "whale" and "other" components, and
%    Ntrue is number of true individual detections across stations.

function [whale,other,Ntrue,first_time,max_ctime,stats_overall]=evaluate_individual_stations(param,manual_indiv,stations,raw_stations,goodDASAR,run_options,min_ctime,max_ctime)

whale=[];
other=[];
Ntrue=[];
first_time=[];
%max_ctime=[];
stats_overall=[];

if size(manual_indiv.ctime,1)==0
    return
end
%Set maximum time to evaluate--are we evaluating only a partial day?
bufferTol=run_options.eval_tol;
first_time=manual_indiv.ctime(1,manual_indiv.ctime(1,:)>0);


feature_names=stations(1).param.feature_name;

%Count total number of call detections across all DASARS, whether used for
%localization or not.
Ntrue=0;
for K=1:length(goodDASAR),
    Ntrue=Ntrue+length(find(manual_indiv.ctime(:,K)>0));
end

Nwhale=zeros(length(goodDASAR),4);
Nmiss=zeros(length(goodDASAR),4);
Nfalse=zeros(length(goodDASAR),4);
Nred=zeros(length(goodDASAR),4);
for K=1:length(goodDASAR)
    disp(sprintf('vvvvvvvvvvvvv Begin Station %s vvvvvvvvv',goodDASAR{K}));
    Nmanual_calls=valid_call_count(manual_indiv.ctime(:,K),min_ctime,max_ctime);
    
    
    Nauto_calls=valid_call_count(raw_stations(K).raw_detections.ctime,min_ctime,max_ctime);
    fprintf('Comparing %i manual calls with %i raw detections, max time  %s\n',Nmanual_calls,Nauto_calls,ctime2str(max_ctime));
    %[Imiss,Ifalse,Itrue_man,Itrue_auto]=find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, max_ctime, auto_features,correct_features)
    
    [Imiss_raw,junk,Itrue_man_raw]=find_similar_elements_ovlap_timelimit(raw_stations(K).raw_detections.ctime, ...
        raw_stations(K).raw_detections.duration, ...
        manual_indiv.ctime(:,K), manual_indiv.duration(:,K),min_ctime,max_ctime);
    
    %Imiss,Ifalse,Itrue_man,Itrue_auto,
    Nwhale(K,1)=Nwhale(K,1)+length(Itrue_man_raw);
    Nmiss(K,1)=Nmiss(K,1)+length(Imiss_raw);
    Nfalse(K,1)=Nfalse(K,1)+length(junk);
    
    Nauto_calls=valid_call_count(raw_stations(K).interval_detections.ctime,min_ctime,max_ctime);
    fprintf('Comparing %i manual calls with %i interval detections, max time  %s\n',Nmanual_calls,Nauto_calls,ctime2str(max_ctime));
    
    [Imiss_interval,junk,Itrue_man_interval]=find_similar_elements_ovlap_timelimit(raw_stations(K).interval_detections.ctime, ...
        raw_stations(K).interval_detections.duration, manual_indiv.ctime(:,K), ...
        manual_indiv.duration(:,K), min_ctime,max_ctime);
    
    Nwhale(K,2)=Nwhale(K,2)+length(Itrue_man_interval);
    Nmiss(K,2)=Nmiss(K,2)+length(Imiss_interval);
    Nfalse(K,2)=Nfalse(K,2)+length(junk);
    
    try
        Nauto_calls=valid_call_count(raw_stations(K).morph_detections.ctime,min_ctime,max_ctime);
        fprintf('Comparing %i manual calls with %i morph detections, max time  %s\n',Nmanual_calls,Nauto_calls,ctime2str(max_ctime));
        
        [Imiss_morph,junk,Itrue_man_morph]=find_similar_elements_ovlap_timelimit(raw_stations(K).morph_detections.ctime, ...
            raw_stations(K).morph_detections.duration, manual_indiv.ctime(:,K), ...
            manual_indiv.duration(:,K), min_ctime,max_ctime);
        
        Nwhale(K,3)=Nwhale(K,3)+length(Itrue_man_morph);
        Nmiss(K,3)=Nmiss(K,3)+length(Imiss_morph);
        Nfalse(K,3)=Nfalse(K,3)+length(junk);
        
    end
    Nauto_calls=valid_call_count(stations(K).ctime_min,min_ctime,max_ctime);
    
    fprintf('Comparing %i manual calls with %i station detections, max time  %s\n',Nmanual_calls,Nauto_calls,ctime2str(max_ctime));
    switch run_options.strict_eval_tolerance
        case 1  %strict, match frequecy, time, and duration
            %[Imiss,Ifalse,Itrue_man,Itrue_auto]=find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, max_ctime, auto_features,correct_features)
            
            [Imiss,Ifalse,Iman_pass,Itrue,Iredundant]=find_similar_elements_ovlap_timelimit(stations(K).ctime_min, ...
                stations(K).Totalduration, manual_indiv.ctime(:,K), ...
                manual_indiv.duration(:,K),min_ctime,max_ctime, ...
                [stations(K).Totalfmin; stations(K).Totalfmax],[manual_indiv.flo(:,K) manual_indiv.fhi(:,K)]');
            
        case 2  %match time only
            
            [Imiss,Ifalse,Iman_pass,Itrue,Iredundant]=find_similar_elements_ovlap_timelimit(stations(K).ctime_min, ...
                stations(K).Totalduration, manual_indiv.ctime(:,K)-bufferTol, ...
                manual_indiv.duration(:,K)+2*bufferTol,min_ctime,max_ctime);
            
            
    end
    
    %%Decide whether to include redundancies as "correct" automated
    %%detections, or just ignore
    if strcmp(run_options.redundant_detections,'count_as_true')
        fprintf('Merging %i redundant auto detections into %i auto detections, ',length(Iredundant),length(Itrue));
        Itrue=unique(sort([Itrue Iredundant]));
        fprintf('%i elements in final true vector \n',length(Itrue))
    elseif strcmp(run_options.redundant_detections,'ignore')
        
    end
    
    %%Update statistics
    Nwhale(K,4)=Nwhale(K,4)+length(Itrue);
    Nmiss(K,4)=Nmiss(K,4)+length(Imiss);
    Nfalse(K,4)=Nfalse(K,4)+length(Ifalse);
    Nred(K,4)=Nred(K,4)+length(Iredundant);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Gather  statistics and features on missed calls...   %%%%
    %%%  not used for neural network processing..          %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %collect_miss_statistics;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Construct output cell matricies 'whale' and 'other'  %%%%
    %%%  used to train neural network                      %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    feature_names=stations(K).param.feature_name;
    %Itrue=setdiff(1:length(stations(K).ctime),Ifalse);
    
    %global_features={'Nsegments','Totalfmin','Totalfmax','Totalduration','Contour_Area','Contour_global_bandwidth','Contour_duration', ...
    % 'Contour_local_bandwidth','Contour_fmin','Contour_fmax','ctime_min','ctime_debug','duration_debug','SNR','SEL'};
    global_features=param.feature.global_names;
    global_features=([global_features {'ctime_min','ctime_debug','duration_debug'}]);
    for I=1:length(global_features)
        whale{K}.(global_features{I})=stations(K).(global_features{I})(1,Itrue);
        other{K}.(global_features{I})=stations(K).(global_features{I})(1,Ifalse);
    end
    
    %%We're changing some global variable names
    whale{K}.SNR_global=stations(K).SNR(1,Itrue);
    other{K}.SNR_global=stations(K).SNR(1,Ifalse);
    
    whale{K}.SEL_global=stations(K).SEL(1,Itrue);
    other{K}.SEL_global=stations(K).SEL(1,Ifalse);
    
    whale{K}.Image=stations(K).Image(Itrue);
    other{K}.Image=stations(K).Image(Ifalse);
    
    whale{K}.param.equalization_freq=stations(K).param.equalization_freq;
    whale{K}.equalization=stations(K).equalization(:,Itrue);
    
    
    other{K}.param.equalization_freq=stations(K).param.equalization_freq;
    other{K}.equalization=stations(K).equalization(:,Ifalse);
    
    %Finally, features of individual segments, first one only (greatest
    %area)
    for I=1:length(feature_names)
        whale{K}.(feature_names{I})=stations(K).feature.(feature_names{I})(1,Itrue);
        other{K}.(feature_names{I})=stations(K).feature.(feature_names{I})(1,Ifalse);
        
    end
    
    disp(sprintf('^^^^^^^^^ End Station %s ^^^^^^^',goodDASAR{K}));
    
end


stats_overall.Nwhale=Nwhale;
stats_overall.Nmiss=Nmiss;
stats_overall.Nfalse=Nfalse;
stats_overall.Nred=Nred;

    function collect_miss_statistics
        
        Imiss_morph_step=setdiff(Imiss_morph,Imiss_raw);
        Imiss_interval_step=setdiff(Imiss_interval,Imiss_raw);
        Nmanual=size(manual_indiv.ctime,1);
        %Iman_pass=setdiff(1:Nmanual,Imiss);
        
        for I=1:length(feature_names),
            switch(feature_names{I}),
                case {'maxSNR','anchor_station','DASARstr'}
                    continue;
                otherwise
                    whale{K}.miss.energy.(feature_names{I})=manual_indiv.(feature_names{I})(Imiss_raw,K);
                    whale{K}.miss.interval.(feature_names{I})=manual_indiv.(feature_names{I})(Imiss_interval,K);
                    whale{K}.miss.morph.(feature_names{I})=manual_indiv.(feature_names{I})(Imiss,K);
                    whale{K}.miss.morph_step.(feature_names{I})=manual_indiv.(feature_names{I})(Imiss_morph_step,K);
                    whale{K}.miss.Itrue.(feature_names{I})=manual_indiv.(feature_names{I})(Iman_pass,K);
                    
            end
            
        end
        whale{K}.miss.energy.diff=miss_diff;
        whale{K}.miss.energy.index=Imiss_raw;
        whale{K}.miss.interval.index=Imiss_interval;
        whale{K}.miss.morph.index=Imiss_morph;
        whale{K}.miss.morph_step.index=Imiss_morph_step;  %detections that do not survive morph processing
        whale{K}.miss.Itrue.index=Iman_pass;
        
        for II=1:length(Imiss_morph_step)
            [junk,Iwant]=min(abs(manual_indiv.ctime(Imiss_morph_step(II),K)-raw_stations(K).raw_detections.ctime));
            whale{K}.miss.morph_step.snips_index(II)=Iwant;
            
        end
        
        %%%Additional data for debug check--useful for confirming that morph
        %%%results can be reproduced...
        for II=1:length(Iman_pass)
            [junk,Iwant]=min(abs(manual_indiv.ctime(Iman_pass(II),K)-raw_stations(K).raw_detections.ctime));
            whale{K}.miss.Itrue.snips_index(II)=Iwant;
            
            [junk,Iwant]=min(abs(manual_indiv.ctime(Iman_pass(II),K)-stations(K).ctime_min));
            whale{K}.miss.Itrue.final_index(II)=Iwant;
        end
        
        %%%%%
        
    end

%%%%%%%plot_3D_data.m
    function plot_3D_data
        for K=1:length(goodDASAR),
            
            subplot(2,1,1);
            plot3(whale{K}.(feature_names{Ichcx}),whale{K}.(feature_names{Ichcy}),whale{K}.(feature_names{Ichcz}),'r.');hold on
            set(gca,'fontweight','bold','fontsize',14);
            xlabel(feature_names{Ichcx});ylabel(feature_names{Ichcy});zlabel(feature_names{Ichcz});grid on;
            %adjust_feature_plot_limits(feature_names,Ichcx,Ichcy);
            
            subplot(2,1,2);
            plot3(other{K}.(feature_names{Ichcx}),other{K}.(feature_names{Ichcy}),other{K}.(feature_names{Ichcz}),'k.');hold on
            set(gca,'fontweight','bold','fontsize',14);
            xlabel(feature_names{Ichcx});ylabel(feature_names{Ichcy});zlabel(feature_names{Ichcz});grid on;
            %[x{1},x{2}]=adjust_feature_plot_limits(feature_names,Ichcx,Ichcy);
        end
    end


%%%histogram_2D_features
    function [Ntot_whale,Ntot_other,F]=histogram_2D_features
        
        for K=1:length(goodDASAR),
            
            for Iname=1:2,
                [N,F{Iname}]=hist(whale{K}.(feature_names{chc(Iname)}),linspace(x{Iname}(1),x{Iname}(2),250));
                if K==1,
                    Ntot_whale{Iname}=N;
                else
                    Ntot_whale{Iname}=Ntot_whale{Iname}+N;
                end
                
                [N,F{Iname}]=hist(other{K}.(feature_names{chc(Iname)}),linspace(x{Iname}(1),x{Iname}(2),250));
                if K==1,
                    Ntot_other{Iname}=N;
                else
                    Ntot_other{Iname}=Ntot_other{Iname}+N;
                end
            end
        end
        
        subplot(2,1,1);
        title(sprintf('%i calls',Nwhale));
        subplot(2,1,2);
        title(sprintf('%i false hits',Nother));
        
    end



%%%%%%%plot_hist.m%%%%%%%
    function plot_hist(N,Iplot)
        
        cum=cumsum(N)/sum(N);
        [junk,F95]=min(abs(0.95-cum));
        [junk,F75]=min(abs(0.75-cum));
        [junk,F10]=min(abs(0.10-cum));
        
        subplot(2,2,Iplot);
        bar(F{Iname},N,'k');ylabel('Number events');xlabel(feature_names{Iname});grid on
        title(sprintf('mean value of %s, %62.f, %6.2f to %6.2f covers 75-95%% range', ...
            feature_names{Iname},F{Iname}(F10),F{Iname}(F75),F{Iname}(F95)));
        
    end
end
function Ncalls=valid_call_count(data,min_ctime,max_ctime)
Ncalls=length(find(data>min_ctime&data<max_ctime));

end

function [xlimm,ylimm]=adjust_feature_plot_limits(names,Ix,Iy)
ind=[Ix Iy];
for I=1:2,
    switch lower(names{ind(I)})
        case 'extent'
            limm=[0 1];
        case 'orientation'
            limm=[-90 90];
        case 'eccentricity'
            limm=[0 1];
        case 'robust_fmin'
            limm=[0 300];
        case 'robust_fmax'
            limm=[50 500];
        case 'time_band_product'
            limm=[0 1e3];
        case 'local_bandwidth'
            limm=[0 300];
        case 'robust_bandwidth'
            limm=[0 300];
        case 'solidity'
            limm=[0 1];
        case 'boundingbox'
            limm=[0 5];
        otherwise
            keyboard
            
    end
    
    if I==1,
        xlim(limm);
        xlimm=limm;
    else
        ylim(limm);
        ylimm=limm;
    end
end


end
