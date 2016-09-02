%%%master_evaluation_neuralnet.m.m%%
%% compare output of automatic detector with greeneridge manual results
%% create roc curve of neural network of choice...
 % stats{Isite}.nmatch(Idate,J,Ithresh1,Ithresh2)--Ithresh1: first neural network threshold. Ithresh2: second
 %  network treshold.  Note that there are two nested loops to permit thresholds to be vectors.
%% updated feb 7, 2011
clear all;close all

%%load case-specific data
[icasestr,date_str_local,site_vector,params_chc,DASAR_str_local]=load_local_runparams(mfilename);

mydir=pwd;
%%%set plot_links to 1 in order to make maps and debug information..
%run_options.plot_links=0;


for Isite=site_vector  % for each site..
    disp(sprintf('Site %i',site_vector(Isite)));
    cd(mydir);
    Icase=sprintf(icasestr,Isite);
    master_setup;
    
    run_options.all_biologics=2;
    %run_options.strict_eval_tolerance:  if one, manual and auto detections must
    %               overlap 50% in time and frequency.  if two, compare manual
    %               and auto detections via time tolerances only, with no overlap
    %               requirement.  if three, try to match
    %               by location...
    run_options.strict_eval_tolerance=1;
    run_options.eval_tol=2;  %time uncertainty estimate for manual selection.
    run_options.match_tol=0.5;  %fractional overlap in time (and frequency) required to declare a match.
    %thresholds1=linspace(-0.9,0.9,19);
    thresholds1=-1:0.2:1.0;
    thresholds2=thresholds1;
    %thresholds2=-1;
    
    for Idate=1:length(date_str{1})
        cd(mydir);
    
        if Idate==1
            stats{Isite}.nmatch=-1*ones(length(Idate),9,length(thresholds1),length(thresholds2));
            stats{Isite}.nmiss=stats{Isite}.nmatch;
            stats{Isite}.nfalse=stats{Isite}.nmatch;
            
        end
        clear station_in station
       
        date_str_new{1}{1}=date_str{1}{Idate};
        date_str_new{2}{1}=date_str{2}{Idate};
        
        [date_range,date_range_cell]=expand_date_range(date_str_new);
        
        dasar_coords=load_dasar_coords(Icase,Isite);
        dasar_coords=dasar_coords(1:7,:);  %only main array..
        if all(dasar_coords(6,:)==0)
            dasar_coords(6,:)=dasar_coords(5,:);
        end
        
        
        %load specific dasar file names, used to access data.
        
        try %%attempt to locate processed data locations first
            
            %%try to locate processed results
            for Idd=1:size(date_range,1)
                [gooddasar{Idd},goodfile{Idd},goodname{Idd},dasar_str_date{Idd}]=find_dasar_dates_processed(date_range(Idd,:), ...
                    Isite,DASAR_str,detectiondir, Icase);
            end
            
            t=goodname{1}{1}(8:22);
            mintime_tabs=datenum(str2num(t(1:4)),str2num(t(5:6)),str2num(t(7:8)),str2num(t(10:11)),str2num(t(12:13)),str2num(t(14:15)));
            mintime=24*3600*(mintime_tabs-datenum(1970,1,1,0,0,0));
            maxtime=mintime+run_options.max_hrs*3600;
            
            
            
        catch  %look at raw data locations
            for Idd=1:size(date_range,1)
                [gooddasar{Idd},goodfile{Idd},goodname{Idd},dasar_str_date{Idd}]=find_dasar_dates(date_range(Idd,:), ...
                    Isite,DASAR_str,rawdatadir, Icase);
            end
            
            %if run_options.reload_results==0,
            %    load temporary_data
            %else
            head=readgsif_header(goodfile{1}{1});
            mintime=head.ctbc;
            maxtime=head.ctbc+run_options.max_hrs*3600;
        end
        
        %%%load manual data...
        
        if run_options.all_biologics==1
            disp('loading all biological signals into manual');
            [manual.localized,manual.individual,istrip]=load_manual_results_tsv(manualdir, ...
                gooddasar,date_range,maxtime,run_options,dasar_coords,Isite);
        elseif run_options.all_biologics==0
            disp('loading pinniped signals into separate variable');
            [manual.localized,manual.individual,istrip,pinniped.localized,pinniped.individual]=load_manual_results_tsv(manualdir, ...
                gooddasar,date_range,maxtime,run_options,dasar_coords,Isite);
            
        elseif run_options.all_biologics==2
            disp('loading all whale signals into manual, 2008 default');
            [manual.localized,manual.individual,istrip]=load_manual_results_tsv(manualdir, ...
                gooddasar,date_range,maxtime,run_options,dasar_coords,Isite);
            
        end
        
        
        %%%load raw automated results from image processing stage
        % [auto.locations,auto.locations_ctime,auto.stations,auto.raw_stations,automated_dir,auto_param,auto.goodname,auto.error_area]= ...
        %       load_automated_results(keyword,detectiondir, goodname,date_range,maxtime,run_options,dasar_coords);
        for istation=1:length(goodname{1})
            tagstr=goodname{1}{istation}(1:22);
            
            search_str=[tagstr '*morph*'] ;
            %automated_dir=[detectiondir '/site_0' Isite '/' keyword.scenario '/' keyword.algorithm];
            automated_dir=sprintf('%s/site_0%i/%s/%s',detectiondir,Isite,keyword.scenario,keyword.algorithm);
            disp(sprintf('loading automated detection results from %s',automated_dir));
            if isempty(automated_dir)
                disp('automated directory results not found in load_automated_results.m');
            end
            fname_auto=dir([automated_dir '/' search_str ]);
            if length(fname_auto)>1
                disp(sprintf('more than one automated detection file for %s',search_str));
                for jjj=1:length(fname_auto)
                    disp(fname_auto(jjj).name);
                end
                jjj=input('which file?');
                fname_auto=fname_auto(jjj);
            end
            
            %%%load automated date from file..
            fname_auto=[automated_dir '/' fname_auto.name];
            disp(['loading ' fname_auto]);
            stored_data=load(fname_auto);
            
            
            best_calls=stored_data.best_calls;
            debug_ctimes=stored_data.debug_ctimes;
            debug_durations=stored_data.debug_durations;
            raw_detections=stored_data.raw_detections;
            interval_detections=stored_data.interval_detections;
            param_stored=stored_data.param;
            %clear stored_data
            
            if ~isempty(best_calls.features)
                station_in(istation)=create_station(best_calls.equalization_freq,best_calls.equalization,best_calls.features, ...
                    best_calls.labeled_image, ...
                    debug_ctimes,debug_durations,param_stored.feature.names,param_stored.feature.global_names, ...
                    param_stored.feature.index ,param_stored.feature.Nsegments);
                
                
            end
            
            
        end  %end station...
        
        %%%%%%%%%%%%run neural network%%%%%%%%%%%%%%
        
        param.net=load_neural_network_path('neuralnetupdated');
        %param.net.name=param.net.name(1);
        
        if isempty(manual.individual)
            continue
        end
        ind=manual.individual{1};
        Nbase=size(ind.ctime,1);
        for J=1:size(ind.ctime,2)
            stats{Isite}.ninputs(Idate,J)=length(station_in(J).ctime_min);
        end
        for Ithresh1=1:length(thresholds1)
            for Ithresh2=1:length(thresholds2)
                disp(sprintf('Net 1 thresh: %6.2f, Net 2 thresh: %6.2f', ...
                    thresholds1(Ithresh1),thresholds2(Ithresh2)));
                station=filter_with_nnet(station_in,param.net.name,param.net.dir,[thresholds1(Ithresh1) thresholds2(Ithresh2)] );
                
                for J=1:size(ind.ctime,2)
                    Igood=find(ind.ctime(:,J)>0);
                    [Imiss,Ifalse,iman_pass,Itrue,Iredundant_auto]=find_similar_elements_ovlap_timelimit(station(J).ctime_min,station(J).Totalduration, ind.ctime(Igood,J), ...
                        ind.duration(Igood,J),mintime,maxtime,run_options.match_tol, [station(J).Totalfmin; station(J).Totalfmax],[ind.flo(Igood,J) ind.fhi(Igood,J)]');
                    
                    %%Note that length(Imiss)+length(Itrue)~=legth(station.ctime) because we strip away times not
                    %%within mintime/maxtime
                    stats{Isite}.nmatch(Idate,J,Ithresh1,Ithresh2)=length(Itrue);
                    stats{Isite}.nmiss(Idate,J,Ithresh1,Ithresh2)=length(Imiss);
                    stats{Isite}.nfalse(Idate,J,Ithresh1,Ithresh2)=length(Ifalse);
                    
                    
                end
            end
        end
        
        
        
    end
    save temp
    
end


save(sprintf('NeuralNet_stats_%s_%s',date_str{1}{1},date_str{2}{end}),'date_str','site_vector', ...
    'stats','thresholds1','thresholds2','run_options');

%%Plot a ROC curve.
total.nmatch=zeros(length(thresholds1),length(thresholds2));
total.nmiss=total.nmatch;
total.nfalse=total.nmatch;
total.ninputs=0;
for Isite=site_vector
    total.nmatch = total.nmatch + squeeze(sum(sum(stats{Isite}.nmatch,1),2));
    total.nmiss = total.nmiss + squeeze(sum(sum(stats{Isite}.nmiss,1),2));
    total.nfalse = total.nfalse + squeeze(sum(sum(stats{Isite}.nfalse,1),2));
    total.ninputs=total.ninputs+sum(sum(stats{Isite}.ninputs));
end

total.miss_fraction=total.nmiss./(total.nmiss+total.nmatch);
total.false_ratio=total.nfalse./(total.nmiss+total.nmatch);
for It1=1:length(thresholds1)
    for It2=1:length(thresholds2)
        disp(sprintf('Net 1 thresh: %6.2f, Net 2 thresh: %6.2f, miss percent: %6.2f, false ratio: %6.2f', ...
            thresholds1(It1),thresholds2(It2),100*total.miss_fraction(It1,It2), ...
        total.false_ratio(It1,It2)));
         
    end
end
% 
% figure
% subplot(2,1,1)
% plot(total.miss_fraction,total.false_ratio,'x');grid on;  %Plots first network ROC
% subplot(2,1,2)
% plot(total.miss_fraction(2,:),total.false_ratio(2,:),'x');grid on;  %Plots first network ROC
% 
% 
% figure
% subplot(2,1,1)
% [cc,hh]=contourf(thresholds1,thresholds2,total.miss_fraction,0:0.05:1);clabel(cc,hh);colorbar
% grid on;ylim([-1 -0.4])
% xlabel('threshold2');ylabel('threshold1');title('miss fraction');
% hold on; plot(0.8,-0.8,'yo')
% subplot(2,1,2)
% [cc,hh]=contourf(thresholds1,thresholds2,total.false_ratio,0:0.5:8);clabel(cc,hh);colorbar
% grid on;ylim([-1 -0.4])
% xlabel('threshold2');ylabel('threshold1');title('false ratio');
% hold on; plot(0.8,-0.8,'yo')


