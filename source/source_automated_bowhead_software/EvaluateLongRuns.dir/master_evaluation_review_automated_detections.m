
%%%master_evaluation_review_automated_detections.m%%
%% Review spectrograms of linked automated detections for quality spot-checking...
%% Aaron Thode
%% Jan 29, 2012

clear;close all

%%Load case-specific data
[Icasestr,date_str_local,Site_vector,params_chc, ...
    DASAR_str_local,rawdatadir_local,detectiondir_local,manualdir_local]=load_local_runparams(mfilename);



%%%Set plot_links to 1 in order to make maps and debug information..
run_options.plot_links=1;
run_options.find_DASAR_dates_rawdata=0;  %If one, get dates from raw data (not processed results).  To plot spectrograms
% need to set to 1.
time_inc=1; %Time interval covered by movie frame


for Isite=Site_vector  % For each site..
    Icase=sprintf(Icasestr,Isite);
    disp(sprintf('Processing %s, Site %i',Icase,Isite));
    
    %[rawdatadir,Icase,detectiondir,param,manualdir]=load_pathnames(Icase);
    master_setup;
    run_options.center_dist_limit=50000;
    
    if ~isempty(rawdatadir_local)
        rawdatadir=rawdatadir_local;
    end
    if ~isempty(detectiondir_local)
        detectiondir=detectiondir_local;
    end
    if ~isempty(manualdir_local)
        manualdir=manualdir_local;
    end
    
    
    for Idate=1:length(date_str{1})
        try
            tic
            date_str_new{1}{1}=date_str{1}{Idate};
            date_str_new{2}{1}=date_str{2}{Idate};
            
            [date_range,date_range_cell]=expand_date_range(date_str_new);
            
            DASAR_coords=load_DASAR_coords(Icase,Isite);
            %DASAR_coords=DASAR_coords(1:7,:);  %Only main array..
            
            
            %Load specific DASAR file names, used to access data.
            
            
            for Idd=1:size(date_range,1)
                [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates(date_range(Idd,:), ...
                    Isite,DASAR_str,rawdatadir, Icase);
            end
            
            %if run_options.reload_results==0,
            %    load temporary_data
            %else
            head=readgsif_header(goodFile{1}{1});
            mintime=head.ctbc;
            maxtime=head.ctbc+run_options.max_hrs*3600;
            
            [auto.locations,auto.locations_ctime,auto.stations,auto.raw_stations,automated_dir,auto_param,auto.goodName,auto.error_area]= ...
                load_automated_results(keyword,detectiondir, goodName,date_range,maxtime,run_options,DASAR_coords);
            
            %Conduct optional filtering of results
            if run_options.auto_location_filter==1
                filter_names={'Contour_global_bandwidth','Contour_fmax','Contour_fmin'};
                filter_criteria=[0 20 5; 300 400 500];
                [Ipass,feature_matrix]=crude_filter_feature_locations(auto.locations{1},filter_names,filter_criteria,run_options.min_stations(1),1);
                Ifail=setdiff(1:length(auto.locations{1}),Ipass);
                disp(sprintf('Filter feature activated, %i out of %i locations pass',length(Ipass), length(auto.locations{1})));
                
                
                auto.locations{1}=auto.locations{1}(Ipass);
                auto.locations_ctime{1}=auto.locations_ctime{1}(Ipass,:);
            end
            
            
            toc
        catch
            disp('Couldn''t load');
            continue
        end
        
        %master_load_manual_data;
        DASAR_coords=load_DASAR_coords(Icase,Isite,goodFile{1});
        
        %%%Permit an individual location manual or automated index to be selected.  Print out raw
        %%%information (time, frequency range, duration) for each selection,
        %%%  and plot spectrograms of signal for each location...
        
        if ~isempty(auto)
            
            %try
            movie_name= sprintf('all_AutomatedOnlySite%i_%s_%s',Isite,date_str{1}{Idate},keyword.stage);
            if run_options.auto_location_filter==1
                movie_name=[movie_name '_filtered']
                
            end
            %plot_movie_all_auto_only(movie_name,DASAR_coords,auto,Isite,run_options.center_dist_limit,param,goodFile);
            plot_movie_all_auto(DASAR_coords,[],auto,[],[],Isite,run_options.center_dist_limit,time_inc,goodFile);
            
            % [loc_index,auto_corrected]=plot_movie_all_auto(DASAR_coords,manual,auto,miss_stats,auto_stats,Isite,manual_limits,time_inc,goodFile)
            %
            %catch
            %  disp('Couldn''t plot');
            %end
            
        else  %No automated results loaded...
            disp('Automated localization results are empty');
            
            total.Ncompare_matches{Isite,Idate}=0;
            total.Nstation_total{Isite,Idate}=0;
            
            total.Nstation_match{Isite,Idate}=zeros(1,length(goodName{1}));
            total.Nstation_false{Isite,Idate}=zeros(1,length(goodName{1}));
            total.Nstation_auto{Isite,Idate}=total.Nstation_false{Isite,Idate};
            
            Ntemp=manual.individual{1}.ctime(:);
            total.Nstation_miss{Isite,Idate}=length(find(Ntemp>0));
            total.Nstation_manual{Isite,Idate}=total.Nstation_miss{Isite,Idate};
            
            
            
            total.goodName{Isite}=goodName{1};
            total.Nlongest_link{Isite,Idate}=0;
            total.Nsplit_links{Isite,Idate}=0;
            total.match_flag{Isite,Idate}=zeros(length(manual.localized{1}),1);
            total.auto_match_flag{Isite,Idate}=0;
        end
        
        
        
    end
    % clear auto manual
    
    
end




