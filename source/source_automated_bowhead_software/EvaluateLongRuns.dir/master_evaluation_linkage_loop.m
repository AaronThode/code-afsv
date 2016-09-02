
%%%master_evaluation_linkage_loop.m%%
%% Compare output of automatic detector with Greeneridge manual results
%% Compare both individual DASAR detection performance, along with linkages.
%%      For the latter a minimum number of DASARS must share the same automated and manual detection
%% Aaron Thode
%%  Jan 21, 2010
%%
%% Updated Feb 2, 2011
clear;close all

%%Load case-specific data
[Icasestr,date_str_local,Site_vector,params_chc,DASAR_str_local]=load_local_runparams(mfilename);


%%%Set plot_links to 1 in order to make maps and debug information..
run_options.plot_links=1;
run_options.find_DASAR_dates_rawdata=0;  %If one, get dates from raw data (not processed results).  To plot spectrograms 
                                            % need to set to 1.


for Isite=Site_vector  % Fore each site..
    Icase=sprintf(Icasestr,Isite);
    disp(sprintf('Processing %s, Site %i',Icase,Isite));
    master_setup;
    run_options.all_biologics=2;
    %run_options.strict_eval_tolerance:  If one, manual and auto detections must
    %               overlap 50% in time AND frequency.  If two, compare manual
    %               and auto detections via time tolerances only, with no overlap
    %               requirement.  If three, try to match
    %               by location...
    run_options.strict_eval_tolerance=1;
    run_options.eval_tol=2;  %Time uncertainty estimate for manual selection.
    run_options.match_tol=0.5;  %Fractional overlap in time (and frequency) required to declare a match.
    %[Icasestr,date_str]=load_local_runparams(mfilename);
    
    for Idate=1:length(date_str{1})
        try
            tic
            date_str_new{1}{1}=date_str{1}{Idate};
            date_str_new{2}{1}=date_str{2}{Idate};
            
            [date_range,date_range_cell]=expand_date_range(date_str_new);
            
            DASAR_coords=load_DASAR_coords(Icase,Isite);
            DASAR_coords=DASAR_coords(1:7,:);  %Only main array..
            if all(DASAR_coords(6,:)==0)
                DASAR_coords(6,:)=DASAR_coords(5,:);
            end
            
            
            %Load specific DASAR file names, used to access data.
            
            if run_options.find_DASAR_dates_rawdata==0 %%Attempt to locate Processed data locations first
                
                %%Try to locate processed results
                for Idd=1:size(date_range,1)
                    [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates_Processed(date_range(Idd,:), ...
                        Isite,DASAR_str,detectiondir, Icase);
                end
                
                t=goodName{1}{1}(8:22);
                mintime_tabs=datenum(str2num(t(1:4)),str2num(t(5:6)),str2num(t(7:8)),str2num(t(10:11)),str2num(t(12:13)),str2num(t(12:14)));
                mintime=24*3600*(mintime_tabs-datenum(1970,1,1,0,0,0));
                maxtime=mintime+run_options.max_hrs*3600;
                
                
                
            else  %look at raw data locations
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
            end
            
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
            
            if run_options.all_biologics==1
                disp('Loading all biological signals into manual');
                [manual.localized,manual.individual,Istrip]=load_manual_results_TSV(manualdir, ...
                    goodDASAR,date_range,maxtime,run_options,DASAR_coords,Isite);
            elseif run_options.all_biologics==0
                disp('Loading pinniped signals into separate variable');
                [manual.localized,manual.individual,Istrip,pinniped.localized,pinniped.individual]=load_manual_results_TSV(manualdir, ...
                    goodDASAR,date_range,maxtime,run_options,DASAR_coords,Isite);
                
            elseif run_options.all_biologics==2
                disp('Loading all whale signals into manual, 2008 default');
                [manual.localized,manual.individual,Istrip]=load_manual_results_TSV(manualdir, ...
                    goodDASAR,date_range,maxtime,run_options,DASAR_coords,Isite);
                
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
        
        if ~isempty(manual.localized)&&~isempty(manual.localized{1}.ctev)&&~isempty(auto.locations_ctime)&&~isempty(auto.locations_ctime{1})
            tic
            %[Nstation_linked,Ncompare_matches,Nstation_match,Nlongest_link,Nsplit_links,Nlocs_linked]
            [miss_stats,auto_stats]=evaluate_linking(param,manual.individual{1},auto.locations{1}, auto.locations_ctime{1},run_options,auto.stations{1});
            
            %[miss_stats,auto_stats]=evaluate_linking(param,manual.individual{1},auto.locations{1}, auto.locations_ctime{1},run_options);
            toc
            
            
            Nmanuals(Isite,Idate)=length(manual.localized{1}.ctev);
            %%Summarize statistics
            %disp(sprintf('Out of %i manual detections, %i had more than 2 detections in stations',Nmanuals(Isite,Idate),length(find(miss_stats.Nstation_match>=2))));
            disp(sprintf('Out of %i manual detections, %i had more than 2 detections in locations',Nmanuals(Isite,Idate),length(find(miss_stats.Nhits>=2))));
            disp(sprintf('Out of %i manual detections, %i had more than 2 detections in a single linkage',Nmanuals(Isite,Idate),length(find(miss_stats.Nlongest_link>=2))));
            disp(sprintf('Out of %i manual detections, %i had components in multiple automated locations',Nmanuals(Isite,Idate),length(find(miss_stats.Nsplit_links>=2))));
            
            %total.Nmanual_sites{Isite,Idate}=Nsites_manual;
            total.Ncompare_matches{Isite,Idate}=miss_stats.Ncompare_matches;
            %total.Nstation_match{Isite,Idate}=miss_stats.Nstation_match;
            total.Nstation_total{Isite,Idate}=miss_stats.Nhits;
            %             for II=1:size(miss_stats.Nstations_index_match,2)
            %                 total.Nstation_match{Isite,Idate}(II)=length(find(miss_stats.Nstations_index_match(:,II)>0));
            %                 total.Nstation_manual{Isite,Idate}(II)=length(find(miss_stats.Nstations_index_match(:,II)>=0));
            %                 total.Nstation_auto{Isite,Idate}(II)=length(find(auto.stations{1}(II).ctime_min>0));
            %             end
            
            total.Nstation_match{Isite,Idate}=miss_stats.Nstation_match;
            total.Nstation_manual{Isite,Idate}=miss_stats.Nstation_manual;
            total.Nstation_false{Isite,Idate}=miss_stats.Nstation_false;
            total.Nstation_miss{Isite,Idate}=miss_stats.Nstation_miss;
            total.Nstation_auto{Isite,Idate}=(miss_stats.Nstation_match)+(miss_stats.Nstation_false);
            
            total.goodName{Isite}=goodName{1};
            total.Nlongest_link{Isite,Idate}=miss_stats.Nlongest_link;
            total.Nsplit_links{Isite,Idate}=miss_stats.Nsplit_links;
            total.match_flag{Isite,Idate}=miss_stats.match_flag;
            total.auto_match_flag{Isite,Idate}=auto_stats.match_flag;
            
            %if ~isempty(intersect(Idate,[2 4 5]))
            if run_options.plot_links==1
                try
                    %plot_linkage_comparison(sprintf('linkmapSite%i_%s_%s',Isite,date_str{1}{Idate},keyword.stage),DASAR_coords, ...
                    %    miss_stats.Nlongest_link,manual,auto,miss_stats.Ncompare_matches,miss_stats.Ilocs_best_match,Isite);
                    movie_name= sprintf('all_automatedSite%i_%s_%s',Isite,date_str{1}{Idate},keyword.stage);
                    if run_options.auto_location_filter==1
                        movie_name=[movie_name '_filtered'];
                        
                    end
                    
                    plot_movie_all_auto(movie_name,DASAR_coords,manual,auto,miss_stats,auto_stats,Isite,run_options.center_dist_limit,param,goodFile);
                catch
                    disp('Couldn''t plot');
                end
            end
            
        elseif isempty(manual.localized)||isempty(manual.localized{1}.ctev)  %If no manual data available
            Nmanuals(Isite,Idate)=0;
            total.Ncompare_matches{Isite,Idate}=0;
            total.Nstation_total{Isite,Idate}=0;
            
            total.Nstation_match{Isite,Idate}=zeros(1,length(goodName{1}));
            total.Nstation_manual{Isite,Idate}=zeros(1,length(goodName{1}));
            total.Nstation_miss{Isite,Idate}=zeros(1,length(goodName{1}));
            
            total.Nstation_false{Isite,Idate}=zeros(1,length(goodName{1}));
            for IJ=1:length(auto.stations{1})
                total.Nstation_false{Isite,Idate}(IJ)=length(auto.stations{1}(IJ).ctime_min);
            end
            total.Nstation_auto{Isite,Idate}=total.Nstation_false{Isite,Idate};
            
            
            total.goodName{Isite}=goodName{1};
            total.Nlongest_link{Isite,Idate}=0;
            total.Nsplit_links{Isite,Idate}=0;
            total.match_flag{Isite,Idate}=0;
            total.auto_match_flag{Isite,Idate}=zeros(length(auto.locations{1}),1);
            
            if run_options.plot_links==1
                try
                    %plot_linkage_comparison(sprintf('linkmapSite%i_%s_%s',Isite,date_str{1}{Idate},keyword.stage),DASAR_coords, ...
                    %    miss_stats.Nlongest_link,manual,auto,miss_stats.Ncompare_matches,miss_stats.Ilocs_best_match,Isite);
                    % plot_movie_manual_match('manual_comparison',DASAR_coords,miss_stats.Nlongest_link,manual,auto,miss_stats.Ncompare_matches,miss_stats.Ilocs_best_match);
                    movie_name= sprintf('all_AutomatedOnlySite%i_%s_%s',Isite,date_str{1}{Idate},keyword.stage);
                    if run_options.auto_location_filter==1
                        movie_name=[movie_name '_filtered']
                        
                    end
                    %plot_movie_all_auto_only(movie_name,DASAR_coords,auto,Isite,run_options.center_dist_limit,param,goodFile);
                    plot_movie_all_auto(movie_name,DASAR_coords,manual,auto,[],[],Isite,run_options.center_dist_limit,param,goodFile);
                    
                catch
                    disp('Couldn''t plot');
                end
            end
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

%%%%%%Sumarize and save results%%%%%%%%%%
for Isite=Site_vector
    %Icase=sprintf(Icasestr,Isite);
    % master_setup;
    total.Nmanuals(Isite)=0;
    total.Nstation_matches(Isite)=0;
    total.Nlocation_matches(Isite)=0;
    total.Nlongest_links(Isite)=0;
    total.manual_pass(Isite)=0;
    total.auto_pass(Isite)=0;
    total.manual_miss(Isite)=0;
    total.auto_extra(Isite)=0;
    Nstations=length(total.Nstation_match{Isite,1});
    total.manual_station_pass{Isite}=zeros(1,Nstations);
    total.manual_station_miss{Isite}=zeros(1,Nstations);
    total.auto_station_extra{Isite}=zeros(1,Nstations);
    total.manual_station_total{Isite}=zeros(1,Nstations);
    
    try
        for Idate=1:length(date_str{1})
            %for Idate=2:2
            %a=length(find(total.Nstation_match{Isite,Idate}>=2));
            b=length(find(total.Ncompare_matches{Isite,Idate}>=2));
            c=length(find(total.Nlongest_link{Isite,Idate}>=2));
            
            manual_pass(Isite,Idate)=sum(total.match_flag{Isite,Idate});
            manual_miss(Isite,Idate)=length(total.match_flag{Isite,Idate})-manual_pass(Isite,Idate);
            
            auto_pass(Isite,Idate)=sum(total.auto_match_flag{Isite,Idate});
            auto_extra(Isite,Idate)=length(total.auto_match_flag{Isite,Idate})-auto_pass(Isite,Idate);
            
            
            total.manual_pass(Isite)=total.manual_pass(Isite)+manual_pass(Isite,Idate);
            total.manual_miss(Isite)=total.manual_miss(Isite)+manual_miss(Isite,Idate);
            total.auto_pass(Isite)=total.auto_pass(Isite)+auto_pass(Isite,Idate);
            total.auto_extra(Isite)=total.auto_extra(Isite)+auto_extra(Isite,Idate);
            
            
            %total.Nmanuals(Isite)=total.Nmanuals(Isite)+Nmanuals(Isite,Idate);
            %total.Nstation_matches(Isite)=total.Nstation_matches(Isite)+a;
            total.Nlocation_matches(Isite)=total.Nlocation_matches(Isite)+b;
            total.Nlongest_links(Isite)=total.Nlongest_links(Isite)+c;
            
            
            %                 total.Nstation_match{Isite,Idate}=miss_stats.Nstation_match;
            %                 total.Nstation_manual{Isite,Idate}=miss_stats.Nstation_manual;
            %                 total.Nstation_false{Isite,Idate}=miss_stats.Nstation_false;
            %                 total.Nstation_miss{Isite,Idate}=miss_stats.Nstation_miss;
            %                 total.Nstation_auto{Isite,Idate}=(miss_stats.Nstation_match)+(miss_stats.Nstation_false);
            %
            
            manual_station_pass=total.Nstation_match{Isite,Idate};
            manual_station_miss=total.Nstation_miss{Isite,Idate};
            manual_station_total=total.Nstation_manual{Isite,Idate};
            auto_station_extra=total.Nstation_false{Isite,Idate};
            
            total.manual_station_total{Isite}=total.manual_station_total{Isite}+manual_station_total;
            total.manual_station_pass{Isite}=total.manual_station_pass{Isite}+manual_station_pass;
            total.manual_station_miss{Isite}=total.manual_station_miss{Isite}+manual_station_miss;
            total.auto_station_extra{Isite}=total.auto_station_extra{Isite}+auto_station_extra;
            
            disp(sprintf('\n%%%%%%%%%%%\n\nSite %i, date %s:',Isite,date_str{1}{Idate}));
            %disp(sprintf('Out of %i manual detections, %i had 2 or more detections in locations',Nmanuals(Isite,Idate),b));
            %disp(sprintf('Out of %i manual detections, %i had 2 or more detections in a single linkage',Nmanuals(Isite,Idate),c));
            diary on
            for II=1:Nstations
                disp(sprintf('Site %i%s, date %s: %i matched manual stations, %i missed manual stations, %i extra auto stations\n', ...
                    Isite,total.goodName{Isite}{II}(5),date_str{1}{Idate},manual_station_pass(II),manual_station_miss(II),auto_station_extra(II)));
            end
            disp(sprintf('Site %i, date %s: %i matched manual locations, %i missed manual locations, %i matched auto locations, %i extra auto locations\n\n%%%%%%%%%%%\n', ...
                Isite,date_str{1}{Idate},manual_pass(Isite,Idate),manual_miss(Isite,Idate),auto_pass(Isite,Idate),auto_extra(Isite,Idate)));
            diary off
        end
    catch
        disp('master_evaluation_linkage_loop line 324');
        lasterror
        %keyboard;
    end
end

diary on

total.all_manual_station_miss=0;
total.all_manual_station_total=0;
total.all_auto_station_extra=0;

for Isite=Site_vector
    
    %disp(sprintf('Manual calls: %6.2f, into linking stage %6.2f, out of linking stage %6.2f',total.Nmanuals(Isite),total.Nstation_matches(Isite),total.Nlongest_links(Isite)));
    %disp(sprintf('%6.4f percent pass before linking, %6.2f pass after linking',100*total.Nstation_matches(Isite)/total.Nmanuals(Isite),100*total.Nlongest_links(Isite)/total.Nmanuals(Isite)));
    
    
    fprintf('**************\n');
    Nstations=length(total.Nstation_match{Isite,1});
    
    for II=1:Nstations
        miss_rate=total.manual_station_miss{Isite}(II)./(total.manual_station_total{Isite}(II));
        false_ratio=total.auto_station_extra{Isite}(II)/(total.manual_station_total{Isite}(II));
        
        fprintf('Site %i%s summary: %6.2f miss rate, %6.2f false ratio, %i matched manual stations, %i missed manual stations,  %i extra auto stations\n', ...
            Isite,total.goodName{Isite}{II}(5),miss_rate*100, false_ratio, total.manual_station_pass{Isite}(II), ...
            total.manual_station_miss{Isite}(II),total.auto_station_extra{Isite}(II));
    end
    
    station_miss_rate=sum(total.manual_station_miss{Isite})./sum(total.manual_station_pass{Isite}+total.manual_station_miss{Isite});
    station_false_ratio=sum(total.auto_station_extra{Isite})/sum(total.manual_station_pass{Isite}+total.manual_station_miss{Isite});
    
    miss_rate=total.manual_miss(Isite)./(total.manual_pass(Isite)+total.manual_miss(Isite));
    false_ratio=total.auto_extra(Isite)/(total.manual_pass(Isite)+total.manual_miss(Isite));
    
    
    fprintf(' Site %i station summary: %6.2f miss rate station, %6.2f false ratio station\n', ...
        Isite,100*station_miss_rate,station_false_ratio);
    
    fprintf(' Site %i summary: %6.2f miss rate, %6.2f false ratio, %i matched manual locations, %i missed manual locations, %i matched auto locations, %i extra auto locations\n*********\n', ...
        Isite,100*miss_rate,false_ratio,total.manual_pass(Isite),total.manual_miss(Isite),total.auto_pass(Isite),total.auto_extra(Isite));
    
    total.all_manual_station_miss=total.all_manual_station_miss+sum(total.manual_station_miss{Isite});
    total.all_manual_station_total=total.all_manual_station_total+sum(total.manual_station_total{Isite});
    total.all_auto_station_extra=total.all_auto_station_extra+sum(total.auto_station_extra{Isite});
    
    
    
end

miss_rate=sum(total.manual_miss)./(sum(total.manual_pass)+sum(total.manual_miss));
false_ratio=sum(total.auto_extra)./sum(total.manual_pass+total.manual_miss);


miss_rate_station=total.all_manual_station_miss./(total.all_manual_station_total);
false_ratio_station=total.all_auto_station_extra/(total.all_manual_station_total);

fprintf(' Total summary: %6.2f miss rate station, %6.2f false ratio station, %6.2f miss rate localization, %6.2f false ratio localaization', ...
    100*miss_rate_station,false_ratio_station,100*miss_rate,false_ratio);

run_options
diary off
%total.DASARstr=manual.individual{1}.DASARstr;  %Make sure we know what DASARS are available...
save(sprintf('Linkage_stats_%s_%s_Option%i',date_str{1}{1},date_str{2}{end},run_options.strict_eval_tolerance),'date_str','Site_vector', ...
    'manual_pass','manual_miss','auto_pass','auto_extra','total','run_options');


