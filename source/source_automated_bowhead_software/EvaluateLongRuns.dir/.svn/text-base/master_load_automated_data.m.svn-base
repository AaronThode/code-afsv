%%%%%%master_load_automated_data.m

date_str_new{1}{1}=date_str{1}{Idate};
date_str_new{2}{1}=date_str{2}{Idate};

[date_range,date_range_cell]=expand_date_range(date_str_new);

DASAR_coords=load_DASAR_coords(Icase,Isite);
DASAR_coords=DASAR_coords(1:7,:);  %Only main array..
if all(DASAR_coords(6,:)==0)
    DASAR_coords(6,:)=DASAR_coords(5,:);
end


%Load specific DASAR file names, used to access data.

try %%Attempt to locate Processed data locations first
    
    %%Try to locate processed results
    for Idd=1:size(date_range,1),
        [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates_Processed(date_range(Idd,:), ...
            Isite,DASAR_str,detectiondir, Icase);
    end
    
    t=goodName{1}{1}(8:22);
    mintime_tabs=datenum(str2num(t(1:4)),str2num(t(5:6)),str2num(t(7:8)),str2num(t(10:11)),str2num(t(12:13)),str2num(t(12:14)));
    mintime=24*3600*(mintime_tabs-datenum(1970,1,1,0,0,0));
    maxtime=mintime+run_options.max_hrs*3600;
    
    
    
catch  %look at raw data locations
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
end

[auto.locations,auto.locations_ctime,auto.stations,auto.raw_stations,automated_dir,auto_param,auto.error_area]= ...
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