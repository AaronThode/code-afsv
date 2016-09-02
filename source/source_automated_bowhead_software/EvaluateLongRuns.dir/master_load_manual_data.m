%master_load_manual_data.m
% needs date_str,Idate,Site_vector,Isite,DASAR_str,rawdatadir,manualdir,Icase,
% run_options.max_hrs,

date_str_new{1}{1}=date_str{1}{Idate};
date_str_new{2}{1}=date_str{2}{Idate};

[date_range,date_range_cell]=expand_date_range(date_str_new);


%Load specific DASAR file names, used to access data.

try 
     %%Try to locate processed results
    for Idd=1:size(date_range,1)
        [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates_Processed(date_range(Idd,:), ...
            Isite,DASAR_str,detectiondir, Icase);
    end
    
    t=goodName{1}{1}(8:22);
    mintime_tabs=datenum(str2num(t(1:4)),str2num(t(5:6)),str2num(t(7:8)),str2num(t(10:11)),str2num(t(12:13)),str2num(t(12:14)));
    mintime=24*3600*(mintime_tabs-datenum(1970,1,1,0,0,0));
    maxtime=mintime+run_options.max_hrs*3600;
    
   
    
catch %%Attempt to locate raw data locations 
    
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
