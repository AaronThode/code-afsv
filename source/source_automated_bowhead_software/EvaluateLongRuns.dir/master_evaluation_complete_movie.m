%%%master_evaluation_complete_movie.m%%
%% Compare output of automatic detector with Greeneridge manual results
%% Aaron Thode
%%  Sept. 18, 2008
clear;close all
mydir=pwd;

[Icasestr,date_str_local,Sites,params_chc, ...
    DASAR_str_local,rawdatadir_local,detectiondir_local]=load_local_runparams(mfilename);

%time_inc=1;
time_inc=4;  %Four hour update
%time_inc=24;  %Four hour update

%%%%%Properties to view all sites
mapmidLon=-147;  % Set middle of latitude and longitude values
mapmidLat=70.7;
maplimits=[-141 69.85]; % Zoom level - this is right edge and bottom of map


%%%%%%Properties to view Site 5 (VLA) only
%mapmidLon=-143;  % Set middle of latitude and longitude values
%mapmidLat=70.33;
%maplimits=[-142 70+1/6]; % Zoom level - this is right edge and bottom of map

%%%%%%Properties to view Site 2 only
%mapmidLon=-149;  % Set middle of latitude and longitude values
%mapmidLat=70.75;
%maplimits=[-147 70.5]; % Zoom level - this is right edge and bottom of map




run_options.shooting_only=0;  %If 1, plot GPS positive for active seismic only
mask_locations.on=0;
mask_locations.lats= [71+1/6 71 70.66 70.5];
mask_locations.longs= -[151 148.5 146 141];
run_options.auto_only=1;  %Only show automated results

if findstr(Icasestr,'Shell08')>0
    seis=load('../ShipTracking/Seismic_tracks_2008_final.mat');
end

%goodDASAR_all=cell(length(Sites)
total_locs=0;
total_calls=0;
for Isite=Sites
    
    if Isite==6
        Isite_index=6;
        Isite=0;
       
    else
        Isite_index=Isite;
    end
    Icase=sprintf(Icasestr,Isite);
    master_setup;
    run_options.max_hrs=24;

    [Icasestr,date_str,Sites]=load_local_runparams(mfilename);
    
    keyword_all{Isite_index}=keyword;
    [date_range,date_range_cell]=expand_date_range(date_str);
    DASAR_coords{Isite_index}=load_DASAR_coords(Icase,Isite,[],'wantL');
    DASAR_coords_UTM{Isite_index}=load_DASAR_coords(Icase,Isite,[]);
    
    %if (DASAR_coords{Isite}(6,2)==0)
    %    DASAR_coords{Isite}(6,:)=DASAR_coords{Isite}(5,:);
    %end
    
    
    %Load specific DASAR file names, used to access data.
    
    
    try %%Attempt to locate Processed data locations first
       
        %%Try to locate processed results
        for Idd=1:size(date_range,1),
            [tmp1,tmp2,tmp3]=find_DASAR_dates_Processed(date_range(Idd,:), ...
                Isite,DASAR_str,detectiondir, Icase);
            
            if isempty(tmp1)
                [tmp1,tmp2,tmp3]=find_DASAR_dates(date_range(Idd,:),Isite,DASAR_str,rawdatadir, Icase);
                
            end
            fprintf('%i DASARS found in processed results for Site %i on %s\n',length(tmp1),Isite,date_range(Idd,:))

            goodDASAR_all{Isite_index,Idd}{1}=tmp1;
            goodFile_all{Isite_index,Idd}{1}=tmp2;
            goodName_all{Isite_index,Idd}{1}=tmp3;
        end
        
          
       
        
    catch  %look at raw data locations
        fprintf('Loading DASAR names from raw data, Site %i date %s\n', ...
            Isite,date_range(Idd,:));
        for Idd=1:size(date_range,1),
            [tmp1,tmp2,tmp3]=find_DASAR_dates(date_range(Idd,:),Isite,DASAR_str,rawdatadir, Icase);
             fprintf('%i DASARS found in processed results for Site %i on %s\n',length(tmp1),Isite,date_range(Idd,:))

            goodDASAR_all{Isite_index,Idd}{1}=tmp1;
            goodFile_all{Isite_index,Idd}{1}=tmp2;
            goodName_all{Isite_index,Idd}{1}=tmp3;
        end
    end
end
param.time_inc=time_inc;

Icount=0;clear FF

for Iguns=1:3
    past_active(Iguns).long=[];
    past_active(Iguns).lat=[];
    past_notactive(Iguns).long=[];
    past_notactive(Iguns).lat=[];
    
end

%%%Generate movie
if run_options.auto_only==0
    movie_name= sprintf('manVsauto_movie_%s',Icase);
else
    movie_name= sprintf('complete_movie_%s',Icase);
    
end

save movie_data_temp

%aviobj = avifile(movie_name);
aviobj=VideoWriter(movie_name,'Motion JPEG AVI');
open(aviobj);
%aviobj.Quality=50;

for Idd=1:size(date_range,1)
    tabs1=datenum([date_range(Idd,5:6) '/' date_range(Idd,7:8) '/' date_range(Idd,1:4)]);
    tabs_range=tabs1:datenum(0,0,0,param.time_inc,0,0):(tabs1+datenum(0,0,1,0,0,0));
    clear auto manual
    
    bad_day=zeros(1,max(Sites));
    for Isite=Sites
        
        if Isite==6
            Isite_index=6;
            Isite=0;
            
        else
            Isite_index=Isite;
        end
        Icase=sprintf(Icasestr,Isite);
        for JJJ=1:length(goodFile_all{Isite_index,Idd})
            if isempty(goodFile_all{Isite_index,Idd}{JJJ})
                bad_day(Isite_index)=bad_day(Isite_index)+1;
            end

        end
        if bad_day(Isite_index)>0
            fprintf('%i bad DASAR dates detected at Site %i on %s, skipping site...\n',bad_day(Isite_index),Isite,date_range(Idd,:));
            auto{Isite_index}.locations{1}=[];
            auto{Isite_index}.locations_ctime{1}=[];
            continue
        end
        %head=readgsif_header(goodFile_all{Isite,Idd}{1}{1});
        %maxtime=head.ctbc+run_options.max_hrs*3600;
        
        t=goodName_all{Isite_index,Idd}{1}{1}(8:end);
        mintime_tabs=datenum(str2num(t(1:4)),str2num(t(5:6)),str2num(t(7:8)),str2num(t(10:11)),str2num(t(12:13)),str2num(t(12:14)));
        mintime=24*3600*(mintime_tabs-datenum(1970,1,1,0,0,0));
        maxtime=mintime+run_options.max_hrs*3600;
        
        try
            [auto{Isite_index}.locations,auto{Isite_index}.locations_ctime,station_check]= ...
                load_automated_results(keyword_all{Isite_index},detectiondir, goodName_all{Isite_index,Idd},date_range(Idd,:),maxtime,run_options,DASAR_coords{Isite_index});
            disp(' ');
            total_locs=total_locs+length(auto{Isite_index}.locations{1});
            for IJK=1:size(auto{Isite_index}.locations_ctime{1},2)
               total_calls=total_calls+length(find(auto{Isite_index}.locations_ctime{1}(:,IJK)>0)); 
            end
            disp(sprintf('Cumulative total: %i',total_locs));
        catch
            disp(sprintf('!!!!! automated %s failed to load automated results',goodName_all{Isite_index,Idd}{1}{1}));
            cd(mydir)
        end
        
        if run_options.auto_only==0
            try
                
                [manual{Isite_index}.localized,manual{Isite_index}.individual]=load_manual_results_TSV(manualdir, ...
                    goodDASAR_all{Isite_index,Idd},date_range(Idd,:),maxtime,run_options,DASAR_coords_UTM{Isite_index},Isite);
                disp(sprintf('\n\n\n'));
            catch
                disp(sprintf('Manual results for %s failed to load',goodName_all{Isite_index,Idd}{1}{1}));
                % keyboard
                cd(mydir)
            end
        end
    end
    
    if all(bad_day>0)||~exist('auto')  %No results for any site on this day...
        fprintf('Skipping %s \n',(date_range(Idd,:)));
        continue
    end
    
    
    
    
    for Ihr=1:(length(tabs_range)-1)
        Icount=Icount+1;
        
        %Plot everything except whale locations
        Iplot=1;
        subroutine_plot_movie_background_and_seismics;
        Nlocs=0;
        for Isite=Sites  %for each site
            
            
            if Isite==6
                Isite_index=6;
                Isite=0;
                
            else
                Isite_index=Isite;
            end
            try
                %Plot DASAR coordinates
                m_line(DASAR_coords{Isite_index}(:,1),DASAR_coords{Isite_index}(:,2),'clip','point','marker','^','MarkerSize',8, 'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');

                %%Plot vertical array too!
                if ~isempty(findstr('Shell10',Icase))&&Isite==5
                    VA.lat=[70+26.848/60];
                    VA.lon=[-143-12.572/60];
                    
                    m_line(VA.lon,VA.lat,'clip','point','marker','d','MarkerSize',12, 'MarkerFaceColor','g','MarkerEdgeColor','g','LineStyle','none');
                    m_ruler(0.8,[0.6 0.9],'tickdir','out')  %Draw scale
                end
                Iplots=plot_movie_all_auto_withmap(tabs_range(Ihr+(0:1)),DASAR_coords{Isite_index},auto{Isite_index},Isite_index,mask_locations);
                Nlocs=Nlocs+ length(Iplots);
                title({'Automated results';title1;title2;sprintf('%i locations',Nlocs)});
            catch
                disp(sprintf('%s-%s failed for Site %i',datestr(tabs_range(Ihr)),datestr(tabs_range(Ihr+1)),Isite));
                cd(mydir)
                keyboard;
            end
            
        end
        
        if exist('manual')&&run_options.auto_only==0
            Iplot=2;
            subroutine_plot_movie_background_and_seismics;
            
            for Isite=Sites
                
                if Isite==6
                    Isite_index=6;
                    Isite=0;
                    
                else
                    Isite_index=Isite;
                end
                if Isite>length(manual)
                    continue
                end
                try
                    plot_movie_all_manual_withmap(tabs_range(Ihr+(0:1)),DASAR_coords{Isite_index},manual{Isite_index},Isite,mask_locations);
                catch
                    disp(sprintf('%s-%s failed',datestr(tabs_range(Ihr)),datestr(tabs_range(Ihr+1))));
                    cd(mydir)
                   % keyboard;
                end
                
            end
              title({'Manual results';title1;title2;sprintf('%i locations',Nlocs)});
      
        end
        figure(100);
        
        FF=getframe(gcf);
        %aviobj = addframe(aviobj,FF);
        writeVideo(aviobj,FF);
        
        %keyboard;
    end   %Ihr
      
end %Idd

close(aviobj);




%movie2avi(FF,movie_name,'fps',12);






