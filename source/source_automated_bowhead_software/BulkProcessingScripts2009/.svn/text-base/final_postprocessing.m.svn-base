%%%%%all_processing.m%%%%% %  Conducts a bulk run of automated call detector to Greeneridge Data...
%   September 18, 2008 working version...


clear ;fclose('all');close all

path('../CommonScripts.dir',path);
path('../ComputerSpecificScripts.dir',path);
path('../BulkProcessingScripts',path);

if isempty(findstr(pwd,'Bulk'))
    error('Go to correct directory');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%CHANGE HERE...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Run options..

%keyword to select time and spatial subset of deployment...
%First field determines output folder and is used by TOC_runs to select
% appropriate site and company, which it then uses to select dates and
% locations.  For example, if first field has 'Shell08' and 'Site5'
% strings, then that site is selected.
%Second field determins subfolder
%Third field determines what string in final filename--prevents possible
%  writeovers while allowing individual DASAR files to be resused
%Site_vec=[1:2];
Site_vec=[1:5];
file_tag='NoDttCheck';
Icase_str='Shell08_Site%i_PeakBulkRunCore2.morph.Final';
params_chc='Shell08_March_optimized';
%params_chc='Shell08_Final_February_archive';
param=TOC_params(params_chc);  %Unpack detection parameters from keyword
%param.energy.nstart=23*60*60*1000;
%param.energy.nstart=17*60*60*1000+20*60*1000;
%param.energy.nstart=4*60*60*1000+23*60*1000;
%param.energy.nstart=5*60*60*1000+4.5*60*1000;  %Seals, site 1, 9/21/08
%param.energy.nsamples=0.25*60*60*1000;

%year_str='07';
run_options.debug.raw=0;    %If one, show debug energy detector information
run_options.debug.interval=0;
run_options.debug.morph=0;  %If one, show debug output for morph processing, and load debug_ctimes to learn when "failure" occures
run_options.debug.snips=0;  %If one, compare snips data with direct download of raw data...
run_options.debug.cross_channel=0;  %If 1, plot intermediate cross-channel output, if 2 plot more detailed input...
run_options.min_stations=2;  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria

run_options.max_time_delay=7;  %Maximum distance in seconds allowed between DASARS
run_options.dt_slop=0.5;  %How much tolerance to give time delays between DASARs
run_options.bearing_alg='sel_ratio';
run_options.plot_locations=0;  %If one, plot locations..
run_options.localization_alg='Huber'; %HuberCalibratedKappa, repeated_median, MedianHuber
run_options.filter_chc='min_weight'; %dt_error, min_weight
run_options.kappa_Nsamples=100;
run_options.auto_location_filter=1;  %Final desperation filter on location features to strip out seismics

run_options.interval_filter.on=0;
run_options.interval_filter.max_pass=0.05;
run_options.interval_filter.freq_tol=50;
run_options.interval_filter.time_span=datenum(0,0,0,4,0,0);
run_options.interval_filter.min_ratio=0.59

run_options.seismic_veto.on=0;
run_options.seismic_veto.std=10;  %Angular standard deviation in degrees.
run_options.cutoff_veto_distance=30000;  %distance in km beyond which non-active seismic source is not vetoed.
run_options.threshold='fixed';  %If 'fixed', use seismic_veto.std fixed angle for rejections threshold.
% Otherwise, 'adaptive' utilizes twice the standard deviation
% of the measured bearing for threshold.
run_options.debug.plot_seismic_veto=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%NEVER TOUCH BELOW UNDER NORMAL CIRCUMSTANCES%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seis=load('../ShipTracking/Seismic_tracks_2008_final.mat');

%Run a very rigid interval filter
if run_options.debug.interval==1;
    interval_debug.names=param.interval_remove.names;
    interval_debug.index=param.interval_remove.debug_index;
    if isfield(param.interval_remove,'debug_time')
        interval_debug.time=param.interval_remove.debug_time;
    end
else
    interval_debug=[];
end

now_str=datestr(now,30);


for Isite=Site_vec
    Icase=sprintf(Icase_str,Isite);
    homedir=pwd;

    [rawdatadir,outputdir,param]=load_pathnames(Icase,param);
    run_options.calculate_airgun_characteristics=0;  %Never change..

    [Site_vector,date_str,DASAR_str,keyword]=TOC_runs(Icase);
    run_options.filtering_stage=keyword.algorithm;  %'contour','morph','both','none': Determines what third stage processing takes place..
    
    date_str{1}={'20080812'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
    date_str{2}={'20081002'};  %End date
    
     %date_str{1}={'20080825'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
    %date_str{2}={'20080825'};  %End date 8/26 other check

   % date_str{1}={'20080822'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
   % date_str{2}={'20080901'};  %End date


    %%Expand date parameters into a string of contiguous dates.
    date_range=expand_date_range(date_str);

    
    cd(homedir)
    finaldir2=sprintf('%s/Site_0%i/',outputdir,Isite);
    finaldir1=sprintf('%s/Site_0%i/%s',outputdir,Isite,keyword.scenario);
    finaldir=sprintf('%s/Site_0%i/%s/%s',outputdir,Isite,keyword.scenario,keyword.algorithm);
    disp(sprintf('Output to be written to %s',finaldir));
    param.calibration_keyword=Icase;


    %%Loop through desired dates
    for Idate=1:size(date_range,1)
        cd(homedir);
        mydir=pwd;
        try
            fclose('all');clear station locations;
            %%For each date, identify all DASARS containing that date,
            %%restricted by DASAR_str
            [goodDASAR,goodFile,goodName,goodDASARstr]=find_DASAR_dates(date_range(Idate,:),Isite,DASAR_str,rawdatadir,Icase);
            if isempty(goodDASAR),
                disp(sprintf('%s was not assigned any DASARS at Site %i',date_range(Idate,:),Isite));
            end


            cd(finaldir);
            fname_out=sprintf('%s_%s_%s_FilteredLocations',goodName{end},keyword.stage,run_options.localization_alg);
            disp(fname_out);
            % keyboard
            %save(fname_out,'locations','locations_ctime','station','raw_station','param','goodName','goodFile','Icase','Isite','run_options');
            run_options_org=run_options;

            load(fname_out);
            run_options=run_options_org;
            cd(mydir);
            [goodDASAR,goodFile,goodName,goodDASARstr]=find_DASAR_dates(date_range(Idate,:),Isite,DASAR_str,rawdatadir,Icase);


            %%%Revised interval detector that uses bearing as a feature...




            %%%Filter out bearings
            Iship=1;
            tabs=zeros(length(station),length(locations));
            bearings=-1*ones(size(tabs));
            sd_bearings=bearings;
            %fmin=bearings;
            %fmax=bearings;
            for I=1:length(locations)
                ttstart1=clock;
                tabs(:,I)=datenum(1970,1,1,0,0,locations{I}.ctime_debug+param.energy.bufferTime);
                bearings(:,I)=locations{I}.bearing;
                sd_bearings(:,I)=locations{I}.bearing_sd;
                %fmin(:,I)=locations{I}.Totalfmin;
                %fmax(:,I)=locations{I}.Totalfmax;
                if isfield(locations{I},'position')
                    locations{I}=rmfield(locations{I},'position');
                end

            end
            %%Identify potential interval

            param.interval_remove.Ndet=20;
            param.interval_remove.Nmiss=5;
            param.interval_remove.ICItol=0.75;

            if run_options.interval_filter.on==1

                 for Istation=1:length(station)
                    clear test

                    disp(sprintf('station: %i, date: %s',Istation,datestr(tabs(Istation,I))));

                    Igood=find(tabs(Istation,:)>datenum(1971,1,1,0,0,0));

                    if run_options.debug.interval==1
                        interval_debug.fname=goodFile{Istation};
                        %interval_debug.time=datenum(2008,8,26,3,20,0);
                        interval_debug.index=880;
                        interval_debug.bearing=bearings(Istation,Igood);

                    end
                    [ICI,trelated]=compute_ici_bothways_bearing(tabs(Istation,Igood),datenum(1970,1,1,0,0,raw_station(Istation).raw_detections.ctime),[9 22],...
                        param.interval_remove.Ndet,param.interval_remove.Nmiss,...
                        param.interval_remove.ICItol,interval_debug);

                    Iint=find(ICI>0);
                    
                    if length(Iint)<=10
                        continue
                    end
                    test.ICI=ICI(Iint);
                    test.bearing=bearings(Istation,Igood(Iint));
                    test.sd=sd_bearings(Istation,Igood(Iint));
                    test.tabs=tabs(Istation,Igood(Iint));

                    for JJ=2:(length(Iint)-1)
                        index1=max([1 JJ-10]);
                        index2=min([length(Iint) (JJ+10)]);
                        index=setdiff(index1:index2,JJ);
                        index=[index1:(JJ-1) (JJ+1):index2];
                        test.medianb(JJ)=median(test.bearing(index));
                        test.std(JJ)=std(test.bearing(index));

                    end
                    test.medianb(length(Iint))=test.medianb(JJ);

                    Ikill=find(abs(test.medianb-test.bearing)<=run_options.seismic_veto.std);

                    Ifix=Igood(Iint(Ikill));
                    
                    for III=Ifix
                        %keyboard
                        locations{III}.bearing(Istation)=NaN;

                    end
                
                    if run_options.debug.interval==1
                        Ifigure=50;comment='bearing_count';
                        
                        figure(Ifigure-1);
                        subplot(3,1,1)
                        plot(bearings(Istation,Igood(Iint)),ICI(Iint),'x',bearings(Istation,Ifix),ICI(Iint(Ikill)),'rx');grid on
                        
                        subplot(3,1,2)
                        plot(bearings(Istation,Igood(Iint)),test.Icount,'x',bearings(Istation,Ifix),test.Icount(Ikill),'rx');xlabel('bearing');ylabel('Number of shared bearings');
                        hold on
                        plot(bearings(Istation,Igood(Iint)),run_options.interval_filter.max_pass*4*60*60./test.ICI.','g.');
                        hold off
                          % subplot(2,1,2)
                        % plot(test.Icount,test.fcount,'rx');xlabel('bearing match');ylabel('frequency match');
                        subplot(3,1,3)
                        plot(bearings(Istation,Igood(Iint)),test.fcount./test.Icount,'x',bearings(Istation,Ifix),test.fcount(Ikill)./test.Icount(Ikill),'rx');xlabel('bearing');ylabel('Ratio of matched freq/number bearings');


                        figure(Ifigure);
                        subplot(2,1,1)
                        plot(tabs(Istation,Igood),ICI,'go',tabs(Istation,Igood(Iint)),ICI(Iint),'bo',tabs(Istation,Ifix),ICI(Iint(Ikill)),'rx');
                        legend('Before ICI','After ICI','After bearing test','Location','best');
                        ylabel('interval(sec)');
                        datetick('x',14);
                        subplot(2,1,2)
                        plot(tabs(Istation,Igood),bearings(Istation,Igood),'go',tabs(Istation,Igood(Iint)),bearings(Istation,Igood(Iint)),'bo',tabs(Istation,Ifix),bearings(Istation,Ifix),'rx');
                        legend('Before ICI','After ICI','After bearing test','Location','best');

                        axis('ij');;ylabel('azimuth (deg)');
                        datetick('x',14);
                        
                        figure(Ifigure+1)
                        subplot(2,1,1)
                        plot(Iint,ICI(Iint),'bo',Iint(Ikill),ICI(Iint(Ikill)),'rx');axis('ij');ylabel('interval(sec)');axis('xy');xlabel('Iint');
                        subplot(2,1,2)
                        plot(Iint,bearings(Istation,Igood(Iint)),'bo',Iint(Ikill),bearings(Istation,Ifix),'rx');xlabel('Iint');ylabel('azimuth (deg)');
                        legend('After ICI','After bearing test','Location','best');
                        axis('ij')
                        
                        figure(Ifigure)
                        orient tall
                        print('-djpeg',sprintf('postprocessing_ICInBearing_review_%s_%s_pt1',goodName{Istation},comment));

                        figure(Ifigure+1)
                        orient tall
                        print('-djpeg',sprintf('postprocessing_ICInBearing_review_%s_%s_pt2',goodName{Istation},comment));


                        
                        pause;
                        close((Ifigure+(-1:1)));
                    end
                end
               

            end
            %disp(sprintf('total comp time: %6.2f',etime(clock,ttstart1)));



            %Vector processing
            disp('Starting GPS veto');
            for Istation=1:length(station)


                %Create station bearing vector
                Igood=find(tabs(Istation,:)>datenum(1971,1,1,0,0,0));
                bearing_seis=interp1(seis.Seismic(Iship).tabs,seis.Seismic(Iship).bearing{Isite,Istation},tabs(Istation,Igood));
                dist_seis=interp1(seis.Seismic(Iship).tabs,seis.Seismic(Iship).dist{Isite,Istation},tabs(Istation,Igood));
                guns_seis=interp1(seis.Seismic(Iship).tabs,seis.Seismic(Iship).guns,tabs(Istation,Igood));
                %OK, now compare bearings
                if strcmp(run_options.threshold,'fixed')
                    Ibad=find(abs(bearing_seis-bearings(Istation,Igood))<=run_options.seismic_veto.std);
                else
                    Ibad=find(abs(bearing_seis-bearings(Istation,Igood))<=3*sd_bearings(Istation,Igood));

                end
                %Of this subset, check if guns are off
                Ikeep=find(guns_seis(Ibad)==0&dist_seis(Ibad)>run_options.cutoff_veto_distance);
                Ireject=setdiff(1:length(Ibad),Ikeep);

                Ifix=Igood(Ibad(Ireject));


                if run_options.debug.plot_seismic_veto==1
                    plot(seis.Seismic(Iship).tabs,seis.Seismic(Iship).bearing{Isite,Istation},'ro');hold on; xlimm=xlim;
                    plot(tabs(Istation,Igood),bearings(Istation,Igood),'bx');
                    plot(tabs(Istation,Ifix),bearings(Istation,Ifix),'go');
                     datetick('x',14);
                
                end
               for I=Ifix
                    %keyboard
                    locations{I}.bearing(Istation)=NaN;

                end

            end

            %Recompute locations with removed bearings.
            locations=compute_position(locations,goodFile,param,Icase,Isite,run_options);
            
           
            [fname_tsv,Iwrite]=write_tsv(locations,goodName,[file_tag '_' now_str]);

            cd(finaldir);
            fname_out=sprintf('%s_%s_%s_FilteredLocations_%s',goodName{end},keyword.stage,run_options.localization_alg,file_tag);
            % keyboard
            save(fname_out,'locations','locations_ctime','station','raw_station','param','goodName','goodFile','Icase','Isite','run_options');
            cd(mydir);


        catch
            disp(sprintf('Crash! date_range(Idate,:) =%s failed',date_range(Idate,:)));
            errormsg=lasterror;
            disp(errormsg.message);
            %disp(length(errormsg.stack))
            if (length(errormsg.stack)>0)
                disp(errormsg.stack(1).file)
                disp(errormsg.stack(1).name)
                disp(errormsg.stack(1).line)
            end
            cd(mydir)

        end
        if exist('locations')
            clear locations
        end
    end %IDate
end %Isite


