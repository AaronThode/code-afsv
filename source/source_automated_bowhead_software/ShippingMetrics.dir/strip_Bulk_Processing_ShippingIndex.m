%%strip_Bulk_Processing_ShippingIndex.m%%
%%  Aaron Thode
%%  April 5, 2013
%%  Load shipping metrics and then use to remove detections from final localizations

function strip_Bulk_Processing_ShippingIndex

persistent ship MYSITE %Once have loaded huge variable, don't do again

close all

%%location of raw acoustic data (needed for plotting debug)
%  /Users/thode/Jonah/Data/Shell_Automated_Results/Arctic_2012/Processed.dir/Site_04
rawdatadir='~/Jonah/Data/Shell2012_GSI_Data/';
%rawdatadir='/Volumes/ThodePortable/Shell2012_GSI_Data/';
%rawdatadir='/Volumes/Field_Backup/Shell2012/Shell2012_GSI_Data/';

%%Location of automated results
results.topdir='~/Jonah/Data/Shell_Automated_Results/Arctic_2012/Processed.dir';
%results.topdir='/Volumes/ThodePortable/Shell_Automated_Results/Arctic_2012/Processed.dir';
%results.topdir='/Volumes/Field_Backup/Shell2012/Arctic_2012/Processed.dir/';

Site_vec=[ 4];  %Set to 0 for site 0
%DASAR_str='*';
date_str{1}={'20120718'};  %
date_str{2}={'20121118'};  %End date

Icase_str='Shell12_Site%i_InitialRun';  %keyword for particular automated run

Npeaks=1;
Ipeak=1;  %Choice of peak detector
chc='avg';  %Choice of consolidating across time window: med or avg


time_window=2; %Minutes

debug.plot=0;  %If one, make lots of debug plots
debug.flag=0;
debug.tabs=datenum(2012,9,18,1,16,2); %16:02:50
debug.dasar=7;

[date_range, date_range_cell]=expand_date_range(date_str);

for I=Site_vec
    
    if ~exist('MYSITE','var')||isempty(MYSITE)||MYSITE~=I
        disp('reloading data');
        %if I==3
        %    ship=load(sprintf('Shell2012_shipping_metrics_combined_0929_Sites3n4'));  %loads metric; param
        %
        %else
        ship=load(sprintf('ShippingMetrics.dir/Shell2012_shipping_metrics_basic_Site%i.mat',I));  %loads metric; param
        %end
        MYSITE=I;
    end
    
    %convert time window into metric index
    tabs=ship.metric{end}{end}.tabs(1:2);
    tabs=datevec(diff(tabs));
    time_window=ceil(time_window*60/tabs(end));
    
    switch I
        case 3
            Filt_str='Nov2012_Huber_FilteredLocations';  % search str for automated results.
        case {0,4}
            Filt_str='FilteredLocations_Huber_FilteredLocations';
        otherwise
            Filt_str='Nov2012_Huber_FilteredLocations';  % search str for automated results.
            
    end
    
    %%%Some previously-tested rules...
    %rule='(dfmid<5 | dfmin<2 | dfmax<2) && (Fmax-Fmin)<30';
    %rule='dfmin<10 | abs(Fmin-144)<15 ';  %144 Hz has a broadband drilling noise activity...
    %rule='abs(Fmin-144)<15';  %A test rule to check that locations are indeed removed...
    %rule='(dfmin<15 || abs(Fmin-135)<10) & Fmax-Fmin<20  ';  %144 Hz has a broadband drilling noise activity...
    %rule='(dfmin<15 || abs(Fmin-135)<10 || abs(Fmin-100)<10) & Fmax-Fmin<20  '; %Most restrictive rule
    %rule='((dfmin<10 ) & Fmax-Fmin<20) || solidity>0.95 ';
    %rule='((dfmin<10 ||dfmax<10) & Fmax-Fmin<20) || solidity>0.95 ';
    %rule='((dfmin<10 ||dfmax<10) & Fmax-Fmin<25) || solidity>0.95 ';
    
    if I==4
        rule='((dfmin<10 ) & Fmax-Fmin<25) || solidity>0.99 ';
    else
        rule='(dfmin<10 ) & Fmax-Fmin<25  ';
    end
    
    
    %Load automated data
    Icase=sprintf(Icase_str,I);
    
    for Idd=1:size(date_range,1)
        fprintf('Date %i checking\n',Idd);
        [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates(date_range(Idd,:), ...
            I,'*',rawdatadir, Icase);
    end
    
    DASAR_coords=load_DASAR_coords(Icase,I,goodFile{1});
    
    %Cycle through dates
    for Idate=1:size(date_range,1)
        
        if I==6
            I=0;
        end
        search_dir=sprintf('%s/Site_%02i/%s/morph/',results.topdir,I,Icase);
        fnames=dir(sprintf('%s/*%s*%s*mat',search_dir,date_range(Idate,:),Filt_str));
        
        if length(fnames)>1
            keyboard
        elseif isempty(fnames)
            continue
        end
        
        
        fprintf('Loading %s...\n',fnames.name);
        data=load(sprintf('%s/%s',search_dir,fnames.name));
        
        
        %%Initialize temporary variables that will mark removed localizations and detections
        locations_ctime=data.locations_ctime;
        %locations_ctime_combined=data.locations_ctime;
        locations_dfmin=1000*ones(size(locations_ctime));
        locations_dfmax=1000*ones(size(locations_ctime));
        locations_fmin=1000*ones(size(locations_ctime));
        locations_fmax=locations_fmin;
        locations_solidity=locations_fmin;
        locations_tabs=datenum(1970,1,1,0,0,data.locations_ctime);  %datenumbers of locations
        %%Make a copy of station_ctimes so we know what to strip
        for II=1:length(data.station)
            station_tabs_min{II}=datenum(1970,1,1,0,0,data.station(II).ctime_min);
        end
        
        %%Construct output names
        [pathstr,namme,extt] = fileparts(fnames.name);
        fname_out=[namme '_industrial_postprocessed' extt];
        fname_out_tsv=[namme '_industrial_postprocessed' ];
        fname_industrial=[namme '_industrial_rejected' extt];
        fname_industrial_tsv=[namme '_industrial_rejected' ];
        
        for Idasar=1:length(data.goodName)
            %for Idasar=6:6
            
            
            DASAR_str=data.goodName{Idasar}(5);
            J=double(lower(DASAR_str))-double('a')+1;  %Index for metrics
            
            
            if I==0
                I=6;
            end
            tabs=ship.metric{I}{J}.tabs;
            tabs_combined=ship.metric{I}{end}.tabs;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Decompose metric into more convenient variables %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %ctime=86400*(tabs-datenum(1970,1,1,0,0,0));
            %entropy=ship.metric{I}{J}.entropy;
            %kurtosis=ship.metric{I}{J}.kurtosis;
            for III=1:Npeaks
                
                %Erase previous values...
                %isi_Nf{III}=[];
                adp_Nf{III}=[];
                %isi_Nf_combined{III}=[];
                adp_Nf_combined{III}=[];
                
                %isi_pwr{III}= ship.metric{I}{J}.isi{III}.(chc).pwr;
                %isi_F{III}= ship.metric{I}{J}.isi{III}.(chc).Fpeaks;
                %isi_Nf{III}=sum(isi_F{III}>0);
                
                adp_pwr{III}= ship.metric{I}{J}.adp{III}.(chc).pwr;
                adp_F{III}= ship.metric{I}{J}.adp{III}.(chc).Fpeaks;
                adp_Nf{III}=sum(adp_F{III}>0);
                
                %isi_F_combined{III}= ship.metric{I}{end}.isi{III}.(chc).Fpeaks;
                %isi_Nf_combined{III}=sum(isi_F_combined{III}>0);
                
                adp_F_combined{III}= ship.metric{I}{end}.adp{III}.(chc).Fpeaks;
                adp_Nf_combined{III}=sum(adp_F_combined{III}>0);
                
            end
            
            
            %%Run through each detection time in localization object
            if ~isempty(data.locations_ctime)
                strip_industrial_from_locations;
            end
            
            %%Run through each detection time in station object
            if ~isempty( station_tabs_min{Idasar})
                strip_industrial_from_stations;
            end
        end %Idasar...
        
        if debug.plot==1
            strip_debug_plotting;
        else
            %%What to keep and what to retain?
            Igood=find(sum(locations_ctime'>0)>1);
            Ibad=find(sum(locations_ctime'>0)<2);
            clear Igood_station Ibad_station
            for Istation=1:length(data.station)
                Igood_station{Istation}=find(station_tabs_min{Istation}>0);
                Ibad_station{Istation}=find(station_tabs_min{Istation}<0);
            end
            write_locations_to_file(Igood,Igood_station, fname_out, fname_out_tsv);
            
            write_locations_to_file(Ibad,Ibad_station, fname_industrial,fname_industrial_tsv);
            
        end  %if debug plot
        
    end  %Idate
end  %I, Isite

    function fname_tsv=write_locations_to_file(Igood,Igood_station,fname_out,fname_out_tsv)
        if ~isempty(data.locations)
            locations=data.locations(Igood);
            locations_ctime=data.locations_ctime(Igood,:);
            linknames=fieldnames(data.linking_index);
            for Ilink=1:length(linknames)
                linking_index.(linknames{Ilink})=data.linking_index.(linknames{Ilink})(:,Igood);
            end
            
            station=data.station;
            for Istation=1:length(station)
                tmp=trim_station(station(Istation),Igood_station{Istation});
                station(Istation)=tmp;
            end
            
            try
                for K=1:length(station)
                    disp(sprintf('Processing singletons on station %i',K));
                    if ~isempty(station(K).indicies)
                        
                        %%Remove detections that have been used in localizations
                        Isingleton=setdiff(1:length(station(K).ctime_min),linking_index.I(K,:));
                        
                        %%Only keep station detections that are NOT localized
                        tmp=trim_station(station(K),Isingleton);
                        
                        Ipass=crude_filter_feature_stations(tmp,data.param.final_filter.names,data.param.final_filter.criteria);
                        disp(sprintf('%i out of %i singletons in station %i pass',length(Ipass),length(Isingleton),K));
                        station_single(K)=trim_station(tmp,Ipass);
                        % station_single(K)=compute_bearings_station(tmp,goodFile{K},param,run_options.bearing_alg,0,run_options.kappa_Nsamples);
                    else
                        disp(sprintf('Empty station at %i',K));
                        station_single(K)=station(K);
                        % station_single(K)= compute_bearings_station(station(K),goodFile{K},param,run_options.bearing_alg,0,run_options.kappa_Nsamples);
                        
                    end
                    
                    
                end
                %tsv_out_single=sprintf('%s_%s_singleton',keyword.stage,run_options.localization_alg);
                [fname_tsv_single,Nwrite]=write_tsv_singleton(station,data.goodName,[fname_out_tsv '_singleton']);
            catch
                fprintf('Crash! date_range(Idate,:) =%s singleton failure\n',date_range(Idate,:));
            end
            
        else
            disp('locations is empty, write an empty file');
            locations=[];locations_ctime=[];linking_index=data.linking_index;
            station=data.station;station_single=data.station_single;
        end
        eval(sprintf('!cp %s/%s %s',search_dir,fnames.name,fname_out));
        save(fname_out,'-append','locations','locations_ctime','linking_index','station','station_single');
        % eval(sprintf('!mv %s %s/',fname_out,search_dir));
        %tsv_out=sprintf('%s_%s',keyword.stage,run_options.localization_alg);
        [fname_tsv,Nwrite]=write_tsv(locations,data.goodName,fname_out_tsv);
    end

    function strip_industrial_from_locations
        %strip_industrial_from_locations;
        Igood=find(data.locations_ctime(:,Idasar)>0);
        for K=1:length(Igood)
            Iloc=Igood(K);
            [dtabs,Imatch]=min(abs(locations_tabs(Iloc,Idasar)-tabs));
            
            if debug.flag==1&&debug.dasar==Idasar
                if abs(debug.tabs-locations_tabs(Iloc,Idasar))<datenum(0,0,0,1,0,5)
                    fprintf('debug time hit in DASAR %i\n',Idasar);
                end
            end
            
            %%Save some time by skipping if no shipping metric available...
            %if Nf==0
            %    continue
            % end
            
            %%Search for matches in location variable in current DASAR
            Fmin=data.locations{Iloc}.Totalfmin(Idasar);
            Fmax=data.locations{Iloc}.Totalfmax(Idasar);
            Fmid=0.5*(Fmin+Fmax);
            
            %                Fmin=data.locations{Iloc}.feature(Idasar).robust_fmin;
            %                Fmax=data.locations{Iloc}.feature(Idasar).robust_fmax;
            %                Fmid=0.5*(Fmin+Fmax);
            %
            solidity=data.locations{Iloc}.feature(Idasar).Solidity;
            locations_fmin(Iloc,Idasar)=Fmin;
            locations_fmax(Iloc,Idasar)=Fmax;
            locations_solidity(Iloc,Idasar)=solidity;
            
            try
                Nf=max(adp_Nf{Ipeak}(Imatch+(-time_window:time_window)));
            catch
                try
                    Nf=adp_Nf{Ipeak}(Imatch);
                catch
                    Nf=adp_Nf{Ipeak}(end);
                end
            end
            
            if Nf>0
                
                try
                    %Nf=adp_Nf{Ipeak}(Imatch);
                    F=adp_F{Ipeak}(1:Nf,Imatch+(-time_window:time_window));
                    F=unique(F(:));
                    %F=F(2:end); %Remove zero
                catch
                    disp('Too close to start or end of file');
                    try
                        F=adp_F{Ipeak}(1:Nf,Imatch);
                    catch
                        F=adp_F{Ipeak}(1:Nf,end);
                    end
                end
                
                
                [dfmid,Ibest]=min(abs(Fmid-F));
                [dfmin]=min(abs(Fmin-F));
                [dfmax]=min(abs(Fmax-F));
                
                %Store stats for later...
                locations_dfmin(Iloc,Idasar)=dfmin;
                locations_dfmax(Iloc,Idasar)=dfmax;
                
            else
                locations_dfmin(Iloc,Idasar)=499;
                locations_dfmax(Iloc,Idasar)=499;
                dfmin=NaN;
                dfmax=NaN;
                dfmid=NaN;
            end
            
            if eval(rule)
                locations_ctime(Iloc,Idasar)=-1;  %Strip this out...
                
            end
            
            
            
            %%Search for matches using industrial frequencies across all
            %%DASARS
            %                 Nf=adp_Nf_combined{Ipeak}(Imatch);
            %                 if Nf>0
            %                     F=adp_F_combined{Ipeak}(1:Nf,Imatch);
            %                     [dfmid,Ibest]=min(abs(Fmid-F));
            %                     [dfmin]=min(abs(Fmin-F));
            %                     [dfmax]=min(abs(Fmax-F));
            %                 end
            %                 if eval(rule)
            %                     locations_ctime_combined(Iloc,Idasar)=-1;  %Strip this out...
            %
            %                 end
            
        end  %Igood
    end %strip_Bulk_Processing_Shipping_index

    function strip_industrial_from_stations
        
        
        for Iloc=1:length(station_tabs_min{Idasar})
            %Iloc=Igood(K);
            
            [~,Imatch]=min(abs(station_tabs_min{Idasar}(Iloc)-tabs));
            
            %%Search for matches in location variable in current DASAR
            Fmin=data.station(Idasar).Totalfmin(Iloc);
            Fmax=data.station(Idasar).Totalfmax(Iloc);
            Fmid=0.5*(Fmin+Fmax);
            solidity=data.station(Idasar).feature.Solidity(Iloc);
            
            try
                Nf=max(adp_Nf{Ipeak}(Imatch+(-time_window:time_window)));
            catch
                try
                    Nf=adp_Nf{Ipeak}(Imatch);
                catch
                    Nf=adp_Nf{Ipeak}(end);
                    
                end
            end
            
            if Nf>0
                try
                    %Nf=adp_Nf{Ipeak}(Imatch);
                    F=adp_F{Ipeak}(1:Nf,Imatch+(-time_window:time_window));
                    F=unique(F(:));
                    
                catch
                    disp('Too close to start or end of file');
                    try
                        F=adp_F{Ipeak}(1:Nf,Imatch);
                    catch
                        F=adp_F{Ipeak}(1:Nf,end);
                        
                    end
                end
                
                
                [dfmid,Ibest]=min(abs(Fmid-F));
                [dfmin]=min(abs(Fmin-F));
                [dfmax]=min(abs(Fmax-F));
            else
                dfmid=NaN;
                dfmin=NaN;
                dfmax=NaN;
                
            end
            
            if eval(rule)
                station_tabs_min{Idasar}(Iloc)=-1;  %Strip this out...
            end
            
            
        end  %Igood
    end %strip_industrial_from_stations

    function strip_debug_plotting
        %%strip_debug_plotting.m%%%
        close all
        %Remove failed locations (only one bearings)
        for JJJ=1:1 %local and combined DASAR peak detections
            
            %Review and cleanse detections...
            if JJJ==1
                Igood=find(sum(locations_ctime'>0)>1);
                Ibad=find(sum(locations_ctime'>0)<2);  %Rejections if enough rejections to fail localization
                fprintf('Pass: %i, Rejected: %i\n',length(Igood),length(Ibad));
                % Igood=find(sum(locations_ctime'<0)==0);  %Reject localization if *any* failures...
                %Ibad=find(sum(locations_ctime'<0)>0);
                
                
                %Igood=setdiff(1:size(locations_ctime,1),Ibad);
                titstr='single DASAR';
                
            else
                %Review and cleanse detections...
                Igood=find(sum(locations_ctime_combined'>0)>1);
                Ibad=find(sum(locations_ctime_combined'>0)<2);
                titstr='combined';
            end
            disp(titstr);
            %plot original and cleaned results
            time_inc=24;
            run_options.center_dist_limit=50000;
            
            auto.locations{1}=data.locations(Ibad);
            auto.locations_ctime{1}=data.locations_ctime(Ibad,:);
            fprintf('Plotting rejects %s, choice %s\n',titstr, chc);
            
            try
                plot_movie_all_auto(DASAR_coords,[],auto,[],[],I,run_options.center_dist_limit,time_inc,goodFile);
            end
            
            title(sprintf('%i Rejected localizations at Site %i, Rule: %s, %s, File: %s',length(Ibad),I,rule,titstr,fnames.name),'interp','none','fontsize',8);
            orient landscape
            Idot=findstr(fnames.name,'.')-1;
            print('-djpeg',sprintf('%s_%s_rejected.jpg',titstr,fnames.name(1:Idot)));
            
            %Statistics on the frequency differences between detections and shipping metric...
            tmp=locations_dfmin(Ibad,:);
            tmp=tmp(:);
            Ipass= (tmp<=500);
            stats.reject.dfmin=tmp(Ipass);
            
            tmp=locations_dfmax(Ibad,:);
            tmp=tmp(:);
            Ipass= (tmp<=500);
            stats.reject.dfmax=tmp(Ipass);
            
            tmp=locations_fmin(Ibad,:);
            tmp=tmp(:);
            Ipass= (tmp<=500);
            stats.reject.fmin=tmp(Ipass);
            
            tmp=locations_fmax(Ibad,:);
            tmp=tmp(:);
            Ipass= (tmp<=500);
            stats.reject.fmax=tmp(Ipass);
            
            
            tmp=locations_solidity(Ibad,:);
            tmp=tmp(:);
            Ipass= (tmp<=500);
            stats.reject.solidity=tmp(Ipass);
            
            %     figure
            %     subplot(3,1,1);
            %     hist(stats.reject.dfmin,0:1:100);grid on
            %     xlabel('dfmin for rejected detections');
            %     subplot(3,1,2);
            %     hist(stats.reject.dfmax,0:1:100);grid on
            %     xlabel('dfmax for rejected detections');
            %     subplot(3,1,3);
            %     hist(stats.reject.fmin,0:1:100);grid on
            %     xlabel('duration for rejected detections');
            
            
            figure
            auto.locations{1}=data.locations(Igood);
            auto.locations_ctime{1}=data.locations_ctime(Igood,:);
            fprintf('Plotting accepted %s, choice %s\n',titstr, chc);
            
            try
                Ifalse_org=plot_movie_all_auto(DASAR_coords,[],auto,[],[],I,run_options.center_dist_limit,time_inc,goodFile);
            end
            title(sprintf('%i Retained localizations at Site %i, Rule: %s, %s File: %s',length(Igood),I,rule,titstr,fnames.name),'interp','none','fontsize',8);
            print('-djpeg',sprintf('%s_%s_retained',titstr,fnames.name(1:Idot)));
            
            %Statistics on the frequency differences between detections and shipping metric...
            
            Ifalse=Igood(Ifalse_org);
            tmp=locations_dfmin(Ifalse,:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.false.dfmin=tmp(Ipass);
            
            tmp=locations_dfmax(Ifalse,:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.false.dfmax=tmp(Ipass);
            
            tmp=locations_fmin(Ifalse,:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.false.fmin=tmp(Ipass);
            
            tmp=locations_fmax(Ifalse,:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.false.fmax=tmp(Ipass);
            
            tmp=locations_solidity(Ifalse,:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.false.solidity=tmp(Ipass);
            
            %     figure
            %     subplot(3,1,1);
            %     hist(stats.false.dfmin,0:1:500);grid on
            %     xlabel('dfmin for accepted detections');
            %     subplot(3,1,2);
            %     hist(stats.false.dfmax,0:1:500);grid on
            %     xlabel('dfmax for accepted detections');
            %     subplot(3,1,3);
            %     hist(stats.false.fmin,0:1:500);grid on
            %     xlabel('duration for accepted detections');
            %
            
            Itrue=setdiff(1:length(Igood),Ifalse_org);
            tmp=locations_dfmin(Igood(Itrue),:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.true.dfmin=tmp(Ipass);
            
            tmp=locations_dfmax(Igood(Itrue),:);
            tmp=tmp(:);
            Ipass= tmp<500;
            stats.true.dfmax=tmp(Ipass);
            
            tmp=locations_fmin(Igood(Itrue),:);
            tmp=tmp(:);
            Ipass=find(tmp<500);
            stats.true.fmin=tmp(Ipass);
            
            tmp=locations_fmax(Igood(Itrue),:);
            tmp=tmp(:);
            Ipass=find(tmp<500);
            stats.true.fmax=tmp(Ipass);
            
            tmp=locations_solidity(Igood(Itrue),:);
            tmp=tmp(:);
            Ipass=find(tmp<500);
            stats.true.solidity=tmp(Ipass);
            
            
            %     figure
            %     subplot(3,1,1);
            %     hist(stats.true.dfmin,0:1:500);grid on
            %     xlabel('dfmin for accepted detections');
            %     subplot(3,1,2);
            %     hist(stats.true.dfmax,0:1:500);grid on
            %     xlabel('dfmax for accepted detections');
            %     subplot(3,1,3);
            %     hist(stats.true.fmin,0:1:500);grid on
            %     xlabel('fmin for accepted detections');
            %
            %
            %     figure
            %     subplot(3,1,1);
            %     plot(stats.true.dfmin,stats.true.dfmax,'kx');grid on
            %     xlabel('dfmin accepted and true');ylabel('dfmax accepted but false');
            %     subplot(3,1,2)
            %     plot(stats.reject.dfmin,stats.reject.dfmax,'ro');grid on
            %     xlabel('dfmin rejected');ylabel('dfmax rejected');
            %     subplot(3,1,3);
            %     plot(stats.false.dfmin,stats.false.dfmax,'gx');grid on
            %     xlabel('dfmin accepted and false');ylabel('dfmax accepted but false');
            
            [N,printname,Ibin]=hist2D([stats.true.dfmin stats.true.dfmax]',1:5:100,1:5:100,{'dfmin','dfmax'});
            title('true');
            
            [N,printname,Ibin]=hist2D([stats.reject.dfmin stats.reject.dfmax]',1:5:100,1:5:100,{'dfmin','dfmax'});
            title('reject');
            
            [N,printname,Ibin]=hist2D([stats.false.dfmin stats.false.dfmax]',1:5:100,1:5:100,{'dfmin','dfmax'});
            title('false');
            
            [N,printname,Ibin]=hist2D([stats.true.dfmin stats.true.fmin]',1:5:100,0:5:500,{'dfmin','fmin'});
            title('true');
            
            [N,printname,Ibin]=hist2D([stats.reject.dfmin stats.reject.fmin]',1:5:100,1:5:500,{'dfmin','fmin'});
            title('reject');
            
            [N,printname,Ibin]=hist2D([stats.false.dfmin stats.false.fmin]',1:5:100,1:5:500,{'dfmin','fmin'});
            title('false');
            
            
            [N,printname,Ibin]=hist2D([stats.true.fmin stats.true.fmax-stats.true.fmin]',1:2:500,1:2:50,{'fmin','fmax-fmin'});
            title('true');
            
            [N,printname,Ibin]=hist2D([stats.reject.fmin stats.reject.fmax-stats.reject.fmin]',1:2:500,1:2:50,{'fmin','fmax-fmin'});
            title('reject');
            
            [N,printname,Ibin]=hist2D([stats.false.fmin stats.false.fmax-stats.false.fmin]',1:2:500,1:2:50,{'fmin','fmax-fmin'});
            title('false');
            
            
            
            [N,printname,Ibin]=hist2D([stats.true.fmin stats.true.solidity]',1:5:200,0:0.01:1,{'fmin','solidity'});
            title('true');
            
            [N,printname,Ibin]=hist2D([stats.reject.fmin stats.reject.solidity]',1:5:200,0:0.01:1,{'fmin','solidity'});
            title('reject');
            
            [N,printname,Ibin]=hist2D([stats.false.fmin stats.false.solidity]',1:5:200,0:0.01:1,{'fmin','solidity'});
            title('false');
            
        end  %strip debug plotting
    end

end  %strip_Bulk_Processing
