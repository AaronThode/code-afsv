%%[Ifinal]=plot_movie_all_auto(DASAR_coords,manual,auto,miss_stats,auto_stats,Isite, ...
%%   manual_limits,time_inc,goodFile)
%
%%% Plot locations of all automated detections, with reference to manual results
%%
%% Inputs:
%%    DASAR_coords:
%%    manual: structure with manual detection information loaded from load_manual_results.m
%%    auto:   structure with auto detection information loaded from load_automated_results.m
%%    miss_stats: output of 'stats' in evalute_linking.  Contains fields
%%       with lengths equal to number of manual localizations.
%%          .match_flag:  if 1, an automated detection matches this manual
%%          result.
%%    auto_stats:  output of 'auto_stats' in evaluate_linking.  Contains
%%      fields with lengths equal to number of automated locations.
%%          .match_flag:  if 1, a manual detection matches this automated
%%          result.
%%    Isite: scalar number of Site location
%%    manual_limits:  plotting limits, lat/long
%%    time_inc: time increment in hours
%%    param:  parameter structure:
%%%             .Nframes:
%%              
%%    goodFile: cell array of strings showing absolute files names of associated DASAR data.  Used for
%%      spectrograms

%% OUtputs:
%%      Ifinal:  indicies of locations that were selected...
            

function Ifinal=plot_movie_all_auto(DASAR_coords,manual,auto,miss_stats, ...
    auto_stats,Isite,manual_limits,time_inc,goodFile)

%strr='ABCDEFGHIJKLMNOP';
Ifinal=[];
%Construct DASAR labels before anything else happens
strr=[];
for II=1:length(goodFile{1})
    Islash=max(findstr(goodFile{1}{II},'/'));
    %disp(goodFile{1}{II}(Islash+6-1));
    strr=[strr goodFile{1}{II}(Islash+6-1)];
end
review_status='look';  %%Can be 'look' or 'write'
auto_corrected=[];  %This is a vector of autocorrected data..


if isinf(manual_limits)
    disp('Changing infinite limits to 50 km')
    manual_limits=50000;
end

loc_index=[];
if Isite==1
    xlimm=[-30 30];
    ylimm=[-30 30];
else
    xlimm=1.25*manual_limits*[-1 1]/1000;
    ylimm=[-30 30];
    ylimm=[-30 20];
    %xlimm=[-30 30];
    
end

%%Find midnight of the day in question, based on time of first detection..
It=find(auto.locations_ctime{1}(1,:)>0);
ctime_min=auto.locations_ctime{1}(1,It(1));
tabs_min=(datenum(1970,1,1,0,0,ctime_min));
tvec=datevec(tabs_min);
tabs_midnt=datenum(tvec(1),tvec(2),tvec(3));

%Set ctime1 to midnight before the first detection...
sec_elaps=datevec(tabs_min-tabs_midnt);
ctime1=ctime_min-sec_elaps(6)-60*sec_elaps(5)-3600*sec_elaps(4);

if sec_elaps(3)>0
    sec_elaps=datevec(-tabs_min+tabs_midnt);
    ctime1=ctime_min+sec_elaps(6)+60*sec_elaps(5)+3600*sec_elaps(4);
end

ctime2=ctime1+24*60*60;  %Ctime2 is the following midnight...

% if ~isfield(param,'Nframes')
%     Nframes=25;
%     Nframes=7;
% else
%     Nframes=param.Nframes;
% end

if ~exist('time_inc','var')
    time_inc=4;
end

ctime_range=ctime1:(time_inc*3600):ctime2;

%%Create a grid of times for location events, even if localization fails...


Nloc_fail=0;
for Iauto=1:size(auto.locations_ctime{1},1);
    if isfield(auto.locations{1}{Iauto}.position,'ctime')
        auto_ctime(Iauto)=auto.locations{1}{Iauto}.position.ctime;
        
    else
        %disp('Warning!! Unsuccesful localizations present in data');
        Nloc_fail=Nloc_fail+1;
    end
end
fprintf('Warning!! %i unsuccesful localizations present in data..',Nloc_fail);

Icount=0;
for Iindex=1:(length(ctime_range)-1)  %Fore every time block
    
    figure(1);set(gcf,'pos',[352         296        1197         747]);
    clf;
    
    %%%Manual plotting%%%%%%%%%
    if ~isempty(manual)&&~isempty(manual.localized)
        Icall=find(manual.localized{1}.ctev>=ctime_range(Iindex)&manual.localized{1}.ctev<=ctime_range(Iindex+1));
        VM=[manual.localized{1}.utmx(Icall) manual.localized{1}.utmy(Icall)];
        Ipass=find(miss_stats.match_flag(Icall)>0);
        Imiss=find(miss_stats.match_flag(Icall)==0);
        
        % keyboard;
        
        %Plot manual results
        subplot(2,2,1)
        plot_location(Isite,DASAR_coords,[],[],VM(Ipass,:)) ;
        %plot_letter_label('a)');
        text(xlimm(1)+5,ylimm(2)-5,'a)','fontweight','bold','fontsize',14);
        set(gca,'xtick',linspace(xlimm(1),xlimm(2),20));
        set(gca,'xticklabel',[]);
        xlabel('');
        xlim(xlimm);ylim(ylimm);
        hold on
        %title(sprintf('%i Manual calls matched %s - %s',Nmanual(1),ctime2str(ctime_range(Iindex)),ctime2str(ctime_range(Iindex+1))));
        text(-3,15,sprintf('%s - %s',ctime2str(ctime_range(Iindex)),ctime2str(ctime_range(Iindex+1))));
        text(-3,13,(sprintf('%i manual matches out of %i',length(Ipass),length(miss_stats.match_flag(Icall)))),'fontweight','bold','fontsize',14);
        
        %title(sprintf('%i manual calls matched out of %i',length(Ipass),length(miss_stats.match_flag(Icall))));
        
        
        subplot(2,2,4);
        Dn=plot_location(Isite,DASAR_coords,[],[],VM(Imiss,:)) ;
        text(-3,13,(sprintf('%i manual calls unmatched',length(Imiss))),'fontweight','bold','fontsize',14);
        %plot_letter_label('d)');
        text(xlimm(1)+5,ylimm(2)-5,'d)','fontweight','bold','fontsize',14);
        
        set(gca,'xtick',-50:5:50)
        ylabel('');set(gca,'yticklabel',[]);
        %set(gca,'xticklabel',[]);
        %title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));
        xlim(xlimm);ylim(ylimm);
        hold on;
    end
    
    %%%%%%March through automated calls%%%%%%%%
    %Icall_auto: top level index of automated calls that fall within time band
    Icall_auto=find(auto_ctime>=ctime_range(Iindex)&auto_ctime<=ctime_range(Iindex+1));
    VM_auto=zeros(length(Icall_auto),2);
    for I=1:length(Icall_auto)
        try
            VM_auto(I,:)=[auto.locations{1}{Icall_auto(I)}.position.location];
        catch
            disp('auto localization failure, revise master_setup to permit succesful localizations only');
            
        end
    end
    
    %Ipass_auto and Ifalse_auto have length(Icall_auto)--only refer to
    %   current time periods
    
    if isempty(manual)||isempty(manual.localized)
        Ifalse_auto=1:length(Icall_auto);
        Dn=plot_location(Isite,DASAR_coords,[],[],VM_auto) ;
        
        set(gca,'xtick',-50:5:50);
        
        text(-3,13,((sprintf('%i automated calls',length(Ifalse_auto)))),'fontweight','bold','fontsize',14);
         text(-3,15,sprintf('%s - %s',ctime2str(ctime_range(Iindex)),ctime2str(ctime_range(Iindex+1))));
       
        %title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));
        xlim(xlimm);ylim(ylimm);
        hold on;
        
        if Isite==1
            line(-6.8+manual_limits*[1 1],[-1 20]);
            line(-6.8-manual_limits*[1 1],[-1 20]);
            line(10.7368+manual_limits*[1 1],[-15 -5]);
            line(10.7368-manual_limits*[1 1],[-15 -5]);
            
        else
            
            line(manual_limits*[1 1]/1000,ylimm);
            line(-manual_limits*[1 1]/1000,ylimm);
            
            set(gca,'xminortick','on','yminorgrid','on','yminortick','on');
            poss=get(gca,'pos');
            
%             if I<3,poss(2)=.55;end
%             poss(3)=.4;
%             poss(4)=0.42;
%             set(gca,'pos',poss);
        end
    else  %if manual and automated data present
        Ipass_auto=find(auto_stats.match_flag(Icall_auto)>0);
        Ifalse_auto=find(auto_stats.match_flag(Icall_auto)==0);
        
        
        subplot(2,2,2)
        Dn=plot_location(Isite,DASAR_coords,[],[],VM_auto(Ipass_auto,:)) ;
        %plot_letter_label('b)');
        text(xlimm(1)+5,ylimm(2)-5,'b)','fontweight','bold','fontsize',14);
        
        set(gca,'xtick',-50:5:50)
        xlabel('');ylabel('');
        set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
        
        
        %title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));
        text(-3,13,sprintf('%i auto calls matched',length(Ipass_auto)),'fontweight','bold','fontsize',14);
        xlim(xlimm);ylim(ylimm);
        hold on;
        
        subplot(2,2,3)
        Dn=plot_location(Isite,DASAR_coords,[],[],VM_auto(Ifalse_auto,:)) ;
        %plot_letter_label('c)');
        text(xlimm(1)+5,ylimm(2)-5,'c)','fontweight','bold','fontsize',14);
        
        set(gca,'xtick',-50:5:50);
        
        text(-3,13,((sprintf('%i automated calls unmatched',length(Ifalse_auto)))),'fontweight','bold','fontsize',14);
        %title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));
        xlim(xlimm);ylim(ylimm);
        hold on;
        disp(sprintf('Matched manual calls: %i, missed manual calls: %i, matched auto calls: %i, extra auto calls: %i', length(Ipass),length(Imiss),length(Ipass_auto),length(Ifalse_auto)));
        if Isite==1
            for I=1:4
                subplot(2,2,I)
                line(-6.8+manual_limits*[1 1],[-1 20]);
                line(-6.8-manual_limits*[1 1],[-1 20]);
                line(10.7368+manual_limits*[1 1],[-15 -5]);
                line(10.7368-manual_limits*[1 1],[-15 -5]);
            end
        else
            for I=1:4
                subplot(2,2,I)
                line(manual_limits*[1 1]/1000,ylimm);
                line(-manual_limits*[1 1]/1000,ylimm);
                
                set(gca,'xminortick','on','yminorgrid','on','yminortick','on');
                poss=get(gca,'pos');
                
                if I<3,poss(2)=.55;end
                poss(3)=.4;
                poss(4)=0.42;
                set(gca,'pos',poss);
            end
        end
    end
    
   
    
    
    Icount=Icount+1;
    FF(Icount)=getframe(gcf);
    
    %%Flagging a particular detection
    Iyes=1;
    while Iyes<3
        
        %If no manual results are present, just review the automated data.
        if isempty(manual)||isempty(manual.localized)
            Iyes=menu('Select a region?','Yes','No, gracias');
            if Iyes==2
                Iyes=3;  %No reviews
            end
            
        else
            Iyes=menu('Select a region?','False Alarms','Misses','No, gracias');
        end
        switch Iyes
            case 1
                disp('Click on two opposite corners');
                tmp=ginput(2);
                
                xg=mean(DASAR_coords(:,1));
                yg=mean(DASAR_coords(:,2));
                tmp=tmp*1000+[1;1]*[xg yg];
                
                xl=min(tmp(:,1));
                xu=max(tmp(:,1));
                yu=max(tmp(:,2));
                yl=min(tmp(:,2));
                %Iextra=setdiff(Icall_auto,miss_stats.Nlocs_all_match);
                %Iplot is referenced to Icall_auto--the local time segment only
                Iplot=[];
                for IJK=1:length(Ifalse_auto)
                    if (VM_auto(Ifalse_auto(IJK),1)>=xl&&VM_auto(Ifalse_auto(IJK),1)<=xu&&VM_auto(Ifalse_auto(IJK),2)>=yl&&VM_auto(Ifalse_auto(IJK),2)<=yu)
                        Iplot=[Iplot Ifalse_auto(IJK)];
                    end
                end
                h=plot((VM_auto(Iplot,1)-xg)/1000,(VM_auto(Iplot,2)-yg)/1000,'s','markersize',5);hold on
                
                
                %%Histogram features associated with these detections(outside subfunction)
                %Icall_auto: indicies within timespan. Iplot: indicies of Icall_auto within location area.
                Ifinal=Icall_auto(Iplot);  %Ifinal indexes original auto object.
                fprintf('Total number of locations in this window: %i, Selected: %i; Not selected: %i\n',length(Ifalse_auto),length(Ifinal),length(Ifalse_auto)-length(Ifinal));
                [~,Ibinn]=histogram_features(auto,auto_ctime, Ifinal, goodFile{1});
                
                %Ibinn are indicies of original Ifinal that become output Ifinal: Ifinal_out=Ifinal_in(Ibinn);
                if ~isempty(Ibinn)
                    Iplot=Iplot(Ibinn);
                    figure(1)
                    h=plot((VM_auto(Iplot,1)-xg)/1000,(VM_auto(Iplot,2)-yg)/1000,'^','markersize',5,'color',[1 0 1]);hold on
                    Ifinal=Icall_auto(Iplot);
                end
                
                %%%Display spectrograms associated with these detections
                Icheck_false=menu('Review spectrograms?','Yes','No');
                if Icheck_false==1
                    review_false_alarms;
                end
                
                %loc_index=[loc_index Iplot];
                
            case 2  %Review misses...
                disp('Click on two opposite corners');
                tmp=ginput(2);
                
                xg=mean(DASAR_coords(:,1));
                yg=mean(DASAR_coords(:,2));
                tmp=tmp*1000+[1;1]*[xg yg];
                
                xl=min(tmp(:,1));
                xu=max(tmp(:,1));
                yu=max(tmp(:,2));
                yl=min(tmp(:,2));
                %Iextra=setdiff(Icall_auto,miss_stats.Nlocs_all_match);
                %Iplot is referenced to Icall_auto--the local time segment only
                Iplot=[];
                for I=1:length(Imiss)
                    if (VM(Imiss(I),1)>=xl&VM(Imiss(I),1)<=xu&VM(Imiss(I),2)>=yl&&VM(Imiss(I),2)<=yu)
                        Iplot=[Iplot Imiss(I)];
                    end
                end
                h=plot((VM(Iplot,1)-xg)/1000,(VM(Iplot,2)-yg)/1000,'s','markersize',5);hold on
                
                %%%Display spectrograms associated with these detections
                review_misses;
                
                
                %loc_index=[loc_index Iplot];
                
            otherwise
        end
    end %while Iyes
   % clf
    
end  %Iindex
%movie2avi(FF,movie_name,'fps',5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%function review_false_alarms.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%display spectrograms of false alarms and classify reasons for mismatch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function review_false_alarms
        Iyes2=menu('Display associated spectrograms?','Yes','No');
        if Iyes2==0
            return
        end
        Iplot_ref=1;
        Nplots=length(Iplot);
        %Nplots=input((sprintf('Enter number of samples(%i), enter -1 to just list sorted times:',Nplots)));
        
        prompt={'Enter number of samples(%i), enter -1 to just list sorted times:'};
        name='Sample selection';
        numlines=1;
        defaultanswer={int2str(Nplots)};
        
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        Nplots=str2num(answer{1});
        
        if isempty(Nplots)
            Nplots=length(Iplot);
        end
        
        if Nplots<0  %Print out time distribution
            %datestr(datenum(1970,1,1,0,0,auto_ctime(Iplot_ref)))
            datetemp=datenum(1970,1,1,0,0,auto_ctime(Iplot));
            hist(datetemp,10);datetick('x',6);
            keyboard
            return
        end
        
        %Randomize detections if desired
        if length(Nplots)==1&Nplots<length(Iplot)
            %plot_index=Iplot(ceil(length(Iplot)*rand(1,Nplots)));  %Random selection
            plot_index=Iplot(1:Nplots);
        elseif length(Nplots)>1
            Nplots=Nplots(Nplots<=length(Iplot));
            plot_index=Iplot(Nplots);
            
        else
            plot_index=Iplot;
        end
        for II=plot_index  %II is a member of Iplot
            Iplot_ref=Icall_auto(II);
            %%Create a new map showing call information
            figure;
            set(gcf,'pos',[ 23          66         815        1011]);
            
            A=auto.locations{1}{Iplot_ref}.position.major;
            B=auto.locations{1}{Iplot_ref}.position.minor;
            ANG=auto.locations{1}{Iplot_ref}.position.ellipse_ang;
            Dn=plot_location(Isite,DASAR_coords,[],[],VM_auto(II,:),A,B,ANG) ;
            xlim(xlimm);ylim(ylimm);hold on
            Istations=(auto.locations{1}{Iplot_ref}.station_indicies>0);
            plot(Dn(Istations,1),Dn(Istations,2),'^b','markerfacecolor',[0 0 0],'markersize',10);hold on
            for JJJ=1:size(Dn,1)
                text(Dn(JJJ,1)-2,Dn(JJJ,2)-2,strr(JJJ),'fontweight','bold','fontsize',10)
            end
            %plot_overlay_location(DASAR_coords(Istations,:),'^b');
            
            if Isite==5&&~isempty(findstr(goodFile{1}{1},'T2010'))
                VA_coords =[4.174860699660919e+05 7.817274204098196e+06];
                rng=sqrt(sum((VA_coords-VM_auto(II,:)).^2));
                title(sprintf('Distance from vertical array is %6.2f km',rng/1000));
            end
            
            %%Display spectrogram of links in question
            param.Nfft=256;param.ovlap=.75; param.Fs=1000;
            [ysound,tabss]=display_automated_crosslink_data(auto.locations{1},Iplot_ref,param,goodFile{1});
            %  set(gcf,'pos',[0.4479    0.0526    0.4495    0.8582]);
            %text(1,400,sprintf('Localization num: %i',Iplot_ref));
            filter_names={'Contour_global_bandwidth','Contour_fmax','Contour_fmin'};
            
            for K=1:length(filter_names)
                fprintf('%s: %s\n',filter_names{K},mat2str(auto.locations{1}{Iplot_ref}.(filter_names{K})));
                
            end
            
            ori=[auto.locations{1}{Iplot_ref}.feature.Orientation];
            Ipss=find(~isnan(ori));
            for K=Ipss
                fprintf('Orientation DASAR %i: %6.2f\n', K,ori(K));
            end
            
            Ifail=1;
            while Ifail==1
                choices={'Pass, Manual review missed','Pass, alternate linkage','Fail, weak or blob','Fail, seal','Fail, seal fragment','Fail, walrus','Fail,airgun','Fail, ship','Fail, bad call linkage', ...
                    'Uncertain (further review)','Play sound','Review closest manual','Print spectrogram','Results NOT SAVED'};
                Nchc=length(choices)-4;
                if strcmp(review_status,'write')
                    choices{end}='Results SAVED';
                end
                Ichc=menu('What is your decision?',choices);
                if Ichc<=Nchc
                    Ifail=0;
                    if ~isempty(findstr(review_status,'write'))&&(~isempty(manual)&&~isempty(manual.localized))
                        name='Optional comment';
                        numlines=1;
                        defaultanswer={'none'};
                        
                        try
                            answer=inputdlg(prompt,name,numlines,defaultanswer);
                            comment=answer{1};
                        catch
                            
                            comment=[];
                        end
                        auto_corrected=update_reviewed_results_false(auto_corrected,auto.locations{1}{Iplot_ref},Iplot_ref,choices{Ichc},comment,Isite);
                        
                    end
                    
                elseif Ichc==Nchc+1;
                    disp('Play sound');
                    Icall=menu('which DASAR?',goodFile{1});
                    soundsc(ysound{Icall},param.Fs);
                elseif Ichc==Nchc+2&&(~isempty(manual)&&~isempty(manual.localized))
                    
                    %Display closest manual result to each link in
                    %question...
                    
                    param.energy.bufferTime=2;
                    myctime=auto.locations{1}{Iplot_ref}.ctime_min;
                    myctime2=auto.locations_ctime{1}(Iplot_ref,:);
                    dt=zeros(1,length(myctime));
                    Icall=dt;
                    Nstations=length(myctime);
                    for K=1:Nstations
                        if myctime(K)>0
                            [dt(K),Icall(K)]=min(abs(myctime(K)-manual.individual{1}.ctime(:,K)));
                            disp(sprintf('Station %s manual best match is %i, dt is %6.2f sec',strr(K),Icall(K),dt(K)));
                            ctimes_out=display_manual_location_data(manual.localized{1},manual.individual{1},Icall(K),param,goodFile{1});
                            subplot(4,2,Nstations-K+1);
                            text(1,400,sprintf('dt is %6.2f sec',dt(K)),'fontweight','bold','fontsize',14);
                            
                            pause;
                            close(max(get(0,'child')));
                            disp(sprintf('\n\n'));
                        end
                    end
                elseif Ichc==Nchc+3
                    figure(gcf)
                    orient landscape
                    tmp=strr(Istations);
                    print('-djpeg',sprintf('D%s%s_allDASARs_specgram',tmp(end),datestr(tabss(end),30)));
                    figure(gcf-1)
                    orient tall
                    print('-djpeg',sprintf('D%s%s_allDASARs_map',tmp(end),datestr(tabss(end),30)));
                    
                elseif Ichc==length(choices)
                    if strcmp(review_status,'look')
                        review_status='write';
                    else
                        review_status='look';
                    end
                else
                    Ifail=1;
                end
                
                
            end
            
            
            Ifig=get(0,'child');
            close(setdiff(Ifig,1));
            
            
        end
    end

    function review_misses
        Iyes2=menu('Display associated spectrograms?','Yes','No');
        if Iyes2==0
            return
        end
        Iplot_ref=1;
        Nplots=length(Iplot);
        %Nplots=input((sprintf('Enter number of samples(%i), enter -1 to just list sorted times:',Nplots)));
        
        prompt={'Enter number of samples(%i), enter -1 to just list sorted times:'};
        name='Sample selection';
        numlines=1;
        defaultanswer={int2str(Nplots)};
        
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        Nplots=str2num(answer{1});
        
        if isempty(Nplots)
            Nplots=length(Iplot);
        end
        
        if Nplots<0
            datestr(datenum(1970,1,1,0,0,manual.localized{1}.ctev(Icall)))
            return
        end
        
        %Randomize detections if desired
        if Nplots<length(Iplot)
            %plot_index=Iplot(ceil(length(Iplot)*rand(1,Nplots)));  %Random selection
            plot_index=Iplot(1:Nplots);
        else
            plot_index=Iplot;
        end
        for II=plot_index  %II is a member of Iplot
            Iplot_ref=Icall(II);
            %%Create a new map showing call information
            figure;
            set(gcf,'pos',[ 23          66         815        1011]);
            Dn=plot_location(Isite,DASAR_coords,[],[],VM(II,:)) ;
            xlim(xlimm/2);ylim(ylimm);hold on
            %Istations=(auto.locations{1}{Iplot_ref}.station_indicies>0);
            Istations=find(manual.individual{1}.ctime(Iplot_ref,:)>0);
            plot(Dn(Istations,1),Dn(Istations,2),'^b','markerfacecolor',[0 0 0],'markersize',10);hold on
            for JJJ=1:size(Dn,1)
                text(Dn(JJJ,1)-2,Dn(JJJ,2)-2,strr(JJJ),'fontweight','bold','fontsize',10)
            end
            %plot_overlay_location(DASAR_coords(Istations,:),'^b');
            
            %%Display spectrogram of links in question
            param.Nfft=256;param.ovlap=.75; param.Fs=1000;
            param.energy.bufferTime=2;
            
            display_manual_location_data(manual.localized{1},manual.individual{1},Iplot_ref,param,goodFile{1});
            
            %ysound=display_automated_crosslink_data(auto.locations{1},Iplot_ref,param,goodFile{1});
            set(gcf,'pos',[0.4479    0.0526    0.4495    0.8582]);
            %text(1,400,sprintf('Localization num: %i',Iplot_ref));
            %             filter_names={'Contour_global_bandwidth','Contour_fmax','Contour_fmin'};
            %
            %             for K=1:length(filter_names)
            %                 disp(sprintf('%s: %s',filter_names{K},mat2str(auto.locations{1}{Iplot_ref}.(filter_names{K}))));
            %
            %             end
            
            
            
            Ifail=1;
            while Ifail==1
                choices={'Mistaken manual','Pass, alternate linkage','Fail, not marked energy','Fail, not marked interval','Fail, not marked image','Fail, not marked net','Fail, not linked','Play sound','Review closest automated','Results NOT SAVED'};
                Nchc=length(choices)-3;
                if strcmp(review_status,'write')
                    choices{end}='Results SAVED';
                end
                Ichc=menu('What is your decision?',choices);
                if Ichc<=Nchc
                    Ifail=0;
                    if ~isempty(findstr(review_status,'write'))&&(~isempty(manual)&&~isempty(manual.localized))
                        name='Optional comment';
                        numlines=1;
                        defaultanswer={'none'};
                        
                        try
                            answer=inputdlg(prompt,name,numlines,defaultanswer);
                            comment=answer{1};
                        catch
                            
                            comment=[];
                        end
                        auto_corrected=update_reviewed_results_missed(auto_corrected,manual.individual{1},Iplot_ref,choices{Ichc},comment,Isite);
                        
                    end
                    
                elseif Ichc==Nchc+1;
                    disp('Play sound not enabled');
                    %Icall=menu('which DASAR?',goodFile{1});
                    %soundsc(ysound{Icall},param.Fs);
                elseif Ichc==Nchc+2
                    
                    %Display closest automated result to each link in
                    %question...
                    
                    %myctime=auto.locations{1}{Iplot_ref}.ctime_min;
                    myctime=manual.individual{1}.ctime(Iplot_ref,:);
                    dt=zeros(1,length(myctime));
                    II=dt;
                    Nstations=length(myctime);
                    for K=1:Nstations
                        if myctime(K)>0
                            [dt(K),II(K)]=min(abs(myctime(K)-auto.locations_ctime{1}(:,K)));
                            
                            
                            disp(sprintf('Station %s manual best match is %i, dt is %6.2f sec',strr(K),II(K),dt(K)));
                            display_automated_crosslink_data(auto.locations{1},II(K),param,goodFile{1});
                            %ctimes_out=display_manual_location_data(manual.localized{1},manual.individual{1},Icall(K),param,goodFile{1});
                            subplot(4,2,Nstations-K+1);
                            text(1,200,sprintf('linking: dt is %6.2f sec',dt(K)),'fontweight','bold','fontsize',10,'color','w');
                            
                            [dt(K),II(K)]=min(abs(myctime(K)-auto.raw_stations{1}(K).raw_detections.ctime));
                            text(1,400,sprintf('energy: dt is %6.2f sec',dt(K)),'fontweight','bold','fontsize',10,'color','w');
                            
                            [dt(K),II(K)]=min(abs(myctime(K)-auto.raw_stations{1}(K).interval_detections.ctime));
                            text(1,350,sprintf('interval: dt is %6.2f sec',dt(K)),'fontweight','bold','fontsize',10,'color','w');
                            
                            [dt(K),II(K)]=min(abs(myctime(K)-auto.raw_stations{1}(K).morph_detections.ctime));
                            text(1,300,sprintf('image: dt is %6.2f sec',dt(K)),'fontweight','bold','fontsize',10);
                            
                            [dt(K),II(K)]=min(abs(myctime(K)-auto.stations{1}(K).ctime_min));
                            text(1,250,sprintf('neural net: dt is %6.2f sec',dt(K)),'fontweight','bold','fontsize',10,'color','w');
                            
                            pause;
                            close(max(get(0,'child')));
                            disp(sprintf('\n\n'));
                        end
                    end
                elseif Ichc==length(choices)
                    if strcmp(review_status,'look')
                        review_status='write';
                    else
                        review_status='look';
                    end
                end
                
            end
            
            
            Ifig=get(0,'child');
            close(setdiff(Ifig,1));
            
            
        end
    end

%%plot_location(Isite,DASAR_coords,bearings,VM,A,B,ANG);
% DASAR_coords UTM, bearings in degrees, VM in UTM
    function  DASAR_coordsn=plot_location(Isite,DASAR_coords,bearings,Igood,VM,A,B,ANG,linel)
        
        %LL=3;
        
        
        if nargin==4
            VM=[];
            A=[];
            B=[];
            ANG=[];
        elseif nargin==5
            A=[];
            B=[];
            ANG=[];
        end
        
        
        %%Should we plot other features, like a vertical array?
        VA_coords=[];
        if Isite==5&&~isempty(findstr(goodFile{1}{1},'T2010'))
            VA_coords =[4.174860699660919e+05 7.817274204098196e+06]/1000;
        end
        
        %Convert to km
        VM=VM/1000;
        DASAR_coords=DASAR_coords/1000;
        A=A/1000;
        B=B/1000;
        
        if ~exist('linel')
            linel=35; %length of bearing lines in km
        end
        %subplot(3,1,LL);
        xg=mean(DASAR_coords(:,1));
        yg=mean(DASAR_coords(:,2));
        DASAR_coordsn(:,1)=DASAR_coords(:,1)-xg;
        DASAR_coordsn(:,2)=DASAR_coords(:,2)-yg;
        
        plot(DASAR_coordsn(:,1),DASAR_coordsn(:,2),'r^','markersize',8,'markerfacecolor',[1 0 0]);hold on
        if ~isempty(VA_coords)
            plot(VA_coords(1)-xg,VA_coords(2)-yg,'gd','markersize',8,'markerfacecolor',[0 1 1]);hold on
            
        end
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Easting (km)');
        ylabel('Northing (km)');
        grid on;
        
        %Convert bearings from nautical to mathematical frame.
        if ~isempty(bearings)
            bearings=(90-bearings(Igood))*pi/180;
            for I=1:length(Igood)
                XX=DASAR_coords(Igood(I),1)+[0 linel*cos(bearings(I))]-xg;
                YY=DASAR_coords(Igood(I),2)+[0 linel*sin(bearings(I))]-yg;
                line(XX,YY);
            end
        end
        if ~isempty(VM)
            plot(VM(:,1)-xg,VM(:,2)-yg,'ko','markerfacecolor',[0 0 0],'markersize',5);
        end
        
        
        %Plot error elipps
        if ~isempty(A)
            %ELLIPSE(ra,rb,ang,x0,y0)
            h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
            set(h,'linewidth',0.5);
        end
        
        hold off;
        
    end


end
