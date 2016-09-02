%%%%%plot_site_histogram.m%%%%%%%%%
%%% interactive program to make histogram of calls detected per site,
%%% location, and detection stage.  Plots calls per day and per hour
%  Run from Arctic_2007/Processed folder

clear all
det_options={'raw','noairgun','classified'};

site_names=dir('Site*');
site_names=site_names(1:5);
hr_inc=datenum(0,0,0,1,0,0);
day_inc=datenum(0,0,1,0,0,0);
current_year=2007;
date_limits=[datenum(2007,8,23,0,0,0) datenum(2007,10,10,0,0,0)];

%Ichc=menu('Select a site:',site_names(:).name);
for Isite=1:length(site_names),
    site_name=site_names(Isite).name;
    cd(site_name);

    loc_names=dir('D*');
    cd('FinalResults');

    %Ichc=menu('Select a location:',loc_names(:).name);
    for Iloc=1:length(loc_names),
        loc_name=loc_names(Iloc).name;
        close all;
        %Ichc=menu('Which detection option?',det_options);
        for Ichc=1:3,
            clear tabs ncount ncount_hr tabs_hr
            tabs_hr=[];ncount_hr=[];
            det_str=det_options{Ichc};
            if strcmp(det_str,'classified')
                eval_str='best_ctimes';
            elseif strcmp(det_str,'raw')
                eval_str='best_ctimes_raw';
            elseif strcmp(det_str,'noairgun')
                eval_str='best_ctimes_noairgun';
            end

            fnames=dir([loc_name '*' det_str '.mat']);
            ncount=zeros(length(fnames),1);tabs=ncount;
            for I=1:length(fnames),
                load(fnames(I).name);
                ctimes=eval(eval_str);
                times=datenum(1970,1,1,0,0,ctimes);
                It=max(findstr('T',fnames(I).name)-4);
                tabs(I)=datenum([fnames(I).name(It+(0:1)) '/' fnames(I).name(It+(2:3)) '/' num2str(current_year)]);
                tabs_hr_day= tabs(I):hr_inc:(tabs(I)+day_inc);
                ncount_hr_day=histc(times,tabs_hr_day);

                tabs_hr_day=tabs_hr_day(1:(end-1));
                ncount_hr_day=ncount_hr_day(1:(end-1));

                tabs_hr=[tabs_hr tabs_hr_day];
                ncount(I)=length(ctimes);
                ncount_hr=[ncount_hr ncount_hr_day];
            end

            try,
                figure(1)
                subplot(3,1,Ichc);
                bar(tabs,ncount);xlim(date_limits);datetick('x',6,'keeplimits');grid on;ylabel('Calls/day');
                title(sprintf('%s %s %s',site_name,loc_name,det_str));
                orient tall
                print('-djpeg',sprintf('../../Plots.dir/%s_%s_day',site_name,loc_name));


                figure(2)
                subplot(3,1,Ichc);
                bar(tabs_hr,ncount_hr);xlim(date_limits);datetick('x',6,'keeplimits');grid on;ylabel('Calls/hr');
                title(sprintf('%s %s %s',site_name,loc_name,det_str));
                orient tall
                print('-djpeg',sprintf('../../Plots.dir/%s_%s_hour',site_name,loc_name));
            catch,
                disp('Could not plot, most likely empty count matrix');
            end
            %keyboard;

        end %Ichc

    end
    cd ../..
end