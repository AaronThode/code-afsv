correct_manual_results.m                                                                            0000677 0002527 0000000 00000024563 11002026333 016240  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        
%%Correct_manual_results.m%%
% Run after running evaluate_detector.m
%% Creates ctime vectors of calls that are actually noise, and vice versa

clear noise
Iloc=1;
Iday=menu('Select a day for detailed examination:',dates_list);
site_number=str2num(DASAR_list{Iloc}((end-1)));
if strcmp(lower(computer),'maci')
    dir_loc=sprintf('/Volumes/macmussel1/Arctic_2007/Site_0%i/%s',site_number,DASAR_list{Iloc});
elseif strcmp(lower(computer),'mac')
    dir_loc='../DataFiles.dir';
end

fnn=dir([dir_loc '/' DASAR_list{Iloc} '*' dates{Iday} '*sio']);
fn=[dir_loc '/' fnn.name];

Ichc=menu('What to review?','Missed calls','False detections','Successes','Airguns','Feature comparison');
%Isave=menu('Save time series?','Yes','No');
switch Ichc,
    case 2,
        %%%Review false detections ...
        I_false_classification=[];
        I_false=[];
        I_ctimes=[];
        false_times=times_false{Iday,Iloc};
        durations=best_calls{Iday,Iloc}.duration(Ix_nomatch)/param.Fs;
        try,
        features_false=best_calls{Iday,Iloc}.features(Ix_nomatch);
        BWfinal_false=best_calls{Iday,Iloc}.BWfinal(Ix_nomatch);
        end
        % false_times_mismatch=times_false_mismatch{Iday,Iloc};
        fsave_noise=[fnn.name(1:(end-4)) '_noise'];
        fsave_true=[fnn.name(1:(end-4)) '_FalseCallReview_' datestr(now,30)];

        prompt={'Enter a starting index:','Number of seconds to show in spectrogram'};
        name='TSV files correction: false detections that are calls';
        numlines=1;
        defaultanswer={'2',num2str(1+2*param.energy.bufferTime)};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        Imin=str2num(answer{1});
        tlen=str2num(answer{2});  %Number of seconds to show in spectrogram


        Isample=1;
        for I=Imin:(length(false_times)-1),
            if rem(I,10)==0,
                disp(['Saving ' fsave_true]);
                save(fsave_true,'I_false_classification','I_false','I_ctimes');

            end

            try,
                tlen=2*param.energy.bufferTime+durations(I);
                [x,t,head]=readsiof(fn,false_times(I)-param.energy.bufferTime,tlen);
           
                
                nearest_manual=min([ abs(false_times(I)-times_missed{Iday,Iloc}); abs(false_times(I)-times_good{Iday,Iloc})]);
                [S,F,T,PP]=spectrogram(x,128,96,128,head.Fs,'yaxis');
                figure(3);
                imagesc(T,F,10*log10(abs(PP)));axis('xy');
                caxis([-20 50]);

                %%Check if notes exist on this false detection.
                ID=-1;
                if ~isempty(false_notes)
                    [junk,Inote]=min(abs(false_times(I)-false_notes.ctimes));
                    if junk<=tol,
                        ID=false_notes.reason(Inote);
                    end
                end
                title(sprintf('%s: %i %s False Detection %i of %i at %6.2f seconds, mismatch time: %6.2f, prev false: %6.2f s, next false: %6.2f s', ...
                    fn,ID,datestr(datenum(1970,1,1,0,0,false_times(I))),I,length(false_times), ...
                    false_times(I)-head.ctbc,nearest_manual, ...
                    false_times(I)-false_times(I-1),false_times(I+1)-false_times(I)));

                %[call,Bmean2]=contour_postprocessor(x, false_times(I),param,0,[]);  %Bmean is reset
                
                if 1==0,
                [features,BWfinal,Bmean2]=contour_postprocessor_morph(x, false_times(I),param,1,[]);  %Bmean is reset

                figure(2)
                subplot(3,1,3);
                imshow(BWfinal_false{I});title('original stored image');
                end
                
                %yes=menu('Sound?','Yes','No');
                %if yes==1,
                %    soundsc(x,head.Fs);
                %end

                tmp=menu('What is the decision? ', 'True call','Song component','Airgun','Weak airgun','Seal','Walrus','Tone','False, nothing distinctive');


                I_false_classification=[I_false_classification tmp];
                I_false=[I_false I];
                I_ctimes=[I_ctimes false_times(I)];


            catch,
                disp('Program has failed, contact Aaron.');
                keyboard
            end
        end


    case 5,  %Feature plot...
        plot_features(best_calls{Iday,Iloc},I_false_alarms{Iday,Iloc},I_auto_match{Iday,Iloc},I_missed_calls{Iday,Iloc});
        %REVIEW MISSED CALLS
    case 1,
        %fn='D07s3bT20070901T000000.sio';
        Icontour=menu('Run contour tracer?','Yes','No');
        I_really_found=[];
        missed_times=times_missed{Iday,Iloc};

        %%perhaps airgun times scooped these up?
        %airgun_times=airgun_ctimes{Iday,Iloc};

        tlen=param.energy.bufferTime;
        fsave_noise=[fnn.name(1:(end-4)) '_noise'];
        fsave_weak=[fnn.name(1:(end-4)) '_weakcalls'];

        prompt={'Enter a starting index:','Leave this blank:'};
        name='TSV files correction: marked calls that are actually noise';
        numlines=1;
        defaultanswer={'1','dummy'};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        Imin=str2num(answer{1});

        for I=Imin:length(missed_times),
            if rem(I,25)==0&length(I_really_found)>1,
                disp(['Saving ' fsave_weak]);
                best_ctimes=times_missed{Iday,Iloc}(I_really_found);
                save(fsave_weak,'missed_times');

            end

            %%Check to see if close to airgun time
            %[near_gun_time(I),Igun_close]=min(abs(missed_times(I)-airgun_times));
            [x,t,head]=readsiof(fn,missed_times(I)-tlen,tlen+4);
            figure(1);
            [S,F,T,PP]=spectrogram(x,128,96,128,head.Fs,'yaxis');
            imagesc(T,F,10*log10(abs(PP)));axis('xy');
            caxis([-20 50]);

            %caxis([0 60]);
            title(sprintf('%s: Missed Call %i of %i at %6.2f seconds into file,  %s', ...
                fn,I,length(missed_times),missed_times(I)-head.ctbc,datestr(datenum(1970,1,1,0,0,missed_times(I)))));

            [call,Bmean2]=contour_postprocessor(x, missed_times(I),param,1,[]);  %Bmean is reset

            tmp=menu('True call? ', 'Yes','No');

            if tmp==2,
                %I_false_classification=[I_false_classification I];
            else
                I_really_found=[I_really_found I];
                if Isave==1,
                    weakcall{Isample}.om=x;
                    Isample=Isample+1;
                end

            end
        end
        %tmp=missed_times(I_really_missed);
        tmp=[times_good{Iday,Iloc}' times_missed{Iday,Iloc}(I_really_found)']';
        tmp=sort(tmp);
        Igood=find(diff(tmp)>tol);
        manual{Iday,Iloc}=tmp(Igood);

        %save D07s3*T2007****T000000_missed_corrections revised times

    case 3,
        %%%And what worked? ...
        good_times=times_good{Iday,Iloc};
        Ix_match=I_auto_match{Iday,Iloc};
        tlen=2*param.energy.bufferTime+best_calls{Iday,Iloc}.duration(1,Ix_match);

        Istrip=round(param.Fs*param.energy.bufferTime);

        for I=1:length(good_times),
            [x,t,head]=readsiof(fn,good_times(I)-param.energy.bufferTime,tlen(I));

            figure(1);
            [S,F,T,PP]=spectrogram(x,128,96,128,head.Fs,'yaxis');
            imagesc(T,F,10*log10(abs(PP)));axis('xy');
            caxis([-20 50]);
            index=1:(length(x)-Istrip);

            %caxis([0 60]);
            title(sprintf('%s: direct load Good Detection %i of %i at %6.2f seconds, %s',fn,I,length(good_times), ...
                good_times(I)-head.ctbc, datestr(datenum(1970,1,1,0,0,good_times(I)))));
            %[call,Bmean2]=contour_postprocessor(x(index), good_times(I),param,1,[]);  %Bmean is reset

            %             snips_name=dir('*snips');
            %             IND=best_calls{Iday,Iloc}.index(Ix_match(I));
            %             [x,nstarts_snips,npts_snips]=readEnergySnips(snips_name.name, IND,'short','cell');
            %             figure(1);
            %             spectrogram(x{1},128,96,128,head.Fs,'yaxis');
            %             index=1:(length(x{1})-Istrip);
            %
            %             %caxis([0 60]);
            %             title(sprintf('%s: snips data Good Detection %i of %i at %6.2f seconds',fn,I,length(good_times),good_times(I)-head.ctbc));
            %             [call,Bmean2]=contour_postprocessor(x{1}(index), good_times(I),param,1,[]);  %Bmean is reset
            pause;

        end
    case 4,
        %%%Review airguns...
        I_really_airgun=[];
        %actual_airgun_ctimes=[];
        airgun_times=airgun_ctime{Iday,Iloc};
        % airgun_times_mismatch=times_airgun_mismatch{Iday,Iloc};
        tlen=4;
        Icount=1;
        for I=2:(length(airgun_times)-1),
            try,
                [x,t,head]=readsiof(fn,airgun_times(I)-tlen/2,tlen);
                nearest_manual=min(abs(airgun_times(I)-times_good{Iday,Iloc}));
                spectrogram(x,128,96,128,head.Fs,'yaxis');
                %caxis([0 60]);
                figure(gcf)
                title(sprintf('%s: Airgun %i of %i at %6.2f seconds,  prev airgun: %6.2f s, next airgun: %6.2f s', ...
                    fn,I,length(airgun_times),airgun_times(I)-head.ctbc,airgun_times(I)-airgun_times(I-1),airgun_times(I+1)-airgun_times(I)));
                tmp=input('Enter 1 if not a airgun alarm:  ');
                if isempty(tmp),
                    airgun{Icount}=x;
                    Icount=Icount+1;
                    %I_really_airgun=[I_really_airgun I];
                else
                    % actual_call_ctimes=[actual_call_ctimes airgun_times(I)];
                end
            end

        end

        %Optional creation of airgun files for optimization
        Icount=1;

        for I=[1:50 2000:2050 5000:5050 7000:7050 10000:10050], %9/18
            %for I=[1:50 1000:1050 2000:2050 3000:3050 4000:4050], %8/30
            try,
                [x,t,head]=readsiof(fn,airgun_times(I)-tlen/2,tlen);
                airgun{Icount}.data=x;
                airgun{Icount}.ctime=airgun_times(I)-tlen/2;
                airgun{Icount}.Fs=1000;
                Icount=Icount+1;
                %I_really_airgun=[I_really_airgun I];

            end
        end

        Islsh=max(findstr(fn,'/'))+1;
        fsave=fn(Islsh:(end-4));
        save(fsave,'airgun');
        % tmp=[times_good{Iday,Iloc}' times_missed{Iday,Iloc}' actual_call_ctimes];
        % tmp=sort(tmp);
        % Igood=find(diff(tmp)>tol);
        % times_good_reviewed=tmp(Igood);
end
%What is sigdb of missed calls?

                                                                                                                                             evaluate_detector.m                                                                                 0000677 0002527 0000000 00000023442 11002022602 015146  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        %%%comparecalls_revised.m%
clear all;close all
path(path,'../CommonScripts.dir');

DASAR_list={'D07s3a','D07s3b','D07s3c','D07s3d','D07s3e','D07s3f','D07s3g', ...
    'D07s4a','D07s4b','D07s4c','D07s4d','D07s4e','D07s4f','D07s4g',...
    'D07s5a','D07s5b','D07s5c','D07s5d','D07s5e','D07s5f','D07s5g'};
dates_list={'0823','0824','0825','0826','0827', ...
    '0828','0829','0830','0831','0901', ....
    '0902', '0903','0904','0905','0906', ...
    '0907','0908','0909', '0910','0911', ...
    '0912','0913','0914','0915','0916'...
    '0917','0918','0919','0920','0923', ...
    '0924','0925','0926','0927','0929', ...
    '1001','1002','1003','1004','1005', ...
    '1006','1007','1008','1009','1010'
    };

Icase='Site4_DASARg';
detection_stage={'raw','classified','random'}; %raw, classififed, crosschecked or random.
%The last makes a list of random ctimes to compare performance.
% Crosschecked uses results used from checking multiple DASARS. Note
%       best_calls will not have same legnth as auto_ctimes
ICI_strip=0; %Try to strip airgun repetitive patterns
feature_strip=0;  %Use morphological features to filter detections after ICI stripping...
tol=4;

%%for ICI analysis
Ndet=5;
Nmiss=1;
ICItol=0.5;

switch Icase,
    case 'Site3_DASARa',
        detection_stage={'classified'}; %raw, classififed, crosschecked,or random.
        %The last makes a list of random ctimes to compare performance

        DASAR_list=DASAR_list(1);
        dates_list=dates_list(1:(end-2));

        dates_list=dates_list(1:10);
        %dates_list=dates_list(7);  %August 29, 2007

        %
        %         Ichc=menu('Select day:','8/27','9/18');
        %         if Ichc==1,
        %             dates_list=dates_list(5); %Augsut 27, 2007
        %         else
        %             dates_list=dates_list(27); %Sept. 18, 2007
        %         end
        %         clear Ichc

        autoCase='Morph1Site3';
        %autoCase='CrossChannelSite3';
    case 'Site4_DASARg',
        detection_stage={'classified'}; %raw, classififed, crosschecked,or random.
        %The last makes a list of random ctimes to compare performance

        DASAR_list=DASAR_list(14)
        %dates_list=dates_list(1:21);
        dates_list=dates_list(1:(end-2));

        Ichc=menu('Select day:','8/30','9/20','all days');
        if Ichc==1,
            dates_list=dates_list(8); %Augsut 27, 2007
        elseif Ichc==2,
            dates_list=dates_list(29); %Sept. 18, 2007
        end
        clear Ichc

        %Ichc=menu('Select type:','Contour','Morphological');
        %if Ichc==1,
        autoCase='Contour1Site4';
        disp('Contour option');
        %else
        %    autoCase='Morph1Site4';
        %    disp('Morph option');
        %end
        %clear Ichc
        %
        %autoCase='CrossChannelSite4';
    case 'Site5_DASARa',
        detection_stage={'classified'}; %raw, classififed, crosschecked,or random.
        %The last makes a list of random ctimes to compare performance

        DASAR_list=DASAR_list(15);
        % dates_list=dates_list(1:21);

        %dates_list=dates_list(1:(end-2))
        Ichc=menu('Select day:','8/27','9/18');

        %dates_list=dates_list(7);  %August 29, 2007
        if Ichc==1,
            dates_list=dates_list(5); %Augsut 27, 2007
        else
            dates_list=dates_list(27); %Sept. 18, 2007
        end
        clear Ichc
        dates_list

        autoCase='CrossChannelSite5';
    case 'Site4_Sept18'
        %DASAR_list=DASAR_list(8:end); %Site 4
        DASAR_list=DASAR_list((end)); %Site 4
        detection_stage={'classified'}; %raw, classififed, crosschecked, or random.  The last makes a list of random ctimes to compare performance
        %dates_list=dates_list(27)
        autoCase='CrossChannelSept18Site4'
        %keyboard;
end
force_manual_reload=1;  %If 1, reload, if 2, load stored.  If in doubt, use 1


mydir=pwd;

[manual,SNR,dates,false_notes]=load_manual_results(force_manual_reload, DASAR_list,dates_list);
dbstop if error
for Istage=1:length(detection_stage),
    disp(sprintf('**Testing detection stage: %s',detection_stage{Istage}));
    [auto_ctimes,best_calls,airgun_calls,param]=load_automated_results(autoCase,DASAR_list,detection_stage{Istage},dates);
    %Ndet=param.interval_remove.Ndet;
    %ICItol=param.interval_remove.ICItol;
    %random_ctimes=rand(size(auto_ctimes));
    %Test whether ICI can be identified

    total.missed=0;
    total.false=0;
    total.correct=0;
    for Iloc=1:length(DASAR_list),  %DASAR loop
        for Iday=1:length(dates_list),  %Day loop

            %%Clean up manual ctimes
            %Igood=find(manual{Iday}.wctype(:,Iloc)>0);
            %manual_ctime=manual{Iday}.ctime(Igood,Iloc);
            disp(sprintf('Day: %s DASAR: %s',dates{Iday},DASAR_list{Iloc}));
            manual_ctime=manual{Iday,Iloc};

            %%Restrict to certain SNR ratios...

            auto_ctime=auto_ctimes{Iday,Iloc};
            if ICI_strip==1,
                disp('Stripping regular ICI values');

                %ICI=compute_ici_forward(datenum(1970,1,1,0,0,auto_ctime),[3 40],Ndet,ICItol);
                %keyboard;
                %auto_ctime=best_calls{Iday,Iloc}.features.Centroid;
                %  feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features,{'Centroid'}, ...
                %    1,1);
                %auto_ctime=auto_ctime+feature_vector{1}.';
                ICI=compute_ici_bothways(datenum(1970,1,1,0,0,auto_ctime),[2 30],Ndet,Nmiss,ICItol);
                
                Iguns=find(ICI>0);
                airgun_ctime{Iday,Iloc}=auto_ctime(Iguns);
                %keyboard;
                % airgun_calls{Iday,Iloc}.p=best_calls{Iday,Iloc}.p;
                Imatch=find(ICI<0);
                auto_ctime=auto_ctime(Imatch);
                best_calls{Iday,Iloc}.features=best_calls{Iday,Iloc}.features(Imatch);
                best_calls{Iday,Iloc}.ctime=best_calls{Iday,Iloc}.ctime(Imatch);
                best_calls{Iday,Iloc}.BWfinal=best_calls{Iday,Iloc}.BWfinal(Imatch);
            end

            if feature_strip==1,
                disp('Stripping features from Morphological');

                param=TOC_params('Feb14_2008_0830Site3a');
                
                disp('filtering AND rules');
                Iand=find(param.feature.operator>0);
                feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features,param.feature.feature_name(Iand), ...
                    param.feature.feature_index(Iand),param.feature.Nshapes);
                [best_filtered, auto_ctime,Imatch]=filter_feature_vector(feature_vector,param.feature.feature_name(Iand), ...
                    auto_ctime, param.feature.optvec(Iand,:),'first','and');
                
                
                 disp('filtering OR rules');
                Ior=find(param.feature.operator==0);
                feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features(Imatch),param.feature.feature_name(Ior), ...
                    param.feature.feature_index(Ior),param.feature.Nshapes);
               
                [best_filtered, auto_ctime,Imatchor]=filter_feature_vector(feature_vector,param.feature.feature_name(Ior), ...
                    auto_ctime, param.feature.optvec(Ior,:),'first','or');
                
                Imatch=Imatch(Imatchor);
                
                
                best_calls{Iday,Iloc}.features=best_calls{Iday,Iloc}.features(Imatch);
                best_calls{Iday,Iloc}.ctime=best_calls{Iday,Iloc}.ctime(Imatch);
                best_calls{Iday,Iloc}.BWfinal=best_calls{Iday,Iloc}.BWfinal(Imatch);
               
                
                %feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features,param.feature.feature_name, ...
                  %  param.feature.feature_index,param.feature.Nshapes);
              
            end


            [Ix_nomatch,Ix_match,Iy_match,Iy_nomatch,Ix_mismatch_time,Nxx]=find_similar_elements(auto_ctime,manual_ctime,tol,tol);
            I_false_alarms{Iday,Iloc}=Ix_nomatch;
            I_auto_match{Iday,Iloc}=Ix_match;
            I_manual_match{Iday,Iloc}=Iy_match;
            I_missed_calls{Iday,Iloc}=Iy_nomatch;
            Nmanual(Iday,Iloc)=length(manual_ctime);

            total.missed(Iloc)=total.missed(Iloc)+length(Iy_nomatch);
            total.false(Iloc)=total.false(Iloc)+length(Ix_nomatch);
            total.correct(Iloc)=total.correct(Iloc)+length(Ix_match);

            missed_call_fraction(Iday,Iloc)=length(Iy_nomatch)/length(manual_ctime);
            false_alarm_fraction(Iday,Iloc)=length(Ix_nomatch)/(length(Ix_nomatch)+length(Ix_match));  %Ignores auto_ctime within tol
            success_rate(Iday,Iloc)=length(Iy_match)/(length(Ix_nomatch)+length(Ix_match));

            times_missed{Iday,Iloc}=manual_ctime(Iy_nomatch);
            times_false{Iday,Iloc}=auto_ctime(Ix_nomatch);
            times_false_mismatch{Iday,Iloc}=Ix_mismatch_time;
            times_good{Iday,Iloc}=manual_ctime(Iy_match);

            %snr_missed{Iday,Iloc}=manual{Iday}.stndb(Igood(Iy_nomatch),Iloc);
            %snr_found{Iday,Iloc}=manual{Iday}.stndb(Igood(Iy_match),Iloc);
            %disp(sprintf('Mean snr of missed call: %6.2f mean snr of found call: %6.2f',median(snr_missed{Iday,Iloc}),median(snr_found{Iday,Iloc})));

        end %Iday
        disp(sprintf('DASAR: %s  date: %s  auto stage: %s',DASAR_list{Iloc},dates{Iday},detection_stage{Istage}));
        disp(sprintf('Tolerance: %6.2f, Number matched calls, number missed calls, number false alarms',tol));
        [dates' I_auto_match I_missed_calls I_false_alarms]
        disp('missed call fraction, false alarm call fraction');
        [ missed_call_fraction(:,Iloc) false_alarm_fraction(:,Iloc)]
    end %J
    total.missed_percent=total.missed./(total.missed+total.correct);
    total.false_percent=total.false./(total.false+total.correct);
    disp('total missed and false fraction for all days');
    total
end %Istage


                                                                                                                                                                                                                              get_feature_names.m                                                                                 0000644 0002527 0000000 00000001206 10772540467 015143  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        %%%Obtain feature names and secondary index if needed..
function [feature_name, feature_index]=get_feature_names(features);

prompt='Number of features:';
answer=inputdlg({'Number of features'},'Feature Name Selector',1,{'4'});
Nfeatures=str2num(answer{1});
names=fieldnames(features{1});
%Nfeatures=4;
for JJ=1:Nfeatures,
    Ichc=menu(sprintf('Select Feature %i',JJ),names);
    feature_name{JJ}=names{Ichc};
    test_value=getfield(features{1},feature_name{JJ});
    if sum(size(test_value))==2,
        feature_index(JJ)=1;
    else
        feature_index(JJ)=input(sprintf('Select index between 1 and %i: ',length(test_value)));
    end


end                                                                                                                                                                                                                                                                                                                                                                                          load_automated_results.m                                                                            0000677 0002527 0000000 00000007422 11002022306 016213  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        
%%%load_automated_results.m%%%
%function auto_ctimes=load_automated_results(Icase,DASAR_list,detect_str,dates),
%   Input: detect_str: 'classified','raw','random'
%     Icase, a number representing a particular scenario to load
%     DASAR_list, a cell matrix of DASAR names
%     dates:  a cell matrix of dates in form '0928'
%     detect_str: 'classified','raw','etc'.
%  Output: auto_ctimes;
%     auto_ctimes(length(dates),length(DASAR_locs));

function [auto_ctimes,best_calls,airgun_calls,param]=load_automated_results(Icase,DASAR_list,detect_str,dates),
airgun_calls=[];
if strcmp(lower(computer),'mac')
    loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection';
elseif strcmp(lower(computer),'maci')
    loc.base='/Volumes/macmussel1/Arctic_2007/Processed';
end

switch Icase,
    case 'Morph1Site4',
        if strcmp(lower(computer),'mac')
            loc.automated=[loc.base '/Detections_2007.dir/Run5.dir'];
        elseif strcmp(lower(computer),'maci')
            loc.automated=[loc.base '/Site_04/FinalResults_morph' ];

        end
    case 'Contour1Site4',
        if strcmp(lower(computer),'mac')
            loc.automated=[loc.base '/Detections_2007.dir/Run5.dir'];
        elseif strcmp(lower(computer),'maci')
            loc.automated=[loc.base '/Site_04/FinalResults_contour' ];

        end
        %keyboard;
    case 'Morph1Site3',
        if strcmp(lower(computer),'mac')
            loc.automated=[loc.base '/Detections_2007.dir/Run5.dir'];
        elseif strcmp(lower(computer),'maci')
            loc.automated=[loc.base '/Site_03/FinalResults_morph' ];

        end
    case 'CrossChannelSite5',
        if strcmp(lower(computer),'mac')
            loc.automated=[loc.base '/Detections_2007.dir/Run5.dir'];
        elseif strcmp(lower(computer),'maci')
            loc.automated=[loc.base '/Site_05/FinalResults2' ];

        end

    case 'CrossChannelSite4',
        if strcmp(lower(computer),'mac')
            loc.automated=[loc.base '/Detections_2007.dir/Run5.dir'];
        elseif strcmp(lower(computer),'maci')
            loc.automated=[loc.base '/Site_04/FinalResults_contour' ];

        end
end
disp(sprintf('loading automated detection results from %s',loc.automated));
if isempty(loc.automated)
    disp('keyboard not found in load_automated_results.m');
end
% fnames_auto=dir([loc.automated '/*' detect_str '.mat']);
for I=1:length(dates),

    if isempty(dates{I}),
        continue;
    end
    for J=1:length(DASAR_list),
        if strcmp(detect_str,'random')|strcmp(detect_str,'crosschecked')
            search_str=[DASAR_list{J} '*2007' dates{I}  '*' 'classified.mat'] ;
        else
            search_str=[DASAR_list{J} '*2007' dates{I}  '*' detect_str '.mat'] ;
        end
        fname_auto=dir([loc.automated '/' search_str ]);
        if length(fname_auto)>1,
            error(sprintf('More than one automated detection file for %s',search_str));
        end
        fname_auto=[loc.automated '/' fname_auto.name];
        auto=load(fname_auto);
        if isfield(auto.best_calls,'level')
            auto.best_calls=rmfield(auto.best_calls,'level');
        end
        if strcmp(detect_str,'classified')
            auto_ctimes{I,J}=auto.best_ctimes;
            best_calls{I,J}=auto.best_calls;
            param=auto.param;
        elseif strcmp(detect_str,'crosschecked')
            auto_ctimes{I,J}=auto.best_ctimes_cross;
            best_calls{I,J}=auto.best_calls;
            keyboard;
            param=auto.param;
        elseif strcmp(detect_str,'raw')
            auto_ctimes{I,J}=auto.best_ctimes_raw;
            param=auto.param;
            best_calls{I,J}=[];

        elseif strcmp(detect_str,'random')
            auto_ctimes{I,J}=sort(min(auto.best_ctimes)+(max(auto.best_ctimes)-min(auto.best_ctimes))*rand(size(auto.best_ctimes)));
        end

    end
end


end                                                                                                                                                                                                                                              manual_select_thresholds.m                                                                          0000644 0002527 0000000 00000000756 10774771344 016554  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        function optvec=manual_select_thresholds(feature_name);
%%%Restrict ranges and recompute successes...
for JJ=1:length(feature_name),
    prompt{(JJ-1)*2+1}=sprintf('%s min value',feature_name{JJ});
    prompt{(JJ-1)*2+2}=sprintf('%s max value',feature_name{JJ});
    defaultanswer{(JJ-1)*2+1}='0';
    defaultanswer{(JJ-1)*2+2}='1';
end
name='Feature pruning';
numlines=1;

answer=inputdlg(prompt,name,numlines,defaultanswer);

for I=1:length(answer);
    optvec(I)=str2num(answer{I});
end

end                  pathdef.m                                                                                           0000644 0002527 0000000 00000012151 10770541726 013077  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        function p = pathdef
%PATHDEF Search path defaults.
%   PATHDEF returns a string that can be used as input to MATLABPATH
%   in order to set the path.

  
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $ $Date: 2007/06/07 14:45:14 $


% DO NOT MODIFY THIS FILE.  IT IS AN AUTOGENERATED FILE.  
% EDITING MAY CAUSE THE FILE TO BECOME UNREADABLE TO 
% THE PATHTOOL AND THE INSTALLER.

p = [...
%%% BEGIN ENTRIES %%%
     matlabroot,'/toolbox/matlab/general:', ...
     matlabroot,'/toolbox/matlab/ops:', ...
     matlabroot,'/toolbox/matlab/lang:', ...
     matlabroot,'/toolbox/matlab/elmat:', ...
     matlabroot,'/toolbox/matlab/elfun:', ...
     matlabroot,'/toolbox/matlab/specfun:', ...
     matlabroot,'/toolbox/matlab/matfun:', ...
     matlabroot,'/toolbox/matlab/datafun:', ...
     matlabroot,'/toolbox/matlab/polyfun:', ...
     matlabroot,'/toolbox/matlab/funfun:', ...
     matlabroot,'/toolbox/matlab/sparfun:', ...
     matlabroot,'/toolbox/matlab/scribe:', ...
     matlabroot,'/toolbox/matlab/graph2d:', ...
     matlabroot,'/toolbox/matlab/graph3d:', ...
     matlabroot,'/toolbox/matlab/specgraph:', ...
     matlabroot,'/toolbox/matlab/graphics:', ...
     matlabroot,'/toolbox/matlab/uitools:', ...
     matlabroot,'/toolbox/matlab/strfun:', ...
     matlabroot,'/toolbox/matlab/imagesci:', ...
     matlabroot,'/toolbox/matlab/iofun:', ...
     matlabroot,'/toolbox/matlab/audiovideo:', ...
     matlabroot,'/toolbox/matlab/timefun:', ...
     matlabroot,'/toolbox/matlab/datatypes:', ...
     matlabroot,'/toolbox/matlab/verctrl:', ...
     matlabroot,'/toolbox/matlab/codetools:', ...
     matlabroot,'/toolbox/matlab/helptools:', ...
     matlabroot,'/toolbox/matlab/demos:', ...
     matlabroot,'/toolbox/matlab/timeseries:', ...
     matlabroot,'/toolbox/matlab/hds:', ...
     matlabroot,'/toolbox/matlab/guide:', ...
     matlabroot,'/toolbox/matlab/plottools:', ...
     matlabroot,'/toolbox/local:', ...
     matlabroot,'/toolbox/shared/controllib:', ...
     matlabroot,'/toolbox/dspblks/dspblks:', ...
     matlabroot,'/toolbox/dspblks/dspmasks:', ...
     matlabroot,'/toolbox/dspblks/dspmex:', ...
     matlabroot,'/toolbox/dspblks/dspdemos:', ...
     matlabroot,'/toolbox/shared/filterdesignlib:', ...
     matlabroot,'/toolbox/gads:', ...
     matlabroot,'/toolbox/gads/gads:', ...
     matlabroot,'/toolbox/gads/gadsdemos:', ...
     matlabroot,'/toolbox/images:', ...
     matlabroot,'/toolbox/images/icons:', ...
     matlabroot,'/toolbox/images/images:', ...
     matlabroot,'/toolbox/images/imdemos:', ...
     matlabroot,'/toolbox/images/imdemos/demosearch:', ...
     matlabroot,'/toolbox/images/imdemos/html:', ...
     matlabroot,'/toolbox/images/imuitools:', ...
     matlabroot,'/toolbox/images/iptformats:', ...
     matlabroot,'/toolbox/images/iptutils:', ...
     matlabroot,'/toolbox/images/medformats:', ...
     matlabroot,'/toolbox/nnet:', ...
     matlabroot,'/toolbox/nnet/nncontrol:', ...
     matlabroot,'/toolbox/nnet/nndemos:', ...
     matlabroot,'/toolbox/nnet/nnet:', ...
     matlabroot,'/toolbox/nnet/nnet/nnanalyze:', ...
     matlabroot,'/toolbox/nnet/nnet/nncustom:', ...
     matlabroot,'/toolbox/nnet/nnet/nndistance:', ...
     matlabroot,'/toolbox/nnet/nnet/nnformat:', ...
     matlabroot,'/toolbox/nnet/nnet/nninit:', ...
     matlabroot,'/toolbox/nnet/nnet/nnlearn:', ...
     matlabroot,'/toolbox/nnet/nnet/nnnetinput:', ...
     matlabroot,'/toolbox/nnet/nnet/nnnetwork:', ...
     matlabroot,'/toolbox/nnet/nnet/nnperformance:', ...
     matlabroot,'/toolbox/nnet/nnet/nnplot:', ...
     matlabroot,'/toolbox/nnet/nnet/nnprocess:', ...
     matlabroot,'/toolbox/nnet/nnet/nnsearch:', ...
     matlabroot,'/toolbox/nnet/nnet/nntopology:', ...
     matlabroot,'/toolbox/nnet/nnet/nntrain:', ...
     matlabroot,'/toolbox/nnet/nnet/nntransfer:', ...
     matlabroot,'/toolbox/nnet/nnet/nnweight:', ...
     matlabroot,'/toolbox/nnet/nnguis:', ...
     matlabroot,'/toolbox/nnet/nnguis/nftool:', ...
     matlabroot,'/toolbox/nnet/nnguis/nntool:', ...
     matlabroot,'/toolbox/nnet/nnobsolete:', ...
     matlabroot,'/toolbox/nnet/nnresource:', ...
     matlabroot,'/toolbox/nnet/nnutils:', ...
     matlabroot,'/toolbox/optim/optim:', ...
     matlabroot,'/toolbox/optim/optimdemos:', ...
     matlabroot,'/toolbox/shared/optimlib:', ...
     matlabroot,'/toolbox/signal/signal:', ...
     matlabroot,'/toolbox/signal/sigtools:', ...
     matlabroot,'/toolbox/signal/sptoolgui:', ...
     matlabroot,'/toolbox/signal/sigdemos:', ...
     matlabroot,'/toolbox/shared/spcuilib:', ...
     matlabroot,'/toolbox/shared/dastudio:', ...
     matlabroot,'/toolbox/wavelet:', ...
     matlabroot,'/toolbox/wavelet/wavedemo:', ...
     matlabroot,'/toolbox/wavelet/wavedemo/demosearch:', ...
     matlabroot,'/toolbox/wavelet/wavedemo/html:', ...
     matlabroot,'/toolbox/wavelet/wavelet:', ...
     matlabroot,'/toolbox/wavelet/wmultisig1d:', ...
     matlabroot,'/work:', ...
     '/Users/Shared/DataPreProcessor/MATLAB_IO:', ...
     '/Users/Shared/at/MATLAB:', ...
     '/Users/Shared/Projects/Arctic_2007/CommonScripts.dir:', ...
     matlabroot,'/toolbox/shared/imageslib:', ...
%%% END ENTRIES %%%
     ...
];

p = [userpath,p];
                                                                                                                                                                                                                                                                                                                                                                                                                       plot_features.m                                                                                     0000644 0002527 0000000 00000016525 11000040324 014320  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        %%%plot and compare detection features against each other%%%%
% function plot_features(best_calls,I_false_alarms,I_auto_match,I_missed_calls);

function plot_features(best_calls,I_false_alarms,I_auto_match,I_missed_calls);

features=best_calls.features;


%%%Obtain feature names and secondary index if needed..
[feature_name, feature_index]=get_feature_names(features);

%%Extract feature vector to level desired...
%keyboard;
Nshapes=3;  %Number of shapes per detection to evaluate...
p.false=extract_feature_vector(features(I_false_alarms),feature_name,feature_index,Nshapes);

p.match=extract_feature_vector(features(I_auto_match),feature_name,feature_index,Nshapes);

p.missed=extract_feature_vector(features(I_missed_calls),feature_name,feature_index,Nshapes);



p.Nfeatures=length(feature_name);
p.feature_name=feature_name;

%%plot comparisons of features...
Nplots=p.Nfeatures-1;
colorstr='kgcym';
plotstr='oxs';
for Iplot=1:Nplots,
    for Ishape=1:Nshapes,
        figure(Ishape)
        ps=[colorstr(Ishape) 'o'];

        subplot(Nplots,1,Iplot);
        plot(p.false{1}(:,Ishape),p.false{Iplot+1}(:,Ishape),ps);hold on;
        plot(p.match{1}(:,Ishape),p.match{Iplot+1}(:,Ishape),['r' plotstr(Ishape)]);
        grid on
        xlabel(feature_name{1});ylabel(feature_name{Iplot+1});
        % keyboard;
    end
    hold off
end

figure
for Iplot=1:(Nplots+1),
    subplot(Nplots+1,2,2*Iplot-1);
    [NN,XX]=hist(p.false{Iplot},200);
    bar(XX,NN);
    grid on
    title(feature_name{Iplot});
    ylabel('False detections');

    subplot(Nplots+1,2,2*Iplot);
    hist(p.match{Iplot},XX);

    grid on
    title(feature_name{Iplot});
    ylabel('True detections');

    hold off
end

figure
if p.Nfeatures>2,
    plot3(p.false{1}(:),p.false{2}(:),p.false{3}(:),'ko',p.match{1}(:),p.match{2}(:),p.match{3}(:),'ro');grid on
    xlabel(feature_name{1});ylabel(feature_name{2});zlabel(feature_name{3});
end

winnow_chc=menu('How set discrimination?','Manually','Optimization');
if winnow_chc==1,
    optvec=manual_select_thresholds(p.feature_name);
    optvec=reshape(optvec,2,[])';
else
    disp('Insert default selection criterica here...');
end

keyboard;
%%Create the evaluation function...
eval_str={'first','all'};
debug=1;
for I=1:2,

    disp(sprintf('Using evaluation criteria: %s',eval_str{I}));
    disp('False alarm filter...');
    [p.false_filter, ctimes_false,p.false_pass_I,p.false_fail_I]=filter_feature_vector(p.false,feature_name,best_calls.ctime(I_false_alarms), ...
        optvec,eval_str{I},'and',debug);
    
   
    disp('True call filter...');
    [p.match_filter, ctimes_match,p.match_pass_I,p.match_fail_I]=filter_feature_vector(p.match,feature_name,best_calls.ctime(I_auto_match), ...
        optvec,eval_str{I},'and',debug);

    Nlost=length(p.match{1})-length(ctimes_match);
    Nmissed_new=length(p.missed{1})+Nlost;
    Nfalse_new=length(ctimes_false);
    Nmatch_new=length(ctimes_match);

    disp(sprintf('After discrimination: %i matched calls %i missed calls, %i false alarms',Nmatch_new,Nmissed_new,Nfalse_new));
    disp(sprintf('%i true calls lost, %i false alarms stripped\n\n\n',Nlost, length(p.false{1})-Nfalse_new));

    missed_fraction=Nmissed_new/(Nmissed_new+Nmatch_new);
    false_ratio=Nfalse_new/Nmatch_new;
    fit=1*missed_fraction+1*false_ratio;
end

keyboard;
% lost_missed=features(p.match_fail_I{1});
%  figure;
%    
% for J=1:length(lost_missed),
%     imshow(lost_missed{J}.Image);
%     pause;
% end
% keyboard;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TEXT TO DEMONSTRATE POLYNOMIAL FITS%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %I_false_alarms =Ix_nomatch;
%         %I_auto_match =Ix_match;
%         %I_manual_match =Iy_match;
%         %I_missed_calls =Iy_nomatch;
%
%         %Compare magnitudes...
%         %         magnitude_calls=10*log10(eps+sum(best_calls .magnitude(1:2,I_auto_match )));
%         %         magnitude_false=10*log10(eps+sum(best_calls .magnitude(1:2,I_false_alarms )));
%         %
%         Ifreq=1:10;
%         magnitude_calls=10*log10(eps+sum(best_calls .magnitude(Ifreq,I_auto_match )));
%         magnitude_false=10*log10(eps+sum(best_calls .magnitude(Ifreq,I_false_alarms )));
%
%
%         figure;
%         subplot(2,1,1);plot(10*log10(eps+sum(magnitude_calls)),'ro');hold on;grid on;ylim([50 140]);title('peak SEL');
%         subplot(2,1,2);plot(10*log10(eps+sum(magnitude_false)),'ko');grid on; ylim([50 140]);
%
%
%         duration_calls=(best_calls .duration(1,I_auto_match ));
%         duration_false=(best_calls .duration(1,I_false_alarms ));
%
%         figure;
%         subplot(2,1,1);plot(duration_calls,'ro');hold on;grid on;title('total duration');
%         subplot(2,1,2);plot(duration_false,'ko');grid on
%
%         freq_calls=(best_calls .median_freq(I_auto_match ));
%         freq_false=(best_calls .median_freq(I_false_alarms ));
%
%         figure;
%         subplot(2,1,1);plot(freq_calls,magnitude_calls,'ro');ylim([40 130]);
%         hold off;grid on;xlabel('median frequency');ylabel('magnitude(dB SEL)');
%
%
%         subplot(2,1,2);plot(freq_false,magnitude_false,'ko');ylim([40 130]);
%         grid on;xlabel('median frequency');ylabel('magnitude(dB SEL)');
%
%
%         %Plot contour polynomial fit
%         for III=1:3, %Contour level..
%             figure;
%
%             subplot(2,2,1);
%             plotcc='.x*sdvphox*sdvph';
%             %for III=1:length(best_calls .p),
%             p_calls=best_calls .p{III}(:,I_auto_match );
%             plot(p_calls(end-1,:),p_calls(end-2,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .p{III}(:,I_false_alarms );
%             plot(p_false(end-1,:),20+p_false(end-2,:),[plotcc(2) 'k']);hold on
%             title('scaled slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('slope (Hz/sec)');ylabel('Curvature');
%             grid on
%
%             subplot(2,2,2);
%             p_calls=best_calls .p{III}(:,I_auto_match );
%             plot(p_calls(end,:),p_calls(end-1,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .p{III}(:,I_false_alarms );
%             plot(p_false(end,:),20+p_false(end-1,:),[plotcc(2) 'k']);hold on
%             title('scaled slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('mean F');ylabel('Slope (Hz/sec)');
%             grid on
%
%             subplot(2,2,3);
%
%             p_calls=best_calls .praw{III}(:,I_auto_match );
%             plot(p_calls(end-1,:),p_calls(end-2,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .praw{III}(:,I_false_alarms );
%             plot(p_false(end-1,:),5000+p_false(end-2,:),[plotcc(2) 'k']);hold on
%             title('raw slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('slope (Hz/sec)');ylabel('Curvature');
%             grid on
%             legend('black: airgun','red: calls')
%
%             subplot(2,2,4);
%
%             p_calls=best_calls .praw{III}(:,I_auto_match );
%             plot(p_calls(end,:),p_calls(end-1,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .praw{III}(:,I_false_alarms );
%             plot(p_false(end,:),5000+p_false(end-1,:),[plotcc(2) 'k']);hold on
%             title('raw slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('mean F');ylabel('Slope (Hz/sec)');
%             grid on
%             legend('black: airgun','red: calls')
%         end                                                                                                                                                                           read_tsv.m                                                                                          0000677 0002527 0000000 00000014555 10743503003 013275  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        %%%read_tsv.m%%%%%
%function [ind,localized]=read_tsv(fname,ctmin,ctmax,DASAR_list);
%  Based on documentation in WCoutFormat03.doc by Bill McLennan
%Input:
%   fname:  full filename (including path) of *.tsv file.  Must have tsv
%       extension in string.
%   ctmin: minimum c-time desired
%   ctmax: maximum c-time desired
%   DASAR_list: cell matrix of DASAR names to extract
%Output:
%   ind: row is call number, column is corresponding station in DASAR_list
%       if particular call not found in DASAR, ind.wctype will be 0
%   localized:  localization information about call, see below.
% localized.ctev(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.atev(Icall)=datenum(tabfield(tline,J));J=J+1;
% localized.utmx(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.utmy(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.wctype(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.area(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.axmajor(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.axminor(Icall)=str2num(tabfield(tline,J));J=J+1;
%
% localized.Baxis(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.Ang(Icall)=str2num(tabfield(tline,J));J=J+1;
% localized.Nused(Icall)=str2num(tabfield(tline,J));J=J+5;
% localized.Outcome{Icall}=(tabfield(tline,J));J=J+1;
% localized.comment{Icall}=(tabfield(tline,J));J=J+1;
% localized.OperatorID{Icall}=(tabfield(tline,J));J=J+1;
% localized.DateTimeProc{Icall}=(tabfield(tline,J));J=J+1;
%
% %count(bin)=count(bin)+1;
% for Idd=1:length(DASAR_list),
%     Igood=findstr(DASAR_list{Idd},tline);
%     if ~isempty(Igood),
%         ind.ctime(Icall,Idd)=str2num(tabfield(tline(Igood:end),5));
%         ind.wctype(Icall,Idd)=localized.wctype(Icall);
%         ind.sigdb(Icall,Idd)=str2num(tabfield(tline(Igood:end),6));
%         ind.stndb(Icall,Idd)=str2num(tabfield(tline(Igood:end),7));
%         ind.flo(Icall,Idd)=str2num(tabfield(tline(Igood:end),8));
%         ind.fhi(Icall,Idd)=str2num(tabfield(tline(Igood:end),9));
%         ind.duration(Icall,Idd)=str2num(tabfield(tline(Igood:end),10));
%         %keyboard;
%     end
% end
function [ind,localized]=read_tsv(fname,ctmin,ctmax,DASAR_list);
fid=fopen(fname,'r');
badbins=0;
nline=0;
Icall=0;
if fid ~= -1    % did file open?
    while feof(fid)==0
        tline=fgets(fid);
        nline=nline+1;
        if tline ==-1
            break
        end
        % read computed time of call
        ctime=str2num(tabfield(tline,1));
        % wctype=str2num(tabfield(tline,5));
        %compute bin number
        %bin=ceil((ctime -ctstart)/binsize);
        %%Aaron addition


        if ctime>=ctmin&ctime<=ctmax %& wctype <=8
            Icall=Icall+1;
            if rem(Icall,500)==0, disp(Icall);end
            J=1;
            localized.ctev(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.atev(Icall)=datenum(tabfield(tline,J));J=J+1;
            localized.utmx(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.utmy(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.wctype(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.area(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.axmajor(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.axminor(Icall)=str2num(tabfield(tline,J));J=J+1;

            localized.Baxis(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.Ang(Icall)=str2num(tabfield(tline,J));J=J+1;
            localized.Nused(Icall)=str2num(tabfield(tline,J));J=J+5;
            localized.Outcome{Icall}=(tabfield(tline,J));J=J+1;
            localized.comment{Icall}=(tabfield(tline,J));J=J+1;
            localized.OperatorID{Icall}=(tabfield(tline,J));J=J+1;
            localized.DateTimeProc{Icall}=(tabfield(tline,J));J=J+1;

            %count(bin)=count(bin)+1;
            for Idd=1:length(DASAR_list),
                Igood=findstr(DASAR_list{Idd},tline);
                if ~isempty(Igood),
                    ind.ctime(Icall,Idd)=str2num(tabfield(tline(Igood:end),5));
                    ind.wctype(Icall,Idd)=localized.wctype(Icall);
                    ind.sigdb(Icall,Idd)=str2num(tabfield(tline(Igood:end),6));
                    ind.stndb(Icall,Idd)=str2num(tabfield(tline(Igood:end),7));
                    ind.flo(Icall,Idd)=str2num(tabfield(tline(Igood:end),8));
                    ind.fhi(Icall,Idd)=str2num(tabfield(tline(Igood:end),9));
                    ind.duration(Icall,Idd)=str2num(tabfield(tline(Igood:end),10));

                    %keyboard;
                end
            end

        else
            % display(tline);
            badbins=badbins+1;
            % display(['bin,wctype,badbins='  num2str(bin) blanks(2) num2str(wctype_all(Icall)) blanks(2) num2str(badbins)])
        end
    end   %of file read loop(while)
    fclose(fid);

    %Just in case the last DASAR in DASAR_list had no calls, make sure
    %number of columns in ind same as DASAR_list
    if exist('ind'),
        ncol_needed=length(DASAR_list)-size(ind,2);
        if ncol_needed>0,
            Icall=size(ind.ctime,1);
            ind.ctime=[ind.ctime zeros(Icall,ncol_needed)];
            ind.wctype=[ind.wctype zeros(Icall,ncol_needed)];
            ind.sigdb=[ind.sigdb zeros(Icall,ncol_needed)];
            ind.stndb=[ind.stndb zeros(Icall,ncol_needed)];
            ind.flo=[ind.flo zeros(Icall,ncol_needed)];
            ind.fhi=[ind.fhi zeros(Icall,ncol_needed)];
            ind.duration=[ind.duration zeros(Icall,ncol_needed)];

        end
    end

end

function test
fname='TSV_files_Shell07/wc090107s3.tsv';
[ind,localized]=read_tsv(fname,0,Inf,{'D07s3a','D07s3b','D07s3c','D07s3d','D07s3e','D07s3f','D07s3g'});

%function fs = TabField(s,n)
%
%  Where s is a string containing Tab-delimited fields,
%  this function will return the nth field, with trailing
%  and leading spaces removed.  If n is an array, then this
%  function returns a cell array, in which each cell is
%  the field indicated by the corresponding element in n.

function fs = tabfield(s,n)
%mbcharvector(s);
%mbint(n);
if isempty(n)
    fs = '';
elseif length(n) == 1
    xx = find((s == char(9)));

    if n > (length(xx) + 1)
        fs = [];
    else
        if isempty(xx)
            fs = s;
        else
            xx = [0;xx(:);length(s) + 1];
            fs = s((xx(n) + 1):(xx(n + 1) - 1));
        end
    end
    fs = strtrim(fs);
else
    fs = cell(size(n));
    for ii = 1:prod(size(n))
        fs{ii} = TabField(s,n(ii));
    end
end
return                                                                                                                                                   tabField.m                                                                                          0000777 0002527 0000000 00000001471 10554031174 013177  0                                                                                                    ustar   thode                           wheel                           0000000 0000000                                                                                                                                                                        %function fs = TabField(s,n)
%
%  Where s is a string containing Tab-delimited fields,
%  this function will return the nth field, with trailing
%  and leading spaces removed.  If n is an array, then this
%  function returns a cell array, in which each cell is
%  the field indicated by the corresponding element in n.

function fs = TabField(s,n)
 %mbcharvector(s);
 %mbint(n);
  if isempty(n)
   fs = '';
  elseif length(n) == 1
   xx = find((s == char(9)));

    if n > (length(xx) + 1)
     fs = [];
    else
      if isempty(xx)
       fs = s;
      else
       xx = [0;xx(:);length(s) + 1];
       fs = s((xx(n) + 1):(xx(n + 1) - 1));
      end
    end
   fs = strtrim(fs);
  else
   fs = cell(size(n));
    for ii = 1:prod(size(n))
     fs{ii} = TabField(s,n(ii));
    end
  end
return                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       