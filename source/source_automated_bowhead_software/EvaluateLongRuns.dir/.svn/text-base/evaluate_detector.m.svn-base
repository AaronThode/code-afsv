%%%evalute_detector.m%%
%% Compare output of automatic detector with Greeneridge manual results
%% Aaron Thode
%%  Sept. 18, 2008


clear all;close all
path(path,'../CommonScripts.dir');
[rawdatadir,Icase,detectiondir,manualdir]=load_pathnames(Icase);

year_str='07';

Icase='Site_4_short.morph.crosschecked';  %Icase

mydir=pwd;
Idebug.val=0;  %If one, show deubug output, and load debug_ctimes to learn when "failure" occures
run_options.load_false_reviews=0;  %If one, load files that contain reviews of false detections that were actually calls.
%The last makes a list of random ctimes to compare performance.
run_options.extract_high_SNR_only=0;  %If one, use a manual detection only if the highest SNR across the site
% for that DASAR
run_options.ICI_strip=0; %Try to strip airgun repetitive patterns
run_options.feature_strip=0;  %Use morphological features to filter detections after ICI stripping...
run_options.duplicate_time_tol=1; %Tolerance time for matching manual and automated detection
params_chc='June6_morph_nofilt';
run_options.tol=3;

%%for ICI analysis
run_options.ICI_strip.Ndet=5;
run_options.ICI_strip.Nmiss=1;
run_options.ICI_strip.ICItol=0.5;

%[detection_stage,date_range,DASAR_list,autoCase,DASAR_list_total]=TOC_data(Icase);
[Site_vector,date_str,DASAR_str,keyword,detection_stage]=TOC_runs(Icase);
date_range=expand_date_range(date_str);
for Idate=1:size(date_range,1),
    [goodDASAR{Idate},goodFile{Idate},goodName{Idate}]=find_DASAR_dates(date_range(Idate,:),Site_vector,DASAR_str,rawdatadir);
end
%if run_options.extract_high_SNR_only==1,
[manual,manual_location,SNR,false_notes,uncertain_notes]=load_manual_results(manualdir,goodDASAR,date_range,run_options);
%else
% [manual,SNR,dates,false_notes,uncertain_notes]=load_manual_results(DASAR_list,date_range,load_false_reviews);
%    [manual,SNR,dates,false_notes,uncertain_notes]=load_manual_results(goodName{Idate},date_range,run_options);
%end
%dbstop if error
keyboard;

for Istage=1:length(detection_stage),
    disp(sprintf('**Testing detection stage: %s',detection_stage{Istage}));
    %[auto_ctimes,best_calls,airgun_calls,param]=load_automated_results(Icase,DASAR_list,detection_stage{Istage},dates);
    %[auto_ctimes,best_calls,airgun_calls,param]=load_automated_results(detectiondir,goodDASAR,date_range,detection_stage{Istage});
    [auto_ctimes,best_calls,airgun_calls,param]=load_automated_results(keyword,detectiondir, ...
        goodName,date_range,detection_stage{Istage});

    %Ndet=param.interval_remove.Ndet;
    %ICItol=param.interval_remove.ICItol;
    %random_ctimes=rand(size(auto_ctimes));
    %Test whether ICI can be identified

    total.missed=0;
    total.false=0;
    total.correct=0;
    for Iloc=1:length(DASAR_list),  %DASAR loop
        for Iday=1:length(date_range),  %Day loop

            %%Clean up manual ctimes
            %Igood=find(manual{Iday}.wctype(:,Iloc)>0);
            %manual_ctime=manual{Iday}.ctime(Igood,Iloc);

            disp(sprintf('Day: %s DASAR: %s',dates{Iday},DASAR_list{Iloc}));
            auto_ctime=auto_ctimes{Iday,Iloc};

            %Study if stripping tol works
            %auto_ctimes_tol=crude_decimate_uneven(auto_ctime,auto_ctime,tol,'median');
            %disp(sprintf('Original detection length: %i, with median filter applied: %i',length(auto_ctime),length(auto_ctimes_tol)));
            %auto_ctime=auto_ctimes_tol;

            %%Load manual detection data, and add SNR information if manual
            %%audits (false_call_reviews) have been conducted
            manual_ctime=manual{Iday,Iloc};
            if load_false_reviews==1,
                Nmanual=length(manual_ctime);

                for I=1:Nmanual,
                    if SNR{Iday,Iloc}(I)<0,
                        [delt,Iwant]=min(abs(manual_ctime(I)-best_calls{Iday,Iloc}.ctime));
                        if delt<tol,
                            SNR{Iday,Iloc}(I)=best_calls{Iday,Iloc}.features{Iwant}.SNR;
                        end

                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%Optional ICI stripping%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%Optional feature stripping%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if feature_strip==1,
                disp('Stripping features from Morphological');
                %keyboard;
                param=TOC_params(params_chc)

                disp('filtering AND rules');
                Iand=find(param.feature.operator>0);
                feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features,param.feature.feature_name(Iand), ...
                    param.feature.feature_index(Iand),param.feature.Nshapes);

                [best_filtered, auto_ctime,Imatch]=filter_feature_vector(feature_vector,param.feature.feature_name(Iand), ...
                    auto_ctime, param.feature.optvec(Iand,:),'all','and',1);


                disp('filtering OR rules');
                Ior=find(param.feature.operator==0);
                if ~isempty(Ior),
                    feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features(Imatch),param.feature.feature_name(Ior), ...
                        param.feature.feature_index(Ior),param.feature.Nshapes);

                    [best_filtered, auto_ctime,Imatchor]=filter_feature_vector(feature_vector,param.feature.feature_name(Ior), ...
                        auto_ctime, param.feature.optvec(Ior,:),'all','or',1);

                    Imatch=Imatch(Imatchor);
                end

                best_calls{Iday,Iloc}.features=best_calls{Iday,Iloc}.features(Imatch);
                best_calls{Iday,Iloc}.ctime=best_calls{Iday,Iloc}.ctime(Imatch);
                best_calls{Iday,Iloc}.BWfinal=best_calls{Iday,Iloc}.BWfinal(Imatch);


                %feature_vector=extract_feature_vector(best_calls{Iday,Iloc}.features,param.feature.feature_name, ...
                %  param.feature.feature_index,param.feature.Nshapes);

            end

            %%Remove detections that match 'uncertain' times from
            %%consideration
            if 1==0,
                [junk,junk2,junk3,Ipossible]=find_similar_elements(uncertain_notes{Iday}.ctimes,auto_ctime,tol);
                disp(sprintf('%i detections associated with uncertain classifications, removed from  consideration',length(junk3)));
                %auto_ctime=auto_ctime(Ipossible);  %Ifalse are times that do not correcspond to uncertain times
                [Ifalse,Ix_match,Iy_match,Iy_nomatch,mindiff]=find_similar_elements(auto_ctime(Ipossible),manual_ctime,tol);
                Ix_nomatch=Ipossible(Ifalse);  %indicies of best_call that are false.
                Ix_match=Ipossible(Ix_match);
                %Iy_match=Ipossible(Iy_match);
                %Iy_nomatch=Ipossible(Iy_nomatch);
                %keyboard;
            else
                [IfalseNuncertain,Ix_match,Iy_match,Iy_nomatch]=find_similar_elements(auto_ctime,manual_ctime,tol);
                false_and_uncertain=auto_ctime(IfalseNuncertain);

                disp('');
                % disp(sprintf('Cumulative missed: %i Cumulative correct: %i',total.missed,total.correct));
                %Strip false alarms that have been classified as 'uncertain')
                false_alarms=false_and_uncertain;

                if ~isempty(uncertain_notes{1}.ctimes),
                    Ixmatch=[0];
                    while length(Ixmatch)>0,
                        [Ix_nomatch,Ixmatch]=find_similar_elements(false_alarms,uncertain_notes{Iday}.ctimes,tol);
                        %keyboard;
                        disp(sprintf('Uncertain signals: %i, Certain false : %i, total false calls was  %i \n',length(Ixmatch),length(Ix_nomatch),length(false_alarms)));
                        false_alarms=false_alarms(Ix_nomatch);

                    end
                    [junk,junk2,Ifalse]=intersect(false_alarms,false_and_uncertain);
                    Ix_nomatch=IfalseNuncertain(Ifalse);

                else
                    Ix_nomatch=IfalseNuncertain;
                end
            end
            %%%Compare SNR estimates
            try,
                for II=1:length(Ix_match),

                    SNR_me(II)=best_calls{Iday,Iloc}.features{Ix_match(II)}.SNR;
                end
                SNR_tsv=SNR{Iday,Iloc}(Iy_match);
                %disp('SNR test');
            end

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
            %times_false_mismatch{Iday,Iloc}=Ix_mismatch_time;
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


