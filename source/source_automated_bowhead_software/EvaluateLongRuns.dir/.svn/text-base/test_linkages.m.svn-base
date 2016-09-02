%Icase='Shell08_Site5_allHuber.morph.crosschecked';  %keyword to select time and spatial subset of deployment...
%Icase='Shell08_Site5_PeakNeuralNetTrain.morph.NoNeuralNet';
clear ;close all
Icase='Shell08_Site5_PeakBulkRunCore2.morph.Final';

master_setup;

total.manual_alltime=0;
total.manual=0;
total.miss=0;
total.toofew=0;
for Idate=1:length(date_str{1})

    master_load_automated_data;

    master_load_manual_data;


    %%%Permit an individual location manual or automated index to be selected.  Print out raw
    %%%information (time, frequency range, duration) for each selection,
    %%%  and plot spectrograms of signal for each location...

    disp(sprintf('\n\n\nUSING PARAMETER VALUES STORED WITH AUTOMATIC RESULTS'));
    param=auto_param;



    switch date_str{1}{Idate}
        case '20080828'
            WEST_locs=csvread('../LinkingEvaluation.dir/Site5_manual/out_0828S508_m.csv',1,0);
        case '20080821'
            WEST_locs=csvread('../LinkingEvaluation.dir/Site5_manual/out_0821S508_m.csv',1,0);
        case '20080906'
            WEST_locs=csvread('../LinkingEvaluation.dir/Site5_manual/out_0906S508_m.csv',1,0);
        case '20080913'
            WEST_locs=csvread('../LinkingEvaluation.dir/Site5_manual/out_0913S508_m.csv',1,0);
        case '20080921'
            WEST_locs=csvread('../LinkingEvaluation.dir/Site5_manual/out_0921S508_m.csv',1,0);
        case '20080929'
            WEST_locs=csvread('../LinkingEvaluation.dir/Site5_manual/out_0929S508_m.csv',1,0);
    end

    %%Locate manuals not tracked
    Imiss=find(WEST_locs(:,2)==0&WEST_locs(:,3)<maxtime);
    missed_times=WEST_locs(Imiss,2);
    Nmanual=length(find(WEST_locs(:,3)<maxtime));
    total.manual=total.manual+Nmanual;
    total.manual_alltime=total.manual_alltime+length(WEST_locs(:,2));
    disp(sprintf('There are %i manual calls in WEST file, %i calls so far ',length(WEST_locs(:,2)),total.manual_alltime));
  
    disp(sprintf('There are %i manual calls in WEST file, %i calls in my TSV processed file less than %s..',Nmanual,length(manual.localized{1}.ctev),ctime2str(maxtime)));
    disp(sprintf('There are %i missing calls out of %i total manual calls',length(Imiss),Nmanual));

    Nauto_maybe=0;
    for JJ=1:length(missed_times)

        %Locate appropriate index in manual results
        Iman=find(WEST_locs(Imiss(JJ),3)==manual.localized{1}.ctev);

        if length(Iman)>1
            disp(sprintf('%s duplicated, index %i,',ctime2str(WEST_locs(Imiss(JJ),3)),JJ));
        end
        Iloc=find(manual.localized{1}.utmx(Iman)==WEST_locs(Imiss(JJ),4));
        Iman=Iman(Iloc);
        if isempty(Iman)
            disp(sprintf('%s not found, index %i,',ctime2str(WEST_locs(Imiss(JJ),3)),JJ));
        end

        %OK, were individual DASARS detected?
        ctime_stations=zeros(1,7);
        Nmanual_hits=length(find(manual.individual{1}.ctime(Iman,:)>0));
        disp('');
        %disp(sprintf('\nThere are %i manual DASARs in this detection',Nmanual_hits));
        for KK=1:7
            ctime_auto=auto.stations{1}(KK).ctime_min;
            duration_auto=auto.stations{1}(KK).Totalduration;
            ctime_man=manual.individual{1}.ctime(Iman,KK);
            duration_man=manual.individual{1}.duration(Iman,KK);

            [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match,Iy_nomatch_diff,Iy_match_diff]=find_similar_elements_ovlap(ctime_auto,duration_auto,ctime_man,duration_man,0.25);

            if ~isempty(Ix_match)
                ctime_stations(KK)=ctime_auto(Ix_match);
                %disp(sprintf('Index %i, DASAR %i, manual time: %s+%6.2f sec, auto time %s+%6.2f, overlap %10.8f',JJ,KK, ...
                %    ctime2str(ctime_man),duration_man,ctime2str(ctime_auto(Ix_match)),duration_auto(Ix_match), Iy_match_diff));
            end
        end
        %disp(sprintf('\n\n\n'))

        if length(find(ctime_stations>0))<2
            Nauto_maybe=Nauto_maybe+1;
        end
    end

    disp(sprintf('%s; Out of %i missed calls, %i arise from insufficient DASARs',date_str{1}{Idate},length(Imiss),Nauto_maybe));
    total.miss=total.miss+length(Imiss);
total.toofew=total.toofew+Nauto_maybe;
end

disp(sprintf('Complete site: Out %i manual calls %i missed calls, %i arise from insufficient DASARs',total.manual,total.miss,total.toofew));