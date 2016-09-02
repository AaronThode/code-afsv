%%%master_evaluation.m%%
%% Compare output of automatic detector with Greeneridge manual results
%% Aaron Thode
%%  Sept. 18, 2008
clear
%Icase='Shell08_Site5_allHuber.morph.crosschecked';  %keyword to select time and spatial subset of deployment...
%Icase='Shell08_Site5_PeakNeuralNetTrain.morph.NoNeuralNet';
Icase='Shell08_Site5_PeakBulkRunCore2.morph.Final';

master_setup;

Idate=menu('Which date?',date_str{1});
master_load_automated_data;

master_load_manual_data;

    
%%%Permit an individual location manual or automated index to be selected.  Print out raw
%%%information (time, frequency range, duration) for each selection,
%%%  and plot spectrograms of signal for each location...

disp('USING PARAMETER VALUES STORED WITH AUTOMATIC RESULTS');
param=auto_param;
Ichc=menu('Which type analysis?','Review manual','Review automated','Missed percentage evaluation for multi-station localiations', ...
    'Individual DASAR evaluation, including neural net','Plot manual/automated overlay','Review manual linkage issues');

%call 985 has high SNR%
switch Ichc
    case 1
        Ncalls=length(manual.localized{Idate}.ctev);
        call_indicies=1:Ncalls;

    case {2,5}
        Ncalls=length(auto.locations{Idate});
        call_indicies=1:Ncalls;
      
    case 3
        %tol=[run_options.tol param.feature.tol(1:2)];
        %tol=0.25;
        %         [Ipass,Imiss,Iloc_match]=evaluate_miss_fraction(manual.individual{1},auto.locations{1}, ...
        %             auto.locations_ctime{1},tol,run_options.min_stations,maxtime(1),  ...
        %             run_options.debug_miss_compare);
        tic
        [Ngood_all,Nstation_match,Nlongest_link,Nsplit_links]=evaluate_linking(param,manual.individual{1},auto.locations{1}, ...
            auto.locations_ctime{1},auto.stations{1},auto.raw_stations{1},date_str_new{1}{1},goodDASAR{1},run_options);
        toc
        
        
        %%Summarize statistics
        disp(sprintf('Out of %i manual detections, %i had more than 2 detections in stations',length(Nstation_match),length(find(Nstation_match>=2))));
        disp(sprintf('Out of %i manual detections, %i had more than 2 detections in locations',length(Nstation_match),length(find(Ngood_all>=2))));
        disp(sprintf('Out of %i manual detections, %i had more than 2 detections in a single linkage',length(Nstation_match),length(find(Nlongest_link>=2))));
        
        if run_options.auto_success_only==1
            area_vector=-4:8;
            [N,X]=hist(log10(auto.error_area/1e6),area_vector);
            [Nm,Xm]=hist(log10(manual.localized{1}.area/1e6),area_vector);
            bar(X,[N' Nm'],'grouped')
            xlabel('log area error, km^2');grid
            legend('automated','manual')
        end


        Ichc2=menu('Examine what?','Missed calls','matched calls','False alarms');
        if Ichc2==1,
            Ncalls=length(Imiss.count);
            call_indicies=Imiss.count;
        else
            Ncalls=length(Ipass.count);
            call_indicies=Ipass.count;
        end


    case 4  %Evalute individual DASARs

        %Compute maximum ctime to permit..

        [whale,other,Ntrue]=evaluate_individual_stations(param,manual.individual{1},auto.stations{1}, ...
            auto.raw_stations{1},date_str_new{1}{1},goodDASAR{1},run_options);

        Ichc=menu('What feature analysis desired?','Fisher discriminant and ROC curve','Hand-picked 2D feature plots','Neural Net','Missed call analysis','None');

        switch Ichc
            case 1
                compute_fisher_discriminant(whale,other,param,Ntrue);


            case 2

                %feature_names=auto.stations{1}(1).feature_name;
                feature_names=fieldnames(whale{1});
                [Ichcx,Ichcy]=plot_2D_features(whale,other,feature_names);
                coords=ginput(1);
                disp(sprintf('Coordinates chosen: %s: %6.2f and %s: %6.2f',feature_names{Ichcx},coords(1),feature_names{Ichcy},coords(2)));
                plot_feature_filtered_spectrogram_and_morph(other,goodFile{1},param,feature_names,Ichcx,Ichcy,coords,'False');
                plot_feature_filtered_spectrogram_and_morph(whale,goodFile{1},param,feature_names,Ichcx,Ichcy,coords,'Whale');

            case 3  %Neural net
                %construct input matrix
                miss_fraction=input('Enter miss fraction (default 10%):');
                if isempty(miss_fraction)
                    miss_fraction=0.1;
                end

                nweights=input('Enter number of hidden cells (default 3):');
                if isempty(nweights)
                    nweights=3;
                end

                run_options.pca=0; %Principal component analysis on inputs?
                run_options.debug_plot=1;

                %Remove assymetric structure from whale and other
                for K=1:length(whale),
                    whale{K}=rmfield(whale{K},{'miss','Image','param','equalization'});
                    other{K}=rmfield(other{K},{'Image','param','equalization'});
                end

                run_options.skip_feature=param.net.skip_feature;

                [pattern,target,historry,net,threshold,preprocess_info,preprocess_info_pca,feature_names]=build_nnetwork(whale,other,miss_fraction,nweights,run_options);

                yes=menu('Save network?','Yes','No');
                if yes==1,
                    if ~isfield(param,'label')
                        param.label='Shell08_nnet';
                    end
                    strr='';
                    if run_options.pca==1,
                        strr='_pca';
                    end
                    network_save_name=sprintf('Nnetwork_%s_%ihidden%s_%s',param.label,size(net.IW{1},1),strr,date_str_new{1});
                    mydir=pwd;
                    cd(automated_dir);
                    cd ..
                    save(network_save_name,'net','param','threshold','preprocess_info','preprocess_info_pca','nweights','miss_fraction','feature_names');
                    cd(mydir);

                    keyboard;
                end
            case 4  %Plot statistics of missed calls
                plot_missed_call_stats_and_spectrogram(whale,manual.localized{1}, param,goodFile{1},goodName{1},date_range(1,:),auto.stations{1});

            case 6  %Use an external csv file to learn what manual calls did not acheive linking stage..
                test_linkages;
                

        end

end

Icall=input(sprintf('Enter desired index [1-%i]:',Ncalls));
while Icall>0,
    switch Ichc
        case 5
            close all
            %Icall_manual=input('Enter manual index:');
            %Icall_manual=14;  Icall_auto=9;
            %Icall_manual=3;Icall_auto=3;  %OK example...
            Icall_manual=24;Icall_auto=15;  %OK example...
            ctimes_out=display_manual_location_data(manual.localized{1},manual.individual{1},Icall_manual,param,goodFile{1});
            Ikeep=find(manual.individual{1}.ctime(Icall_manual,:)>0);
            figure(10);
            subplot(1,2,1);
            plot_location(DASAR_coords{1},manual.individual{1}.bearing(Icall_manual,:),Ikeep, ...
                [manual.localized{1}.utmx(Icall_manual,:) manual.localized{1}.utmy(Icall_manual,:)] , ...
                manual.localized{1}.axmajor(Icall_manual,:),manual.localized{1}.axminor(Icall_manual,:),manual.localized{1}.Ang(Icall_manual,:));
           
            figure(1);
            %Icall_auto=input('Enter auto index:');
            %Icall_auto=9;
            display_automated_crosslink_data(auto.locations{1},Icall_auto,param,run_options,goodFile{1},ctimes_out);
            figure(10);
            subplot(1,2,2);
            pos=auto.locations{1}{call_indicies(Icall_auto)}.position;
            plot_location(DASAR_coords{1},auto.locations{1}{call_indicies(Icall_auto)}.bearing,pos.Ikeep,pos.location,pos.major,pos.minor,pos.ellipse_ang);

            
        case {1,3},

            display(sprintf('Analyzing manual index %i',call_indicies(Icall)));
            display_manual_location_data(manual.localized{1},manual.individual{1},call_indicies(Icall),param,goodFile{1});
           
            Ikeep=find(manual.individual{1}.ctime(Icall,:)>0);
            figure;
            plot_location(DASAR_coords{1},manual.individual{1}.bearing(Icall,:),Ikeep, ...
                [manual.localized{1}.utmx(Icall,:) manual.localized{1}.utmy(Icall,:)] , ...
                manual.localized{1}.axmajor(Icall,:),manual.localized{1}.axminor(Icall,:),manual.localized{1}.Ang(Icall,:));
           

            pause;%close all;

         %  identify_location_failure(Site_vector,manual.localized{1},manual.individual{1}, ...
         %      auto.locations{1},auto.locations_ctime{1},auto.stations{1},auto.raw_stations{1}, ...
         %      call_indicies(Icall),param,run_options,goodFile{1},goodName{1},automated_dir);

        case 2
            display_automated_crosslink_data(auto.locations{1},call_indicies(Icall),param,run_options,goodFile{1});
            pos=auto.locations{1}{call_indicies(Icall)}.position;
            figure;
            plot_location(DASAR_coords{1},auto.locations{1}{call_indicies(Icall)}.bearing,pos.Ikeep,pos.location,pos.major,pos.minor,pos.ellipse_ang);
            dt_err=NaN*ones(length(pos.Ikeep),1);
            dt_true=dt_err;
            dt_est=dt_true;
            Nunits=size(DASAR_coords{1},1);
            ranges=sqrt(sum((ones(Nunits,1)*pos.location-DASAR_coords{1}).^2,2));
            Ianchor=find(auto.locations{1}{call_indicies(Icall)}.dt==0);
            dt_est=(ranges(Ianchor)-ranges(pos.Ikeep))./1490;
            dt_true=auto.locations{1}{call_indicies(Icall)}.dt(pos.Ikeep);
            %cos_match=dt_est'*dt_true./(norm(dt_est).*norm(dt_true));
            dt_err=abs(dt_true-dt_est);
            disp(sprintf('Measured dt: %s\nModeled dt:  %s\nAbs diff: %s',mat2str(dt_true,4),mat2str(dt_est,4), ...
                mat2str(dt_err,4)));
            %dt_err(Ianchor)=0;

    end

    Icall2=input(sprintf('Enter desired index [1-%i], return to increment by 1, or -1 to quit:',Ncalls));
    if isempty(Icall2),
        Icall=Icall+1;
    else
        Icall=Icall2;
    end
    if Icall~=-1,
        close all
    end
end



