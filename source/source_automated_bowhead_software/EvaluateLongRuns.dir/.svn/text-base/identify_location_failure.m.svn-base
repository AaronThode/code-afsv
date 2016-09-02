
function param_out=identify_location_failure(Isite,manlocal,manind, autolocal,auto_ctime,auto_station,auto_raw_stations, ...
    Icall,param,run_options,filename,auto_template_name,automated_dir)

%persistent auto

%%Unpack manual information
%strr_org='abcdefghijklmnop';
strr_org=manind.DASARstr;
Nstation=size(manind.duration,2);
disp(' ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%% identifty_location_failure%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%%What manual stations do not have call
Istation=manind.anchor_station(Icall);
Inot=find(manind.ctime(Icall,:)==0);  %Stations where call is not registered.
manind.ctime(Icall,Inot)=-Inf;

%What manual stations detected call?
Igood=find(manind.ctime(Icall,:)>0);
dt=-manind.ctime(Icall,Istation)+manind.ctime(Icall,Igood);

%Extract manual parameters
duration=(manind.duration(Icall,Igood));
flo=manind.flo(Icall,Igood);
fhi=manind.fhi(Icall,Igood);
SNR=manind.stndb(Icall,Igood);
level=manind.sigdb(Icall,Igood);
wctype=manind.wctype(Icall,Igood);
ctime=manind.ctime(Icall,Igood);
tabs=datenum(1970,1,1,0,0,ctime);
strr=strr_org(Igood);
disp(sprintf('Analyze  manually localized call %i at %s',Icall,datestr(datenum(1970,1,1,0,0,manlocal.ctev(Icall)))));
disp(sprintf('Call type %6.2f, stations used %s, comment %s', manlocal.wctype(Icall),strr,manlocal.comment{Icall}));

Nfig=0;Iplot=0;myfig=0;Nrows=3;

%%Review Each station
%for I=1:length(Igood),  %For every station used in manual result...
for I=1:length(Igood),
    disp(sprintf('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv'));

    start_time=[];
    clear y tlen start_ctime 
    [ctimes,durations,Index,loc]=first_stage_match;

    final(I).Icand=Index.final;
    %Does station have raw material?
    disp(sprintf('Station %s call ctime is %16.12f or %s',strr(I), ctime(I),ctime2str(ctime(I))));
    disp(sprintf('Station %s has %i raw, %i interval, %i morph detections and %i localizations within %6.2f sec of desired call.', ...
        strr(I),length(Index.raw),length(Index.int),length(final(I).Icand),length(loc.Icand),run_options.tol));


    %%Review features of detected calls...

    if ~isempty(final(I).Icand),  %If a detection is present in final results
        names=auto_station(Igood(I)).param.feature_name;
        for Ifea=1:length(names);
            morph_feature.(names{Ifea})=auto_station(Igood(I)).feature.(names{Ifea})(:,final(I).Icand);
        end
        final(I).SEL=auto_station(Igood(I)).SEL(final(I).Icand);
        final(I).flo=auto_station(Igood(I)).Totalfmin(final(I).Icand);
        final(I).fhi=auto_station(Igood(I)).Totalfmax(final(I).Icand);
        final(I).duration=auto_station(Igood(I)).Totalduration(final(I).Icand);


        disp(sprintf(' %s: Manual detection %i has flo %6.2f fhi %6.2f duration %6.2f SNR: %6.2f',strr(I),Icall, flo(I),fhi(I),duration(I),level(I)));
        disp(sprintf('     Morph  detection %i has flo %6.2f fhi %6.2f duration %6.2f SNR: %6.2f \n', [final(I).Icand; final(I).flo; final(I).fhi; final(I).duration; final(I).SEL]));
        if ~isempty(loc.Icand),
            disp(sprintf('     Loc index %i flo %6.2f fhi %6.2f duration %6.2f dt %6.2f \n', ...
                [ loc.Icand; loc.flo; loc.fhi; loc.duration; loc.dt]));
            for JJJ=1:length(loc.present),
                disp(sprintf('Loc index %i has ''%s'' stations indicies %s with ''%s'' anchoring.\n', ...
                    loc.Icand(JJJ),loc.present{JJJ},mat2str(loc.indicies{JJJ}'),loc.anchor(JJJ)));
            end
        else
            disp(sprintf('!!!!! Manual call not in localized data, closest index is %i ofset %6.2f s!!!!',loc.closest_index,loc.closest_time));
        end


        %Download raw time series for raw morph check,

        for JJ=1:length(final(I).Icand)
            [y{JJ},start_ctime(JJ),tlen(JJ),head]=extract_signal_from_station(auto_station(Igood(I)),final(I).Icand(JJ),filename{Igood(I)},'debug');

        end
        final(I).cstart= ctimes.final-start_ctime;
        tit_lab='final';
    elseif ~isempty(Index.morph)  %if detection was present in morph detector but failed neural network.
        for JJ=1:length(durations.morph)
            %tlen(JJ)=4*param.energy.bufferTime;
            tlen(JJ)=durations.morph(JJ)+2*param.energy.bufferTime;
            start_ctime(JJ)=ctimes.morph(JJ)-param.energy.bufferTime;
            [y{JJ},t,head]=readfile(filename{Igood(I)},start_ctime(JJ),tlen(JJ),1,'ctime','calibrate');
        end
        tit_lab='morph';
    elseif ~isempty(Index.int)  %if detection was present in interval detector.
        for JJ=1:length(Index.int)
          
            tlen(JJ)=durations.int(JJ)+2*param.energy.bufferTime;
            start_ctime(JJ)=ctimes.int(JJ)-param.energy.bufferTime;
            
            [y{JJ},t,head]=readfile(filename{Igood(I)},start_ctime(JJ),tlen(JJ),1,'ctime','calibrate');
        end
        tit_lab='interval';
    elseif ~isempty(Index.raw) %in energy detector only...
        for JJ=1:length(Index.raw)
            tlen(JJ)=durations.raw(JJ)+2*param.energy.bufferTime;
            start_ctime(JJ)=ctimes.raw(JJ)-param.energy.bufferTime;
             [y{JJ},t,head]=readfile(filename{Igood(I)},start_ctime(JJ),tlen(JJ),1,'ctime','calibrate');
        end
        tit_lab='raw';
    else  %never detected
        tlen=4*param.energy.bufferTime;
        start_ctime=ctime(I)-2*param.energy.bufferTime;
        [y{1},t,head]=readfile(filename{Igood(I)},start_ctime,tlen,1,'ctime','calibrate');
        tit_lab='none';
    end  %isempty Icand

    %%Actually process raw data with morphological processor to view
    %%performance...
    Iorg=I;
    for JJ=1:length(y),
        disp(sprintf('%i of %i possibilities for %s being examined',JJ,length(y),tit_lab));
        man_ctime=ctime(I)-start_ctime(JJ);
        if ~iscell(y),
            keyboard;
        end

        [Iplot,Nfig,myfig]=plot_spectrogram_and_overlay(Iplot,Nfig,y{JJ},param,man_ctime,start_ctime(JJ), myfig,final(I),duration(I),flo(I),fhi(I), Nrows,strr(I));

        
        disp(sprintf('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n'));
        switch tit_lab
            case 'final'
            case 'morph'
                title(sprintf('moprh pass but no network pass: center time of %s',ctime2str(start_ctime(JJ))));

            case 'interval'
                title(sprintf('interval pass but no morph pass: center time of %s',ctime2str(start_ctime(JJ))));
            case 'raw'
                title(sprintf('raw pass but no interval pass: center time of %s',ctime2str(start_ctime(JJ))));
            case 'none'
                title(sprintf('no pass: center time of %s',ctime2str(start_ctime(JJ))));
        end
        hold off

        if ~isempty(final(I).Icand)
            station_index=final(I).Icand(JJ);
        else
            station_index=[];
        end
        
        [param,TT,FF,index_offset_t,stored_image,labeled_image]=adjust_morph(auto_station,filename,param,Igood(I), station_index,param.morph.eq_time,0,start_ctime(JJ),tlen(JJ));
 
        figure(myfig);
        subplot(Nrows,2,2*Iplot);
        plot_morph_image(TT,FF,station_index,index_offset_t,labeled_image);
        
        %if isempty(final(I).Icand),
            xlimm=xlim;
            xlimm(1)=-param.morph.eq_time;
            xlim(xlimm);
        %end
         Iplot=Iplot-1;
        hold off
        if strcmp(tit_lab,'none')|isempty(Index.raw)
            yes=menu('Energy detector review?: ','Yes','No');

            if yes==1,
                tbuffer=30;
                tstart=datenum(1970,1,1,0,0,head.ctbc);
                twant=datenum(1970,1,1,0,0,start_ctime(JJ)-tbuffer);
                tlenn=tbuffer+10;
                tview=tbuffer+2*param.energy.bufferTime+[-5 5];
                Energy_detector_review(param,filename{Igood(I)},tstart,twant,tlenn,tview);
                yes=menu('Energy detector review?: ','Yes','No');

            end
        end
        
         
    end
    Iplot=Iplot+1;


end  %plots of manual review...


%%%Plot spectrograms of stations with no data
plot_no_morph(auto_station,Igood,Nstation,param,Nrows,strr);


%%%%Test cross-channel match...

yes=menu('Test cross-channel match?','Yes','No');
if yes==1,
    test_cross_channel_match;
else
    %keyboard;
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%');

%%%%%%%%%%%%%%%%%test_cross_channel_match.m%%%%%%%%%%%
    function test_cross_channel_match
        
        Isite=min(findstr(filename{1},'Site_'))+6;
        if isempty(Isite)
           Isite=findstr(filename{1},'/GSI_08')+9; 
        end
        Isite=str2num(filename{1}(Isite));
        
        for K=1:length(Igood)
            SEL=final(K).SEL;
            if ~isempty(SEL)
                [local_max(K),loc_ind(K)]=max(SEL);
            else
                local_max(K)=-1;loc_ind(K)=-1;
            end
        end
        [junk,Imax]=max(local_max);

        debug_params.J_anchor=input('Enter anchor station:');
        if isempty(debug_params.J_anchor),
            debug_params.J_anchor=Igood(Imax);
        end
        
        debug_params.J_anchor_start=input('Enter anchor index:');
        if isempty(debug_params.J_anchor_start)
            debug_params.J_anchor_start=final(Imax).Icand(loc_ind(Imax));
        else
            local_max(Imax)=auto_station(debug_params.J_anchor).SEL(debug_params.J_anchor_start);
        end


        debug_params.cross_channel=2;
        debug_params.Nfft=param.Nfft;
        debug_params.ovlap=param.ovlap;
        debug_params.Fs=param.Fs;
        debug_params.morph=param.morph;
        debug_params.merge=param.merge;

        disp(sprintf('Will use station %i, index %i for anchor, with SEL %6.2f',debug_params.J_anchor,debug_params.J_anchor_start,local_max(Imax)));
        tmp1.tol=param.feature.tol;
        tmp1.weight=param.feature.weight;
        [tmp,change_flag]=alter_parameters(tmp1);
        param.feature.tol=tmp.tol;
        param.feature.weight=tmp.weight;


        while debug_params.cross_channel==2,
            %cross_channel_match(Isite,station,feature_params,goodFile,run_options,debug_params),

            
            [locations, locations_ctime]=cross_channel_match(Isite,auto_station,param.feature,filename,run_options,debug_params);
            if ~isempty(locations),
                disp('');disp('!!!!!!!!!!!!!');disp('Succesful localization....');
                namess=fieldnames(locations{1});
                for JJ=1:(length(namess)-1),  %skip features...
                    disp(sprintf('%s: %s\n',namess{JJ},mat2str(locations{1}.(namess{JJ})',6)));
                end

                namess=fieldnames(locations{1}.feature);
                for JJ=1:(length(namess)),  %skip features...
                    disp(sprintf('%s: %s\n',namess{JJ},mat2str([locations{1}.feature.(namess{JJ})],6)));
                end

                locations=compute_bearings(locations,filename,param,'sel_ratio_FFT',0);

                %%Compute final location..
                run_options.plot_locations=1;
                run_options.localization_alg='Andrews';
                locations=compute_position(locations,filename,param,filename{1},Isite,run_options);

            end
            debug_params.cross_channel=input('Enter ''2'' to test a different anchor station:');
            if debug_params.cross_channel==2,
                debug_params.J_anchor=input('Enter anchor station number:');
                debug_params.J_anchor_start=input('Enter start index number:');

                tmp1.tol=param.feature.tol;
                tmp1.weight=param.feature.weight;
                [tmp,change_flag]=alter_parameters(tmp1);
                param.feature.tol=tmp.tol;
                param.feature.weight=tmp.weight;
            end
        end
    end


%%%%%%%%%%%%%%%%%%%first_stage_match.m%%%%%%%%%%%%%%%%%%%%%%%
    function [ctimes,durations,Index,loc]=first_stage_match


        ctimes.raw=auto_raw_stations(Igood(I)).raw_detections.ctime;
        ctimes.int=auto_raw_stations(Igood(I)).interval_detections.ctime;
        ctimes.morph=auto_raw_stations(Igood(I)).morph_detections.ctime;
        ctimes.final=auto_station(Igood(I)).ctime_min;
        
        durations.raw=auto_raw_stations(Igood(I)).raw_detections.duration;
        durations.int=auto_raw_stations(Igood(I)).interval_detections.duration;
        durations.morph=auto_raw_stations(Igood(I)).morph_detections.duration;
        durations.final=auto_station(Igood(I)).Totalduration;

        ovlap_tol=0.25;
        [Ifalse,Imiss,Index.raw]=find_similar_elements_ovlap(ctimes.raw,durations.raw, ctime(I),duration(I),ovlap_tol);
        [Ifalse,Imiss,Index.int]=find_similar_elements_ovlap(ctimes.int,durations.int, ctime(I),duration(I),ovlap_tol);
        [Ifalse,Imiss,Index.morph]=find_similar_elements_ovlap(ctimes.morph,durations.morph, ctime(I),duration(I),ovlap_tol);
        [Ifalse,Imiss,Index.final]=find_similar_elements_ovlap(ctimes.final,durations.final, ctime(I),duration(I),ovlap_tol);
          
        ctimes.raw=ctimes.raw(Index.raw);
        ctimes.int=ctimes.int(Index.int);
        ctimes.morph=ctimes.morph(Index.morph);
        ctimes.final=ctimes.final(Index.final);
        
        durations.raw=durations.raw(Index.raw);
        durations.int=durations.int(Index.int);
        durations.morph=durations.morph(Index.morph);
        durations.final=durations.final(Index.final);
        
        %Index.raw=find(abs(ctimes.raw-ctime(I))<=run_options.tol);
        %Index.int=find(abs(ctimes.int-ctime(I))<=run_options.tol);
        %Index.morph=find(abs(ctimes.morph-ctime(I))<=run_options.tol);
        %Index.final=find(abs(ctimes.final-ctime(I))<=run_options.tol);


        %Searching auto location is more complex..
        loc.Icand=[];
        loc.flo=[];
        loc.fhi=[];
        loc.duration=[];
        loc.dt=[];
        loc.anchor=[];
        loc.present=[];
        loc.indicies=[];
        loc.closest_time=Inf;
        loc.closest_index=-1;
        
        for J=1:length(autolocal),

            [Ifalse,Imiss,Iauto_match,Imatch,junk,match_ovlap]=find_similar_elements_ovlap(autolocal{J}.ctime_min(Igood(I)), ...
                autolocal{J}.Totalduration(Igood(I)), ctime(I),duration(I),ovlap_tol);
            test_time=abs(autolocal{J}.ctime_min(Igood(I))-ctime(I));

            if ~isempty(Iauto_match),
                loc.Icand=[loc.Icand J];
                loc.flo=[loc.flo autolocal{J}.Totalfmin(Igood(I))];
                loc.fhi=[loc.fhi autolocal{J}.Totalfmax(Igood(I))];

                loc.duration=[loc.duration autolocal{J}.Totalduration(Igood(I))];
                loc.dt=[loc.dt autolocal{J}.dt(Igood(I))];
                loc.anchor=[loc.anchor strr_org(find(autolocal{J}.dt==0))];
                loc.present{length(loc.dt)}=strr_org(find(autolocal{J}.dt>=0));
                loc.indicies{length(loc.dt)}=autolocal{J}.station_indicies;

            end
            if test_time<loc.closest_time,
                loc.closest_index=J;
            end
            loc.closest_time=min([loc.closest_time test_time]);

            if length(loc.Icand)>2,
                break
            end
        end


    end



end


    function plot_no_morph(auto_station,Igood,Nstation,param,Nrows,strr);

        Inot=setdiff(1:Nstation,Igood);
        Nfig=1;
        Iplot=0;
        myfig=-1;
        set(gcf,'units','normalized','pos', [0.05+(Nfig-1)*0.1 0.05 0.45 0.9]);
        %
        for I=1:length(Inot),  %For every station not in manual result...
            [junk,Ifit]=min(abs(Igood-Inot(I)));
            tlen=20;
            start_ctime=[];
            Itry=Ifit;

            while isempty(start_ctime),
                try,
                    start_ctime=(auto_station(Igood(Ifit)).ctime_debug(min(final(Itry).Icand)))-10;

                end
                Itry=Itry-1;
                if Itry==0,
                    Itry=length(Igood);
                end
                if Itry==Ifit,
                    break;
                end
            end
            try,
                disp(sprintf('loading time: %s',ctime2str(start_ctime)));
                [y,t,head]=readfile(filename{Inot(I)},start_ctime,tlen,1,'ctime','calibrate');

                %%Actually process raw data with morphological processor to view
                %%performance...
                [Iplot,Nfig,myfig]=plot_spectrogram_and_overlay(Iplot,Nfig,y,param,-1,start_ctime,myfig,[],[],[],[],Nrows,strr(I));
                set(gcf,'units','normalized','pos', [0.05+(Nfig-1)*0.1 0.05 0.45/2 0.9]);

                title(sprintf('No manual call: Station %s',strr_org(Inot(I))));
            catch,
                disp('start_ctime failed');
            end
        end

    end
%%%%%%plot_morph_image.m%%%%%%%%%%
%%%plot morphological image output

    function plot_morph_image(TT,FF,station_index,index_offset_t,my_image)

        imagesc(TT,fliplr(FF),sum(my_image,3))
        axis('xy')
        set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');

        for Index=1:length(station_index),
            line(index_offset_t(Index)*[1 1],[0 max(FF)],'color',[1 1 1]);
            text(index_offset_t(Index),10,num2str(station_index(Index)),'color',[1 1 1]);
        end
    end

%%%%%%%%%%%%%%%%%%%%adjust_morph%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%plot and adjust morphological image
    function [param_out,TT,FF,index_offset_t,stored_image,labeled_image]=adjust_morph(auto_station,filename,param,Istation,station_index,bufferTime,show_spec,start_ctime,tlen)
      %  function [param_out,TT,FF,index_offset_t,stored_image,labeled_image]=adjust_morph(Istation,station_index,bufferTime,show_spec,start_ctime,tlen)
   
         stored_image=[];labeled_image=[];
        % TT=[];FF=[];index_offset_t=[];stored_image=[];
        if ~isempty(station_index), %morph result does exist
            [xx,start_ctime,tlen]=extract_signal_from_station(auto_station(Istation),station_index,filename{Istation},'debug');
            index_offset_t=auto_station(Istation).ctime_min(station_index)-start_ctime-param.morph.eq_time;
            %%Download original image
            %fname=dir([automated_dir '/' auto_template_name{Istation} '*morph*']);
            %auto=load([ automated_dir '/' fname.name]);
            %end
            %I_image=min(auto_station(Istation).indicies(1,station_index));
            %stored_image=full(auto.best_calls.features{I_image}(1).final_image);
            for III=1:length(station_index)
                stored_image{III}=full(auto_station(Istation).Image{station_index});
            end

            if III>1,
                disp('multiple station_indicies');
                keyboard;
            end
            param.morph.equalization(:,1)=auto_station(Istation).param.equalization_freq;
            param.morph.equalization(:,2)=mean(auto_station(Istation).equalization(:,station_index),2);

        elseif ~exist('start_ctime');
            disp('no original morphological process result');

            keyboard
        else  %either exists as interval or energy detection..
            disp('morph, Interval or energy detection branch:')
            index_offset_t=bufferTime;
            [xx,t,head]=readfile(filename{Istation},start_ctime,tlen,1,'ctime','calibrate');

        end

        param_out=param;

        change_flag=1;
        while change_flag==1,

            [stats,labeled_image]=extract_image_features(xx,start_ctime,param_out,2);
            if ~isempty(stats),
                FF=stats(1).dF*(0:(size(labeled_image,1)-1));
                TT=stats(1).dT*(0:(size(labeled_image,2)-1));
                tstr='features survive';
            else
                FF=(0:(size(labeled_image,1)-1));
                TT=(0:(size(labeled_image,2)-1));
                tstr='no features survived';
            end

            figure(5);
            subplot(2,1,1)

            if ~isempty(stored_image)
                plot_morph_image(TT,FF,station_index,index_offset_t,stored_image{1});
                title(sprintf('Stored image at %ith index in station',station_index));
            else
                title('Stored image does not exist, no candidate survived...');
                cla;
            end

            if ~isempty(labeled_image),
                subplot(2,1,2)
                plot_morph_image(TT,FF,station_index,index_offset_t,labeled_image);
                title(sprintf('direct from sio, %s, start time (including buffer): %s:%i', ...
                    tstr,datestr(datenum(1970,1,1,0,0,start_ctime),31),round(1000*(start_ctime-floor(start_ctime)))));


                hold off
            else
                disp('No morph returned');
                cla

            end
            if exist('show_spec')&&show_spec>0,
                figure(6);
                spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
                %caxis([-60 80])
                pause;
                close(6)
            end

            [tmp,change_flag]=alter_parameters(param_out.morph);
            param_out.morph=tmp;

        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot_spectrogram_and_overlay.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot a spectrogram and morphological image output side by side.
    function [Iplot,Nfig,myfig]=plot_spectrogram_and_overlay(Iplot,Nfig,y,param,man_ctime,start_ctime,myfig, final,duration,flo,fhi,Nrows,strr)

        tabs=datenum(1970,1,1,0,0,start_ctime);
        Iplot=Iplot+1;
        if Iplot==Nrows+1|Iplot==1,
            Iplot=1;Nfig=Nfig+1;
            figure;myfig=gcf;
            set(gcf,'units','normalized','pos', [0.05+(Nfig-1)*0.1 0.05 0.45 0.9]);

        end
        %keyboard;
        figure(myfig);
        subplot(Nrows,2,2*Iplot-1);
        [SS,TT,FF,PP]=spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
        imagesc(FF,TT,10*log10(PP));axis('xy');
        set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');
        caxis([40 120])
        hold on

        if man_ctime>0,
            %Plot manual detection info
            hh=line(man_ctime+[0 duration ],[flo  flo ]);set(hh,'Color',[0 0 1])
            hh=line(man_ctime + [0 duration ],[fhi  fhi ]);set(hh,'Color',[0 0 1])
            hh=line(man_ctime + [0 0],[flo  fhi ]);set(hh,'Color',[0 0 1])
            hh=line(man_ctime + duration *[1 1],[flo  fhi ]);set(hh,'Color',[0 0 1])
            if ~isempty(final)&&~any(isempty(final.Icand))&&~any(isempty(final.cstart))&&~any(isempty(final.duration))

                for J=1:length(final.flo),
                    hh=line(final.cstart(J)+[0 final.duration(J)],[final.flo(J) final.flo(J)]);set(hh,'Color',[1 0 1])
                    hh=line(final.cstart(J) + [0 final.duration(J)],[final.fhi(J) final.fhi(J)]);set(hh,'Color',[1 0 1])
                    hh=line(final.cstart(J) + [0 0],[final.flo(J) final.fhi(J)]);set(hh,'Color',[1 0 1])
                    hh=line(final.cstart(J) + final.duration(J)*[1 1],[final.flo(J) final.fhi(J)]);set(hh,'Color',[1 0 1])
                    disp(sprintf('Automated call %i, station %s, start time %s, SEL %6.2f',final.Icand(J),strr ,datestr(tabs)',final.SEL(J)));
                end
                title(sprintf('Automated call %i, station %s, start time %s, SEL %6.2f',final.Icand(J),strr ,datestr(tabs)',final.SEL(J)));
            end
        end
        
    end


