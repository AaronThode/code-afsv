function plot_missed_call_stats_and_spectrogram(whale,localized,param,goodFile,goodName,date_range,station)

K=menu('Which DASAR?',goodName);

%for K=1:length(goodFile)
%data=whale{K}.miss.Itrue;
data=whale{K}.miss.morph_step;

%plot_stats;

plot_specgram;

%end



    function plot_stats
        
        
        %%%Plot data statistics...
        
        names={'flo','fhi','duration','stndb'};
        unitts={'Hz','Hz','sec','dB'};
        maxval=[500 500 3 20];
        
        figure
        clear tmp
        for Iname=1:length(names),
            
            
            subplot(length(names),1,Iname)
            [N,F]=hist(data.(names{Iname}),linspace(0,maxval(Iname),50));
            cum=cumsum(N)/sum(N);
            [junk,F95]=min(abs(0.95-cum));
            [junk,F75]=min(abs(0.75-cum));
            
            bar(F,N,'k');
            xlim([0 maxval(Iname)]);
            
            set(gca,'fontweight','bold','fontsize',14);
            ylabel('Number events');xlabel(unitts{Iname});grid on
            title(sprintf(' value of %s on DASAR %s, %s, %6.2f to %6.2f covers 75-95%% range', ...
                names{Iname},goodName{K},date_range,F(F75),F(F95)));
            
            
            
        end
        pause;
        orient landscape
        print('-djpeg',sprintf('MissedCallStats_%s_%s',goodName{K},date_range));
        
        
    end


    function plot_specgram
        figure;
        linewidth=1;
        duration=data.duration;
        flo=data.flo;
        fhi=data.fhi;
        
        for I=1:length(data.ctime),
            % wgt=localized.wgt(data.index(I),K);
            %disp(localized.wgt(data.index(I),:));
            %distok=localized.distok(data.index(I));
            tlen=2*param.energy.bufferTime+duration(I);
            ctime=data.ctime(I);
            
            [y,t,head]=readfile(goodFile{K},ctime-param.energy.bufferTime,tlen,1,'ctime','calibrate');
            
            
            [SS,FF,TT,PP]=spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
            figure(18);set(gcf,'pos',[175         566        1229         493]);
            subplot(2,1,1);
            imagesc(TT,FF,10*log10(PP));axis('xy')
            set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');
            hh=line(param.energy.bufferTime + [0 duration(I)],[flo(I) flo(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
            hh=line(param.energy.bufferTime + [0 duration(I)],[fhi(I) fhi(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
            hh=line(param.energy.bufferTime + [0 0],[flo(I) fhi(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
            hh=line(param.energy.bufferTime + duration(I)*[1 1],[flo(I) fhi(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
            
            %title(sprintf('DASAR %s, %s, flo: %6.2f, fhi: %6.2f , type %i, snr %6.2f dB, Huber weight %6.2f, distok: %i, missed call %i of %i',goodName{K},ctime2str(ctime), ...
            %    data.flo(I),data.fhi(I), data.wctype(I), data.stndb(I), wgt, distok, I,length(data.ctime)));
            
            title(sprintf('DASAR %s, %s, flo: %6.2f, fhi: %6.2f , type %i, snr %6.2f dB,  missed call %i of %i',goodName{K},ctime2str(ctime), ...
                data.flo(I),data.fhi(I), data.wctype(I), data.stndb(I),  I,length(data.ctime)));
            
            
            
            %%%%Reproduce detsum file if needed
            %tic
            %[rawdatadir,outputdir,param,manualdir]=load_pathnames('Shell08',param);
            param.energy.dir_out='.';
            detsum=dir([param.energy.dir_out '/' goodName{K} '*detsum']);
            if I==1
                if length(detsum)==0
                    run_energy_detector(goodFile{K},param.energy);
                    tic;disp(sprintf('JAVA program run in %6.2f seconds',toc));
                    
                    detsum=dir([param.energy.dir_out '/' goodName{K} '*detsum']);
                end
                [data_all,head]=readEnergySummary([param.energy.dir_out '/' detsum.name], Inf);
            end
            
            snips_name=dir([param.energy.dir_out '/' goodName{K} '*.snips']);
            %fclose('all');
            [x,nstarts_snips,npts_snips,feq,eq]=readEnergySnips([param.energy.dir_out '/' snips_name.name], data.snips_index(I),'double','cell');
            
            [SS,FF,TT,PP]=spectrogram(x{1},hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
            subplot(2,1,2);
            imagesc(TT,FF,10*log10(PP));axis('xy');title(sprintf('Index %i from snips file',data.snips_index(I)));
            
            Istrip=round(param.Fs*param.energy.bufferTime);  %Remove tail buffer of raw detection
            
            index=1:(length(x{1})-Istrip);
            cbegin=data_all.ctime(data.snips_index(I))-param.energy.bufferTime;
            
            param.morph.equalization=[feq 10*log10(eq{1})];
            yes=[];
            param_org=param;
            while isempty(yes)
                param=param_org;
                [param,change_flag]=alter_parameters(param_org);
                [features,final_image]=extract_image_features(x{1}(index), cbegin,param,2);  %Bmean is reset
                figure(19)
                subplot(2,1,1)
                try,
                    imagesc(final_image);title(sprintf('recomputed using ctime of %s',ctime2str(cbegin)))
                catch
                    imagesc(sum(final_image,3));title(sprintf('recomputed and summed using ctime of %s',ctime2str(cbegin)))
                    
                end
                try
                    subplot(2,1,2)
                    imagesc(station(K).Image{data.final_index(I)});title(sprintf('Stored value  ctime of %s',ctime2str(ctime)))
                end
                pause;close all;
                
                
                yes=input('Enter 1 to stop: ');
            end
            
            %             %Run energy detector
            %             head=readgsif_header(goodFile{K});
            %             startup_sec=120;
            %             twant=datenum(1970,1,1,0,0,ctime-startup_sec);
            %             tview=[startup_sec-3 startup_sec+3];
            %
            %             yes=input('Review energy detector?');
            %
            %             while ~isempty(yes)
            %                 Energy_detector_review(param,goodFile{K},head.tabs_start,twant,startup_sec+20,tview);
            %                 yes=input('Re-review?');
            %             end
            %%Define target file
            %fname='/Volumes/macmussel1/Arctic_2008/Site_05/S508A0/S508A0T20080821T000000.gsi';
            %
            
            %Define desired time
            %twant=datenum(2008,8,21,0,0,0);
            %twant=datenum(2008,8,21,0,50,0);
            
        end
    end
end
