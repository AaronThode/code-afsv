%%%%%%process_one_unit.m%%%%%%
%% Given a file name and run parameters, extract calls from file
%function [best_calls,best_ctimes,best_durations,raw_detections,interval_detections]=process_one_unit(fname,short_fname,param,run_options);
%Input:
%   fname: Complete filename, including pathname, to be processed (sio
%       file, aka Greeneridge format file)
%   short_fname: short filename, without complete pathname or extension.
%       Used to name energy detector files...
%   param: structure of automated detection parameters, with a param.energy and param.interval_remove substructure.
%       See TOC_param.m for details and template.
%   run_options:  Option flags that include
%       force_energy_detector_run:  If 1, force JAVA energy detector to
%           produce a new 'detsum' file
%       Ncalls_to_sample: number of detections to read into memory for
%           detailed processing at any given time.
%       filtering_stage: 'morph','contour','both','none'
%       debug: includes fields
%
%Output:
%      best_calls{I}(J) is a structure array of features of shape J from detection I.
%      best_ctimes is a vector of start times of when 'snips' file data
%       start.  Note that these times are the detection times minus the
%       buffer time inserted in the snips file.
%           To get precise ctime of detection use
%           best_calls{I}(J).ctime;
%           best_ctimes should be most useful for
%           debugging (reconstructing morphological processing.
%      durations: a vector of durations (data_all.npt/Fs) that represent
%           time window used for morphological processing.  Output for
%           potential debugging purposes.
%   raw_detections: First stage output from energy detector containing
%       .ctime
%       .duration
%   interval_detections: Second stage output, same structure as
%   'raw_detections'

%Revised Sept. 6, 2008
% Extensive revision Jan 21, 2010

function [best_calls]=process_one_unit(fname, short_fname,param,run_options,manual,Istation)

debug_params=run_options.debug;
best_calls=[];
best_ctimes=[];
best_durations=[];
raw_ctimes=[];
raw_detections=[];
interval_detections=[];

run_interval_stage=(isfield(param,'interval_remove')&&param.interval_remove.on==1);
if isempty(fname)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Stage one:  run CFAR cell detector, aka 'energy' detector using JAVA
%%%program, appending bearing calculation at end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_options.force_energy_detector_run==1
    tic
    
    if run_interval_stage&&strcmp(param.interval_remove.ICI_removal_alg,'BearingAndTiming')
        %%Force snips_chc to be 2, getting bearing information
        param.energy.snips_chc=2;
    else
        param.energy.snips_chc=1;
    end
    
    run_energy_detector(fname,param.energy);
    disp(sprintf('JAVA program run in %6.2f seconds',toc));
    
end
tic

%%Upload all raw detections from JAVA *detsum files into memory%%
detsum=dir([param.energy.dir_out '/' short_fname '*detsum']);
[data_all,head]=readEnergySummary([param.energy.dir_out '/' detsum.name], Inf);


%%Compare detections with manual results...


auto.tabs=datenum(1970,1,1,0,0,data_all.ctime);
auto.duration=data_all.features(9,:);
auto.ctime=data_all.ctime;
ovlap_tol=0.25;

Iok=find(manual.tabs(:,Istation)>0&manual.tabs(:,Istation)<=max(auto.tabs));
tabs_m=manual.tabs(Iok,Istation);
ctime_m=manual.ctime(Iok,Istation);
dur_m=manual.duration(Iok,Istation);
[Iauto_nomatch,Imanual_nomatch,Iauto_match,Imanual_match]=find_similar_elements_ovlap(auto.ctime,auto.duration, ...
    ctime_m,dur_m,ovlap_tol);


%%Remove repeats in Iauto_match;
[Iauto_match,Iunique]=unique(Iauto_match);
Igood=Iauto_match;
Imanual_match=Imanual_match(Iunique);


fprintf('Out of %i manual detections, %i were missed and %i were matched\n',length(tabs_m),length(Imanual_nomatch),length(Imanual_match));
fprintf('Out of %i auto detections, %i were missed and %i were matched\n',length(auto.tabs),length(Iauto_nomatch),length(Iauto_match));

Imanual_nomatch_highSNR=find(manual.stndb(Iok(Imanual_nomatch),Istation)>=param.energy.threshold);
Imanual_nomatch_highSNR2=find(manual.stndb(Iok(Imanual_nomatch),Istation)>=param.energy.threshold+3);

fprintf('Out of %i missing manual detections, %i were missed that were greater than the CFAR %6.2f dB threshold \n', ...
    length(Imanual_nomatch),length(Imanual_nomatch_highSNR), param.energy.threshold);
fprintf('Out of %i missing manual detections, %i were missed that were greater than the CFAR %6.2f dB +3 dB threshold \n', ...
    length(Imanual_nomatch),length(Imanual_nomatch_highSNR2), param.energy.threshold);

subplot(3,1,1)
hist(manual.stndb(Iok,Istation),0:25); title(sprintf('Original SNR of all manual detections: %s',short_fname));xlim([0 25])
subplot(3,1,2)
hist(manual.stndb(Iok(Imanual_nomatch),Istation),0:25);title('SNR of missed manual detections');xlim([0 25])
subplot(3,1,3)
hist(manual.wctype(Iok(Imanual_nomatch),Istation),0:8);title('Call type of missed manual detections');xlim([0 8])
pause(2);close

manual.ctime=manual.ctime(Iok(Imanual_match),Istation);
manual.tabs=manual.tabs(Iok(Imanual_match),Istation);
manual.duration=manual.duration(Iok(Imanual_match),Istation);
manual.wctype=manual.wctype(Iok(Imanual_match),Istation);
manual.fhi=manual.fhi(Iok(Imanual_match),Istation);
manual.flo=manual.flo(Iok(Imanual_match),Istation);
manual.sigdb=manual.sigdb(Iok(Imanual_match),Istation);
manual.stndb=manual.stndb(Iok(Imanual_match),Istation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Conduct image processing and feature extraction%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%raw_ctimes=data_all.ctime(Igood);
debug_params.rawfile=fname;  %Pass input filename to routine for debug purposes
debug_params.short_fname=short_fname;
Ncalls=run_options.Ncalls_to_sample;
tic

% save temp
% keyboard

process_detection_spectrogram;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Nested functions%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%eval(sprintf('!mv *detsum %s/',detsumdir));
%!rm *snips *detsum


    function process_detection_spectrogram
        %function process_detection_spectrogram(Igood,data_all,head,Ncalls,param, alg_chc,debug_params)
        %% Load snippets from energy detection output and pass to morphological
        %%      processing algorithm...
        %
        %   Igood: vector of indicies corresponding to elements in data_all
        %   data_all: structure array of raw detections, output of
        %           readEnergySummary.m detector
        %   head: structure summarizing summary file data, output of
        %           readEnergySummary.m
        %   Ncalls: number of detections to read into RAM per processing loop
        %   param: detection parameters as provided by TOC_params.
        %   alg_chc:  'contour','morph','both'; c
        %   debug_params:  If exists, structure of debug flags.  Has fields
        %   'rawfile','snips','tol','ctimes_debug', 'short_fname', and 'morph'
        %   Output:
        
        
        errormsg=[];
        
        best_calls.param=param;
        best_calls.run_options=run_options;
        best_calls.station=Istation;
        best_calls.detsum_name=detsum.name;
        best_calls.short_fname=short_fname;
        Ngood=length(Igood);
        best_calls.wctype=zeros(1,Ngood);
        best_calls.ctime=zeros(1,Ngood);
        best_calls.tabs=zeros(1,Ngood);
        best_calls.flo=zeros(1,Ngood);
        best_calls.fhi=zeros(1,Ngood);
        best_calls.duration=zeros(1,Ngood);
        best_calls.stndb=zeros(1,Ngood);
        best_calls.nstart=zeros(1,Ngood);
        best_calls.manual_index=zeros(1,Ngood);
        best_calls.auto_index=zeros(1,Ngood);
        
        Icall=0;
        Fs=head.Fs;
        
        try
            best_ctimes=zeros(1,Ngood);
            best_nstart=best_ctimes;
            best_npt=best_ctimes;
            
            %Istrip=round(Fs*param.energy.bufferTime);  %Remove tail buffer of raw detection
            Istrip=1;
            if Ngood<Ncalls
                Ncalls=Ngood;
            end
            
            for I=1:ceil(Ngood/Ncalls)
                if rem(I,5)==0,disp(sprintf('%6.2f percent of calls processed',100*(I-1)*Ncalls/Ngood));end
                Iabs=Ncalls*(I-1)+(1:Ncalls);
                Iabs=Iabs(Iabs<=Ngood);  %In case desired calls less than Ncalls
                Iref=Igood(Iabs);  %%Index number of snips file to call
                snips_name=dir([param.energy.dir_out '/' debug_params.short_fname '*.snips']);
                
                %%%Read in time series 'x', equalization spectrum 'eq' and
                %%%associated  frequencies 'feq'
                [x,~,~,feq,eq]=readEnergySnips([param.energy.dir_out '/' snips_name.name], Iref,'double','cell','keep_open');
                
                if isempty(x{1}),
                    disp('Snips file not read, empty matrix returned');
                    return;
                end
                for II=1:length(Iref)
                    
                    if run_interval_stage&strcmp(param.interval_remove.ICI_removal_alg,'BearingAndTiming')  %Keep only omni channel
                        x{II}=x{II}(1,:);
                    end
                    index=1:(length(x{II})-Istrip);
                    
                    %%%Define start time of snips--including buffer time
                    cbegin=data_all.ctime(Iref(II))-param.energy.bufferTime;
                    
                    %%%Check for a debug trigger....
                    if debug_params.morph&isempty(debug_params.ctimes_debug)
                        debug_flag=1;
                    elseif min(abs(cbegin-debug_params.ctimes_debug))<=debug_params.tol
                        disp('debug flag triggered');
                        debug_flag=1;
                    else
                        debug_flag=0;
                        
                    end
                    
                    if run_options.plot_manual_data==1
                        wctypes={'Upsweep','Downsweep','Constant','U-shaped','N-Shaped','Undulation','Complex','Bearded Seal'};
                        %tlen=param.energy.bufferTime+data_all.npt(1,Iref(II))/Fs;
                        Nfft=256;
                        ovlap=round(0.9*Nfft);
                        
                        [S,F,T,PP2]=spectrogram(x{II}(1,index),Nfft,ovlap,Nfft,head.Fs,'yaxis');
                        imagesc(T,F,10*log10(abs(PP2)));axis('xy');
                        Itype=manual.wctype(Iabs(II));
                        title(sprintf('Manual Time: %s Auto Time: %s wctype: %s:%6.2f, flo: %6.2f, fhi: %6.2f, duration %6.2f, SNR: %6.2f',datestr(manual.tabs(Iabs(II))), ...
                            ctime2str(cbegin),wctypes{Itype},Itype,manual.flo(Iabs(II)),manual.fhi(Iabs(II)), ...
                            manual.duration(Iabs(II)), manual.stndb(Iabs(II))));
                        caxis([40 100])
                        tplot=manual.ctime(Iabs(II))-cbegin+[0 manual.duration(Iabs(II))];
                        linewidth=2;
                        hh=line(tplot,manual.flo(Iabs(II))*[1 1]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
                        hh=line(tplot,manual.fhi(Iabs(II))*[1 1]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
                        hh=line(tplot(1)*[1 1],[manual.flo(Iabs(II)) manual.fhi(Iabs(II))]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
                        hh=line(tplot(2)*[1 1],[manual.flo(Iabs(II)) manual.fhi(Iabs(II))]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
                        text(1,300,wctypes{Itype},'fontweight','bold','fontsize',24,'color',[1 1 1 ]);
                        
                        colorbar;
                        disp(sprintf('I is %i, II is %i, Iabs(II) is %i, Iref is %i,',I,II,Iabs(II),Iref(II)));
                        pause;
                    end
                    
                    %if not empty, store data.
                    Icall=Icall+1;
                    best_calls.x{Icall}=x{II}(1,index);
                    best_calls.wctype(Icall)=manual.wctype(Iabs(II));
                    best_calls.ctime(Icall)=cbegin;
                    best_calls.tabs(Icall)=manual.tabs(Iabs(II));
                    best_calls.flo(Icall)=manual.flo(Iabs(II));
                    best_calls.fhi(Icall)=manual.fhi(Iabs(II));
                    best_calls.duration(Icall)=manual.duration(Iabs(II));
                    best_calls.stndb(Icall)=manual.stndb(Iabs(II));
                    best_calls.sigdb(Icall)=manual.sigdb(Iabs(II));
                    best_calls.bearing(Icall)=manual.bearing(Iabs(II));
                    best_calls.wgt(Icall)=manual.wgt(Iabs(II));
                   
                    best_calls.nstart(Icall)=data_all.nstart(1,Iref(II));
                    best_calls.manual_index(Icall)=Iok(Imanual_match(Iabs(II)));
                    best_calls.auto_index(Igood)=Iref(II);
                    best_calls.equalization_freq=feq;
                    if Icall==1
                        best_calls.equalization=zeros(length(eq{II}),Ngood);
                    end
                    best_calls.equalization(:,Icall)=10*log10(eq{II});
                    
                    
                    
                    
                end %II, call mark, call
            end  %I -- calls
            
            fclose('all');
            if Icall~=0
                %best_npt=best_npt(1:Icall);
                disp(sprintf('Out of %i detections reviewed, %i pass image processing, or %5.4f percent pass',Ngood,Icall,100*Icall/Ngood));
                
            end
            
        catch
            errormsg=lasterror;
            disp(errormsg.message);
            errormsg.stack.file
            errormsg.stack.name
            errormsg.stack.line
           
            
            fclose('all');
            
            return
            
        end
        
    end


end

