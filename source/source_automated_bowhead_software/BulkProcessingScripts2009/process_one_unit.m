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

function [best_calls,best_ctimes,best_durations,raw_detections,interval_detections]=process_one_unit(fname, ...
    short_fname,param,run_options)

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


%% If bearing information to be used in interval detector,
%%%compute bearing of all detections using readEnergySnips

if run_interval_stage&&strcmp(param.interval_remove.ICI_removal_alg,'BearingAndTiming')
    
    %%%Compute bearing of all imports using readEnergySnips
    thet=extract_snips_bearing(data_all,run_options,param,short_fname);
    %%Add bearing information to data_all and update data_all.names
    data_all.features=[data_all.features; thet];
    Nfeature=size(data_all.features,1);
    data_all.names{Nfeature}='bearing';
    
    %Provide an opportunity to plot the bearing, ICI, and frequency range of
    %raw detections
    if run_options.debug.interval==2
        plot_energy_result(1:length(data_all.ctime),'x');
    end
    
end
tmin=datestr(datenum(1970,1,1,0,0,data_all.ctime(1)));
tmax=datestr(datenum(1970,1,1,0,0,data_all.ctime(end)));

disp(sprintf('Raw data loaded in %6.2f seconds, min time %s max time %s',toc, tmin, tmax));
toc

if isempty(data_all.ctime),
    return
end

%%Create raw_detections output variable..
raw_ctimes=data_all.ctime;
raw_detections.ctime=raw_ctimes;
raw_detections.duration=data_all.npt(1,:)/param.Fs;

%%Check if debug ctimes were detected by energy detector...
Imatch_raw=debug_raw_times;
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Stage two:  run interval detector to remove regular intervals if desired%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(param,'interval_remove')||param.interval_remove.on==0
    Igood=1:length(data_all.ctime);
else
    %%Prepare debug parameters if desired.
    if run_options.debug.interval==1;
        interval_debug.fname=fname;
        interval_debug.names=param.interval_remove.names;
        interval_debug.index=param.interval_remove.debug_index;
        if isfield(param.interval_remove,'debug_time')
            interval_debug.time=param.interval_remove.debug_time;
        end
        
    else
        interval_debug=[];
    end
    
    %%Identify what features to use when filtering ICI
    Iwant=[];
    for JJ=1:length(param.interval_remove.names),
        for KK=1:length(data_all.names),
            if strcmp(data_all.names{KK},param.interval_remove.names{JJ})
                Iwant=[Iwant KK];
            end
        end
    end
    
    %%Compute ICI (regular intervals) from raw data, using selected features
    %%for assistance.
    
    ICI=compute_ici_bothways_feature(datenum(1970,1,1,0,0,data_all.ctime),[5 42],data_all.features(Iwant,:), param.interval_remove.names,...
        param.interval_remove.Ndet,param.interval_remove.Nmiss,...
        param.interval_remove.tol_feature,param.interval_remove.ICItol,interval_debug);
    
    
    %%Added April 5, 2010:  Sometimes walrus sequences and heavy bowhead whale
    %%sequences can have an ICI, but the ICI is inconsistent between calls.
    %  Thus here we march through each ICI detection and check whether detections
    %  nearby share the same ICI.
    
    Iguns=find(ICI>0);
    
    if run_options.debug.interval==2
        hold on
        plot_energy_result(Iguns,'go');
    end
    ICI_score=ones(size(ICI));
    for I2=1:length(Iguns)
        current_time=data_all.ctime(Iguns(I2));
        current_ICI=ICI(Iguns(I2));
        Itest=find(abs(current_time-data_all.ctime(Iguns))<=0.5*param.interval_remove.Ndet*current_ICI);
        Ipass=0;
        for J=1:2  %Harmonic loop:  checks for possibility that a 10 s ICI may have been assigned a 20 s ICI.
            Ipass=Ipass+length(find(abs(ICI(Iguns(Itest))/J-current_ICI)<=param.interval_remove.ICI_std));
        end
        
        %%Are there enough matching ICIs close to the value of the current ICI?
        if (Ipass-1)<param.interval_remove.Nstd
            ICI_score(Iguns(I2))=0;
        else
            %disp('good');
        end
        
    end
    
    
    %%Assign airgun labels to signals with ICIs that are consistent with other
    %%ICI signals... note that the ICI variable still retains its original
    %%output..
    
    Iguns=find(ICI_score.*ICI>0);
    Igood=setdiff(1:length(ICI),Iguns);
    
    %%Assign output variable interval_detections..
    
    interval_ctimes=raw_ctimes(Igood);
    interval_detections.ctime=interval_ctimes;
    interval_detections.duration=data_all.npt(1,Igood)/param.Fs;
    
    fprintf('After %6.2f sec interval check %i regular signals, %i out of %i detections remain.\n',toc, ...
        length(Iguns),length(Igood),length(data_all.ctime));
    Igood=debug_interval_times;
    
end
%%%%%%%End interval detector%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Conduct image processing and feature extraction%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%raw_ctimes=data_all.ctime(Igood);
debug_params.rawfile=fname;  %Pass input filename to routine for debug purposes
debug_params.short_fname=short_fname;
Ncalls=run_options.Ncalls_to_sample;
alg_chc=run_options.filtering_stage;
tic



process_detection_spectrogram;

fprintf('Image processing stage finished in %6.2f seconds or %6.2f detections/minute\n',toc,length(Igood)*60/toc);tic

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
        
        best_calls.features=[];
        best_calls.ctimes=[];
        best_ctimes=[];
        best_durations=[];
        errormsg=[];
        if ~isfield(debug_params,'ctimes_debug')
            debug_params.ctimes_debug=[];
            debug_params.tol=0;
        end
        disp(sprintf('%s method chosen',alg_chc));
        Icall=0;
        Fs=head.Fs;
        
        try
            best_ctimes=zeros(1,length(Igood));
            best_nstart=best_ctimes;
            best_npt=best_ctimes;
            
            Istrip=round(Fs*param.energy.bufferTime);  %Remove tail buffer of raw detection
            
            if length(Igood)<Ncalls,
                Ncalls=length(Igood);
            end
            
            for I=1:ceil(length(Igood)/Ncalls)
                if rem(I,5)==0,disp(sprintf('%6.2f percent of calls processed',100*(I-1)*Ncalls/length(Igood)));end
                Iabs=Ncalls*(I-1)+(1:Ncalls);
                Iabs=Iabs(Iabs<=length(Igood));  %In case desired calls less than Ncalls
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
                    
                    
                    %%%%check snip files vs direct download...
                    %debug_flag=1;
                    %disp(II);
                    %disp(abs(debug_params.ctimes_debug-cbegin));
                    if debug_params.snips==1&&debug_flag==1  %check snip files vs direct download...
                        disp(sprintf('Plotting snips of debug time %s',ctime2str(debug_params.ctimes_debug)));
                        tlen=param.energy.bufferTime+data_all.npt(1,Iref(II))/Fs;
                        
                        [y,t,head]=readfile(debug_params.rawfile,cbegin,tlen,1,'ctime','calibrate');
                        format long e
                        disp(data_all.nstart(Iref(II)));
                        [S,F,T,PP1]=spectrogram(y,128,96,128,head.Fs,'yaxis');
                        figure(3);set(gcf,'units','norm','pos',[6.666e-02     4.858333e-01     2.91666e-01     3.5000e-01]);
                        subplot(2,1,1);
                        imagesc(T,F,10*log10(abs(PP1)));axis('xy');
                        set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');
                        colorbar
                        title(sprintf('direct from sio, start time (including buffer): %s:%i', datestr(datenum(1970,1,1,0,0,cbegin),31),round(1000*(cbegin-floor(cbegin)))));
                        subplot(2,1,2);
                        [S,F,T,PP2]=spectrogram(x{II}(1,index),128,96,128,head.Fs,'yaxis');
                        imagesc(T,F,10*log10(abs(PP2)));axis('xy');
                        title(sprintf('snips file, II is %i, length of detection %6.2f sec',Iref(II),tlen-param.energy.bufferTime));
                        colorbar;
                        keyboard;
                    end
                    
                    
                    %%%Set up inputs to morphological processor, including
                    %%% equalization function...
                    param.morph.equalization=[feq 10*log10(eq{II})];
                    best_calls.equalization_freq=feq;
                    if debug_params.morph|debug_flag,
                        disp(sprintf('cbegin is %s, %16.15f',ctime2str(cbegin),cbegin));
                    end
                    
                    %%Extract features from workhorse algorithm
                    [features,final_image]=extract_image_features(x{II}(1,index), cbegin,param,2*(debug_flag));
                   
                    if ~isempty(features),
                        Icall=Icall+1;
                        best_nstart(Icall)=data_all.nstart(1,Iref(II));
                        best_npt(Icall)=data_all.npt(1,Iref(II));
                        best_calls.tstart(:,Icall)=data_all.nstart(Iref(II))/Fs;
                        best_ctimes(Icall)=cbegin;
                        best_calls.ctime(Icall)=best_ctimes(Icall);
                        %best_calls.duration(:,Icall)=data_all.npt(Iref(II))/Fs;
                        best_calls.duration(:,Icall)=param.energy.bufferTime+data_all.npt(Iref(II))/Fs;  %%AARON RECENT CHANGE
                        best_durations(Icall)=best_calls.duration(1,Icall);
                        best_calls.index(Icall)=Iref(II);
                        best_calls.equalization{Icall}=param.morph.equalization(:,2);
                        %Speciality to morph..
                        for Im=1:size(final_image,3);
                            best_calls.labeled_image{Icall}{Im}=sparse(squeeze(final_image(:,:,Im)));
                        end
                        best_calls.features{Icall}=features;
                        if run_interval_stage
                            best_calls.ICI(Icall)=ICI(Iref(II));
                            best_calls.bearing(Icall)=thet(Iref(II));
                        end
                        
                    elseif  debug_params.morph|debug_flag
                        disp(sprintf('process_detection_spectrogram: No features passed extract_image_features'));
                        
                        
                    end
                    
                    
                    
                end %II, call mark, call
            end  %I -- calls
            
            fclose('all');
            if Icall~=0,
                best_ctimes=best_ctimes(1:Icall);
                %best_nstart=best_nstart(1:Icall);
                %best_npt=best_npt(1:Icall);
                disp(sprintf('Out of %i detections reviewed, %i pass image processing, or %5.4f percent pass',length(Igood),Icall,100*Icall/length(Igood)));
            else
                best_calls.features=[];
                best_ctimes=[];
            end
            
        catch
            errormsg=lasterror;
            disp(errormsg.message);
            errormsg.stack.file
            errormsg.stack.name
            errormsg.stack.line
            keyboard;
            
            fclose('all');
            
            return
            
        end
        
    end


%%%%%%%plot_energy_result.m%%%%%%%%%%
    function plot_energy_result(Irange,ssymbol)
        tabs=datenum(1970,1,1,0,0,data_all.ctime(Irange));
        figure(1)
        subplot(4,1,1)
        plot(tabs,thet(Irange),ssymbol);hold on
        xlabel('Time');ylabel('bearing');
        datetick('x',14)
        grid on
        Iwant=[1 5];
        F=param.energy.f_low:(param.energy.bandwidth/2):250;
        Th=(0:10:360);
        
        %     for JJ=1:length(Iwant);
        %         subplot(3,1,JJ+1)
        %         HH=zeros(length(Th)-1,length(F)-1);
        %         for It=1:(length(Th)-1)
        %             Igood=find(thet>=Th(It)&thet<=Th(It+1));
        %             for If=1:(length(F)-1)
        %                 Igood2=find(data_all.features(Iwant(JJ),Igood)>=F(If)&data_all.features(Iwant(JJ),Igood)<=F(If+1));
        %                 HH(It,If)=log10(length(Igood2));
        %             end
        %         end
        %         imagesc(F,Th,HH);
        %         axis('xy');
        %
        %     end
        subplot(4,1,2)
        plot(tabs,data_all.features(1,Irange),ssymbol);ylabel('Min freq (Hz)');xlabel('Time');hold on
        grid on;datetick('x',14)
        subplot(4,1,3)
        plot(tabs, data_all.features(5,Irange),ssymbol);ylabel('Peak freq (Hz)');xlabel('Time');hold on
        grid on;datetick('x',14)
        if exist('ICI')
            subplot(4,1,4)
            plot(tabs, ICI(Irange),ssymbol);ylabel('ICI(sec)');xlabel('Time');hold on
            grid on;datetick('x',14)
        end
        
        %keyboard
    end

%%%%%%%debug_interval_times.m%%%%%%%%%%
    function Igood_out=debug_interval_times
        
        if ~isempty(interval_debug)
            figure;
            tabs=datenum(1970,1,1,0,0,data_all.ctime(Iguns));
            plot(tabs,ICI(Iguns),'kx');
            datetick('x',14);
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Local time');
            ylabel('Interval (sec)');
            grid on;
            ylim([0 40]);
            keyboard;
            
        end
        
        Igood_out=Igood;
        if debug_params.raw>0&&~isempty(debug_params.ctimes_debug),
            disp('');
            for JJ=1:length(Imatch_raw),
                disp(sprintf('ICI for raw detection %i is %6.2f',Imatch_raw(JJ),ICI(Imatch_raw(JJ))));
                %
            end
            
            [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match]=find_similar_elements(debug_params.ctimes_debug, ...
                data_all.ctime(Igood),debug_params.tol);
            disp(sprintf('%i ctimes present after ICI stripping',length(Iy_match)));
            
            
            if length(Iy_match)==0,
                
                yes=input('Activating debug ctime to see if it survives further processing?');
                
                if ~isempty(yes)
                    %                     ICI_temp=compute_ici_bothways_feature(datenum(1970,1,1,0,0,data_all.ctime),[5 40],data_all.features(Iwant,:), ...
                    %                         param.interval_remove.Ndet,param.interval_remove.Nmiss,...
                    %                         param.interval_remove.tol_feature,param.interval_remove.ICItol,interval_debug);
                    
                    ICI(Imatch_raw(JJ))=-1;  %permit this call to be checked in future
                    Igood_out=find(ICI<0);
                end
                
            end
        end
        
    end

%%%%%%%%%%debug_raw_times.m%%%%%%%%%%
    function Imatch_raw=debug_raw_times
        Imatch_raw=[];
        if debug_params.raw>0&&~isempty(debug_params.ctimes_debug),
            [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match]=find_similar_elements(debug_params.ctimes_debug, ...
                raw_ctimes,debug_params.tol);
            disp(sprintf('%i ctimes present out of %i in raw energy data',length(Iy_match),length(debug_params.ctimes_debug)));
            Imatch_raw=Iy_match;
            for JJ=1:length(Iy_match),
                disp(sprintf('Raw index %i at time %s',Iy_match(JJ),ctime2str(raw_ctimes(Iy_match(JJ)))));
            end
            %keyboard;
        end
    end
end

