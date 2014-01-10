%%%%%%extract_transient_levels.m%%%%%%%%%%%%%%%%%%
% Routine used to extract level information from clips in *.snips file.
%   Key feature is the precomputation of a series of FIR hi and lo pass
%   filters to allow quick filtering of raw acoustic data before extracting
%   the metrics.
%
%  Aaron Thode
%  Dec. 29, 2013
%
% Inputs:
%       snips_name: full pathname to snips file.
%       Igood: indicies of data_all that contain desired signals.
%           Typically contains indicies of signals that have passed ICI test.
%       data_all: Output of a *detsum file.  Contains fields
%               .ctime: vector of ctimes of detections...
%               .npt:  row vector of duration of signals in samples.
%           Note: data_all contains all data from *.snips file, and Igood
%               is used to select subset.
%       Fs (Hz): sampling rate of raw data.
%       bufferTime(s): typically value of param.energy.bufferTime
%               Note that enough bufferTime is needed to account for FIR
%               filtering...
%       filter_transition_band (Hz):  Scalar: width of transition band for
%               precomputed FIR filters. Must be less than half of the
%               minimum spacing between filter_passband_frequencies
%               discussed below.
%       filter_passband_frequencies (Hz): An array of frequencies to
%               compute hi and lopass filters:  Define frequencies over
%               which signal can be filtered before computing levels.
%       run_options:
%               .debug:  If exists, start showing plots at tranisent
%                       detection number given by .debug.
%               Ncalls_to_sample:number of signals to load into RAM memory for
%                   processing.
%               calibration_keyword:  Used by calibrate_GSI_signal.m to
%                   select calibration scheme for raw data...
%
% Outputs:
%        transient_params:
%               success: '1' or '-1';
%function [transient_params]=extract_transient_levels(snips_name,Igood,data_all,Fs,bufferTime, ...
%filter_transition_band,filter_passband_frequencies,run_options)

function [transient_params]=extract_transient_levels(snips_name,Igood,data_all,Fs,bufferTime, ...
    filter_transition_band,filter_passband_frequencies,run_options)
persistent  B_dfilt_hipass B_dfilt_lopass Imaxx Iminn Ncoef

transient_params.success=-1;
temp=dir(snips_name);
if isempty(temp)
    errordlg(sprintf('%s does not exist!',snips_name));
    return
end

bufferTime_org=bufferTime;

if isempty(Igood)
    return
end
%max_clip_count=run_options.max_clip_count;  %Number of clips allowed before data marked as tainted..

debug_snips_check=sign(run_options.debug);

%transient_params.f_center=[10 20 32 40 50 63 80 100 125 160 200 250 315 395]; %1/3 octave levels  %Hz

%%Error checking of input variables
if ~exist('filter_passband_frequencies')||length(filter_passband_frequencies)<2
    disp('extract_transient_levels:  Need to enter filter_passband_frequencies');
    return
end

df=min(diff(filter_passband_frequencies));
if 0.5*df<filter_transition_band
    errordlg('extract_transient_levels:  transition_band too large relative to passband frequency spacing');
    return
end

%%Decide whether to decimate if signal sampling rate is much higher than
%%  desired measuring bandwidth
%%
%if ~isempty(findstr(run_options.short_fname,'AURAL'))
if 10*max(filter_passband_frequencies)<Fs/2
    yes_decimate=1;
    Fs_filt=Fs/10;
else
    Fs_filt=Fs;
    yes_decimate=0;
end


transient_params.f_passband=filter_passband_frequencies;
%transient_params.f_passband(end)=(Fs_filt/2)-1;
%transient_params.f_passband=unique(transient_params.f_passband);

Npass=length(transient_params.f_passband);


if length(Igood)<run_options.Ncalls_to_sample
    run_options.Ncalls_to_sample=length(Igood);
end

%scale_factor=get_calibration_scale_factor(run_options.calibration_keyword);
scale_factor=1;

%% prefiltering loop if hasn't been defined yet.
%param.airgun.filter_bandwidth
yes_filter_recompute=isempty(B_dfilt_lopass)&&~isempty(transient_params.f_passband); 
yes_filter_recompute= yes_filter_recompute|| length(B_dfilt_lopass)~=length(transient_params.f_passband); %If new import
if yes_filter_recompute
    fprintf('Recomputing filters: size of B_dfilt_lopass: %i, length of input parameters: %i\n', ...
        length(B_dfilt_lopass),length(transient_params.f_passband));
    
    %Flag feature for max frequency
    for Ifea=1:length(data_all.names)
        if strcmp(data_all.names{Ifea},'max_freq')
            Imaxx=Ifea;
        elseif strcmp(data_all.names{Ifea},'min_freq')
            Iminn=Ifea;
        end
    end
    
    hh	=	waitbar(0,sprintf('Creating filters ...'));
    
    
    for Ifilter=1:length(transient_params.f_passband)
        waitbar(Ifilter/length(transient_params.f_passband),hh);
        filter_bandwidth(1)=transient_params.f_passband(Ifilter);
        filter_bandwidth(2)=filter_bandwidth(1)+ filter_transition_band;
        %filter_bandwidth(3)=transient_params.f_passband(Ifilter);
        %filter_bandwidth(4)=filter_bandwidth(3)+filter_bandwidth(2)-filter_bandwidth(1);
        disp(sprintf('Building %.0f-%.0f Hz transition band filters',filter_bandwidth(1),filter_bandwidth(2)));
        %disp(sprintf('Transition zone bandwidth: %6.2f',filter_bandwidth(2)-filter_bandwidth(1)));
        
        B_dfilt_lopass{Ifilter}=brickwall_lpf(filter_bandwidth,Fs_filt,0);
        B_dfilt_hipass{Ifilter}=brickwall_hpf(filter_bandwidth,Fs_filt,0);
        Ncoef.lopass(Ifilter)=length(B_dfilt_lopass{Ifilter}.Numerator);
        Ncoef.hipass(Ifilter)=length(B_dfilt_hipass{Ifilter}.Numerator);
        
        filt_time=(Ncoef.lopass(Ifilter)+Ncoef.hipass(Ifilter))/Fs_filt;
        if (filt_time>bufferTime_org) 
            errordlg(sprintf('Filter length %6.2f sec too large for snips bufferTime of %6.2f sec: increase transition band and/or bufferTime:', ...
                filt_time,bufferTime_org));
            close(hh);
            return
        else
            fprintf('Total filter time is %6.2f sec\n',filt_time);
        end
        
    end
    close(hh);
else
    disp('Filters already built');
end
fclose('all');

Icall=0;
transient_params.success=1;  %OK, made it to processing loop

hh	=	waitbar(0,sprintf('Importing snips file...'));


%Preallocatie
Nsnips=length(Igood);
AA=zeros(1,Nsnips);
transient_params.level.max_freq=AA;
transient_params.level.min_freq=AA;
transient_params.ctime =AA;
transient_params.index =AA;
transient_params.clipped =AA;
transient_params.level.peak =AA;
%transient_params.level.peakF =features.peakF;
transient_params.level.t_Malme =AA;
transient_params.level.SE_Malme =AA;
transient_params.level.rms_Malme =AA;
transient_params.noise.rms =AA;
transient_params.noise.SE =AA;
transient_params.noise.duration =AA;
transient_params.level.SNR =AA;
         
        
for I=1:ceil(Nsnips/run_options.Ncalls_to_sample)
    disp(sprintf('Batch %i of %i, Icall %i',I,ceil(length(Igood)/run_options.Ncalls_to_sample),Icall));
    Iabs=run_options.Ncalls_to_sample*(I-1)+(1:run_options.Ncalls_to_sample);
    Iabs=Iabs(Iabs<=length(Igood));
    Iref=Igood(Iabs);
    %snips_name=dir([dir_out '/' run_options.short_fname '*.snips']);
    
    %Load data from JAVA program.  Note that a 4th-order Butterworth filter has already been applied, but
    %clipping events have been retained...
    [x,~,~,~,~,~,head]=readEnergySnips(snips_name, Iref,'double','cell','keep_open');
    %toc(tt1)
    
    
    % disp(sprintf('%6.2f percent done',100*I/(length(Igood)/run_options.Ncalls_to_sample)));
    cstart=data_all.ctime(Iref);
    dstart=diff(cstart);dstart=[0 dstart cstart(end)];
    
    for II=1:length(Iref)
        Icall=Icall+1;
        
         if rem(II,100)==0
              waitbar(Icall/length(Igood),hh);
    
         end
         
        max_freq=data_all.features(Imaxx,Iref(II));
        min_freq=data_all.features(Iminn,Iref(II));
        transient_params.level.max_freq(Icall)=max_freq;
        transient_params.level.min_freq(Icall)=min_freq;
        
        [junk,Ifilt]=min(abs(max_freq-transient_params.f_passband));
        if max_freq>transient_params.f_passband(Ifilt)
            Ifilt=Ifilt+1;
            if Ifilt>Npass
                fprintf('extract_transient_levels: maximum frequency of transient is greater than maximum passband frequency\n');
                Ifilt=Npass;
            end
        end
        
        [junk,Ifilt2]=min(abs(min_freq-transient_params.f_passband));
        if min_freq<transient_params.f_passband(Ifilt2)
            Ifilt2=Ifilt2-1;
            if Ifilt2<1
                fprintf('extract_transient_levels: minimum frequency of transient is less than minimum passband frequency\n');
                Ifilt=1;
            end
            
        end
        
        %%If bandwidth of target signal is less than spacing of
        %%passband_frequencies, force them to be different.
        if Ifilt==Ifilt2
            Ifilt2=Ifilt-1;
        end
        
        if length(x{II}(1,:))<round(Fs* bufferTime)
            disp('signal shorter than buffer time.');
            % keyboard
            continue;
        end
        
        %%%%Read raw data to check for clipping
        %%%Note that even a filtered file, like a DASARC, will assign a
        %%%clipped value to output even original input data are clipped.
        
        yraw=(x{II}(1,:)/scale_factor);
        %clip_count=length(find((yraw>=65536-1)|(yraw<65536)));  %Counts number of times signal pegs out.
        clip_count=length(find(abs(yraw-65536)<2));
        transient_params.clipped(Icall)=clip_count;
        
         
        format long e
        
        
        
        %best_ctimes(Icall)=cstart(II);
        transient_params.ctime(Icall)=cstart(II);
        transient_params.index(Icall)=Iref(II);
        
        
        
        %if ~isempty(findstr(run_options.short_fname,'AURAL'))
        if yes_decimate
            yy=decimate(x{II}(1,:),10,'FIR');
        else
            yy=x{II}(1,:);
        end
        
        yy=yy-mean(yy); %Smooths filtering
        Hd = dfilt.cascade(B_dfilt_lopass{Ifilt},B_dfilt_hipass{Ifilt2});
        yy=filter(Hd,yy);
        %yy=filter(B_dfilt_lopass{Ifilt},yy);
        %yy=filter(B_dfilt_hipass{Ifilt2},yy);
        
         
        %%WARNING!  UP TO JUNE 21 this was active--cut most of signal...
        %% and messed up dt measurement in get_level_metrics, because buffertime
        %%      not compensated...
       
        %Nf=round(length(Bfilt_lopass)/2);
        Nf=round(0.5*(Ncoef.lopass(Ifilt)+Ncoef.hipass(Ifilt2)));
        %yy=yy(Nf:(end-Nf));
        %The 'filter' command only destroys the first part of the signal...
        yy=yy(Nf:(end-Nf));
        bufferTime=bufferTime_org-Nf/Fs_filt;
        
        
        if length(yy)<round(Fs* bufferTime)
            disp('signal shorter than buffer time.');
            % keyboard
            continue;
        end
        
        try
            if debug_snips_check&run_options.debug<Icall  %check snip files
                ctime2str(cstart(II))
                features=get_level_metrics_simple(yy,Fs_filt, bufferTime,1);
                
                pause;
                
            else
                features=get_level_metrics_simple(yy,Fs_filt, bufferTime);
            end
            
            %%%Use below to compare bandpass filtered result with original result...
            
            %                         features=get_level_metrics(x{II}(1,:),Fs, bufferTime, bandwidth,transient_params.f_center,cstart(II));
            %                         subplot(3,1,1)
            %                         title(sprintf('Filter span: %i to %i Hz',filter_bandwidth(1),filter_bandwidth(2)));
            %                         set(gcf,'pos',[760   657   560   420]);
            %                         pause;close all
            
        catch
            disp(sprintf('Get_level_metrics_simple failure at Iref %i\n',Iref(II)));
            keyboard
            continue
            
        end
        
        transient_params.level.peak(Icall)=features.peak;
        %transient_params.level.peakF(Icall)=features.peakF;
        transient_params.level.t_Malme(Icall)=features.t_Malme;
        transient_params.level.SE_Malme(Icall)=features.SE_Malme;
        transient_params.level.rms_Malme(Icall)=features.rms_Malme;
       
        
        %%%Estimate noise levels preceeding call...
        transient_params.noise.rms(Icall)=features.noise.rms;
        transient_params.noise.SE(Icall)=features.noise.SE;
        transient_params.noise.duration(Icall)=features.noise.duration;
       
        transient_params.level.SNR(Icall)=(features.rms_Malme./features.noise.rms).^2 ;
        
        
    end %II, call mark, call
end  %I -- calls
close(hh)

end %entire function
