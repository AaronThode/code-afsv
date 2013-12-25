%%%%%%extract_transient_levels.m%%%%%%%%%%%%%%%%%%
% Routine called by process_one_unit.m to extract level information from
% clips in *.snips file.  
%  Aaron Thode
%  September, 2008
%
% Inputs:
%       Igood: indicies of data_all that contain desired signals.
%           Typically contains indicies of signals that have passed ICI test.
%       data_all: Output of a *detsum file.  Contains fields
%               .ctime: vector of ctimes of detections...
%               .npt:  row vector of duration of signals in samples.
%           Note: data_all contains all data from *.snips file, and Igood
%               is used to select subset.
%       Fs: sampling rate of raw data.
%       bufferTime(s): typically value of param.energy.bufferTime
%               Note that this value may be adjusted if FIR filtering takes place
%       bandwidth:  For use with estimating noise PSD levels
%       filter_bandwidth: filter passband over which to estimate pulse duration...
%       Nfft,ovlap, used for estimating noise PSD levels
%       run_options:
%               .debug.image_signals:
%               datatype: string, either 'short' or 'double', depending on
%                       data format of snips file...
%               Ncalls:number of signals to load into RAM memory for
%               processing.
%               short_fname: used to ID snips file
%               calibration_keyword:  Used by calibrate_GSI_signal.m to
%                   select calibration scheme for raw data...
%       RawFileName:  Full pathname of raw data file, used to check for
%           clipping
%       dir_out: filter passband over which to estimate pulse duration...
%
%function [transients]=extract_transient_levels(Igood,data_all,Fs,bufferTime,bandwidth,  ...
%        Nfft,ovlap,filter_bandwidth,run_options,RawFileName,dir_out)
%
function [transient_params]=extract_transient_levels(Igood,data_all,Fs,bufferTime,bandwidth,Nfft,ovlap,filter_bandwidth,run_options,RawFileName,dir_out)
persistent Bfilt B_dfilt Imaxx

transient_shots=[];
bufferTime_org=bufferTime;

if isempty(Igood)
    return
end
%max_clip_count=run_options.max_clip_count;  %Number of clips allowed before data marked as tainted..

debug_snips_check=run_options.debug.snips;

transient_params.f_center=[10 20 32 40 50 63 80 100 125 160 200 250 315 395]; %1/3 octave levels
f_band=500;  %Hz
transient_params.f_upper=2000:f_band:5000;
transient_params.f_upper(end)=5000;

%f_low = round(0.707* f_center);
%f_high = round(1.414* f_center);


Icall=0;
best_ctimes=zeros(1,length(Igood));
Bmean=[];
% Option 1:

if length(Igood)<run_options.Ncalls_to_sample
    run_options.Ncalls_to_sample=length(Igood);
end

%scale_factor=get_calibration_scale_factor(run_options.calibration_keyword);
scale_factor=1;

if ~isempty(findstr(run_options.short_fname,'AURAL'))
    Fs_filt=Fs/10;
else
    Fs_filt=Fs;
end
%%Optional prefiltering
%param.airgun.filter_bandwidth
if isempty(Bfilt)&&~isempty(filter_bandwidth)
    %Bfilter = brickwall_bpf(filter_bandwidth,Fs,0);
    
    %Flag feature for max frequency
    for Ifea=1:length(data_all.names)
       if strcmp(data_all.names{Ifea},'max_freq')
          Imaxx=Ifea; 
       end
    end
    for Ifilter=1:length(transient_params.f_upper)
        
        filter_bandwidth(3)=transient_params.f_upper(Ifilter);
        filter_bandwidth(4)=filter_bandwidth(3)+filter_bandwidth(2)-filter_bandwidth(1);
        disp(sprintf('Building %.0f-%.0f Hz filter',filter_bandwidth(2),filter_bandwidth(3)));
        disp(sprintf('Transition zone bandwidth: %6.2f',filter_bandwidth(2)-filter_bandwidth(1)));
        
        [B_dfilt{Ifilter},Bfilt]=brickwall_bpf(filter_bandwidth,Fs_filt,0);
    end
    disp('Finished building');
else
    disp('Filter already built');
end
tt1=tic;
%disp(sprintf('missing %i calls at end',floor(length(I
fclose('all');


for I=1:ceil(length(Igood)/run_options.Ncalls_to_sample)
    disp(sprintf('Batch %i of %i',I,ceil(length(Igood)/run_options.Ncalls_to_sample)));
    Iabs=run_options.Ncalls_to_sample*(I-1)+(1:run_options.Ncalls_to_sample);
    Iabs=Iabs(Iabs<=length(Igood));
    Iref=Igood(Iabs);
    snips_name=dir([dir_out '/' run_options.short_fname '*.snips']);
    
    
    %Load data from JAVA program.  Note that a 4th-order Butterworth filter has already been applied, but
    %clipping events have been retained...
    x=readEnergySnips([dir_out '/' snips_name.name], Iref,'double','cell','keep_open');
    %toc(tt1)
    
    
    % disp(sprintf('%6.2f percent done',100*I/(length(Igood)/run_options.Ncalls_to_sample)));
    cstart=data_all.ctime(Iref);
    dstart=diff(cstart);dstart=[0 dstart cstart(end)];
    
    for II=1:length(Iref)
        Icall=Icall+1;
        
        
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
        
        %         if clip_count>0
        %             subplot(2,1,1)
        %             plot(yraw);title(num2str(clip_count));
        %             subplot(2,1,2)
        %             spectrogram(yraw,hanning(Nfft),ovlap,Nfft,Fs,'yaxis')
        %             pause
        %         end
        format long e
        
        
        if debug_snips_check==1  %check snip files
            plot_images;
        end
        %best_ctimes(Icall)=cstart(II);
        transient_params.ctime(Icall)=cstart(II);
        transient_params.index(Icall)=Iref(II);
        
        %         if abs(datenum(1970,1,1,0,0,cstart(II))-7.333041599074074e+05)<datenum(0,0,0,0,0,10)
        %             keyboard
        %         end
        
        if ~isempty(findstr(run_options.short_fname,'AURAL'))
            yy=decimate(x{II}(1,:),10,'FIR');
        else
            yy=x{II}(1,:);
        end
        
        if ~isempty(Bfilt)
            
            %yy=filter(B_dfilt,[x{II}(1,:) zeros(1,length(Bfilt))]);
            max_freq=data_all.features(Imaxx,Iref(II));
            transient_params.level.max_freq(Icall)=max_freq;
            Ifilt=ceil(max_freq/f_band);
            yy=filter(B_dfilt{Ifilt},yy);
            
            
            %yy=filtfilt(B_dfilt,yy);
            
            
            %%WARNING!  UP TO JUNE 21 this was active--cut most of signal...
            %% and messed up dt measurement in get_level_metrics, because buffertime
            %%      not compensated...
            Nf=round(length(Bfilt)/2);
            %yy=yy(Nf:(end-Nf));
            yy=yy(Nf:(end-Nf));
            bufferTime=bufferTime_org-Nf/Fs_filt;
        
            
            %%An exercise to compare different kinds of fitering
%             figure(100);
%             subplot(2,1,1);
%             plot(x{II}(1,:));
%             subplot(2,1,2);
%             plot(filter(B_dfilt{Ifilt},x{II}(1,:)));
%             subplot(3,1,3);
%             plot(filtfilt(B_dfilt{Ifilt},x{II}(1,:)));
            
        end
        
        if length(yy)<round(Fs* bufferTime)
            disp('signal shorter than buffer time.');
            % keyboard
            continue;
        end
        
        try
            %features=get_level_metrics(yy,Fs_filt, bufferTime, bandwidth,transient_params.f_center,cstart(II)); %Debug version
            features=get_level_metrics(yy,Fs_filt, bufferTime, bandwidth,transient_params.f_center);
            
            %%%Use below to compare bandpass filtered result with original result...
            
%                         features=get_level_metrics(x{II}(1,:),Fs, bufferTime, bandwidth,transient_params.f_center,cstart(II));
%                         subplot(3,1,1)
%                         title(sprintf('Filter span: %i to %i Hz',filter_bandwidth(1),filter_bandwidth(2)));
%                         set(gcf,'pos',[760   657   560   420]);
%                         pause;close all
            
        catch
            disp(sprintf('Get_level_metrics failure at Iref %i\n',Iref(II)));
             keyboard
            continue
           
        end
        
        transient_params.level.peak(Icall)=features.peak;
        transient_params.level.peakF(Icall)=features.peakF;
        transient_params.level.t_Malme(Icall)=features.t_Malme;
        transient_params.level.SEL_Malme(Icall)=features.SEL_Malme;
        transient_params.level.rms_Malme(Icall)=features.rms_Malme;
        
        %transient_params.level.t20dB(Icall)=features.t20dB;
        %transient_params.level.SEL20dB(Icall)=features.SEL20dB;
        transient_params.level.SEL_FFT(Icall)=features.SEL_FFT;
        transient_params.level.SEL_FFT_band(1:length(features.SEL_FFT_band),Icall)=features.SEL_FFT_band;
        transient_params.level.rms_FFT_band(1:length(features.rms_FFT_band),Icall)=features.rms_FFT_band;
        transient_params.freq_bandwidth=features.freq_bandwidth;
        
        Nbands=length(features.SEL_FFT_third_octave);
        if Nbands>0
            transient_params.level.SEL_FFT_third_octave(1:Nbands,Icall)=features.SEL_FFT_third_octave;
            transient_params.level.rms_FFT_third_octave(1:Nbands,Icall)=features.rms_FFT_third_octave;
        end
         
        %%%Estimate noise levels preceeding call...
        transient_params.noise.rms(Icall)=features.noise.rms;
        transient_params.noise.SEL(Icall)=features.noise.SEL;
        transient_params.noise.duration(Icall)=features.noise.duration;
        transient_params.noise.SEL_FFT_band(1:length(features.noise.SEL_FFT_band),Icall)=features.noise.SEL_FFT_band;
        transient_params.noise.rms_FFT_band(1:length(features.noise.rms_FFT_band),Icall)=features.noise.rms_FFT_band;
        transient_params.noise.freq_bandwidth=features.noise.freq_bandwidth;
        
        %%%Compute SNR
        % %Must square rms to get intensity units.
        % Subtract one in order to remove noise from signal if rms_Malme does not already subtract noise
        % estimate.
        %transient_params.level.SNR(Icall)=(features.rms_Malme./features.noise.rms).^2 -1;  
        transient_params.level.SNR(Icall)=(features.rms_Malme./features.noise.rms).^2 ;  
        
        
    end %II, call mark, call
end  %I -- calls
%Closes file...
%[x,nstarts_snips]=readEnergySnips([dir_out '/' snips_name.name], Igood(Iabs),'short','cell');

    function plot_images  %May need to be debugged..
        %[S,F,T,PP1]=spectrogram(y,128,96,128,head.Fs,'yaxis');
        %Debug code: test IIR filter...
        %tlen=2* bufferTime+data_all.npt(1,Iref(II))/Fs;
        %[y{1},t,head]=readsiof(RawFileName,cstart(II)-2* bufferTime,tlen);
        %y{1}=y{1}-mean(y{1});
        
        %y{2}=calibrate_GSI_signal(y{1},'filter');
        %iin=round(Fs* bufferTime)+data_all.npt(1,Iref(II));
        %y{3}=x{II}-mean(x{II});
        y{1}=xorg{II}(1:(length(xorg{II})-round( bufferTime*Fs)));  %snips file
        %y{2}=y0((end-iin):end)-(2^16)/2;
        %y{2}=y0((end-iin):end);  %Direct Siofile download...
        y{3}=x{II};  %MATLAB filtered...
        if ~isempty(y{1}),
            for J=1:3,
                
                %if J<3,y{J}=y{J}((end-iin):end);end
                [S,F,T,PP{J}]=spectrogram(y{J},128,96,128,head.Fs,'yaxis');
                
            end
            
            for J=1:3,
                figure(4);
                
                subplot(3,1,J);
                imagesc(T,F,10*log10(abs(PP{J})));axis('xy');
                %caxis([-20 50]);
                colorbar;
                switch J
                    case 2
                        caxis([0 70])
                        title(sprintf('Raw Data: %s time to previous %6.2f time to next: %6.2f',datestr(datenum(1970,1,1,0,0,cstart(II))),dstart(II),dstart(II+1)));
                    case 3
                        title(sprintf('MATLAB filtered data, Date: %s time to previous %6.2f time to next: %6.2f',datestr(datenum(1970,1,1,0,0,cstart(II))),dstart(II),dstart(II+1)));
                    case 1
                        title(sprintf('filtered snips index is %i',Iref(II)));
                        caxis([0 70]);
                end
                
                figure(2)
                subplot(3,1,J)
                plot(y{J});grid on;
                
            end
            %keyboard
            yes=input('Continue? (-1 to stop):');
            if yes==-1,
                debug_snips_check=0;
            end
        end
    end  % inner function

end %entire function
