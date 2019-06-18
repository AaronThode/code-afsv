%%%%MultipleBandEnergyDetector.m
%%%  Given a time series x and a series of band parameters,
%%%  Detect transients over multiple bands
%  burntime: Time in minutes to build up average stats.
%   flo_det, fhi_det:  The absolute minimum and maximum frequencies
%   	to monitor for signals of interest,
%
%   eq_time: an equalization time in seconds.  If eq_time is zero, the equalization is not updated after the first 20 samples.
%
%   bandwidth:  The bandwidth of a subdetector--must be less than fhi_det-flo_det
%
%   threshold:  dB threshold SNR needed to exceed to start detection.
%
%   buffertime: How much time to stuff before and after detection start and stop when writing time series output.
%
%   TolTime: How much time in seconds must pass before new detection started?
%
%   MinTime: Minimum time in seconds a signal needs to exceed threshold SNR to register a detection.
%
%   debug: If 0, no subdetector output.  If 1, write subdetector SEL.  If 2, write suddetector equalization value.
%   	If 3, write ratio of current SEL to equalization SEL (SNR estimate).
%

%    param.eq_time='10';   param_desc{K}='Equalization time (s): should be roughly twice the duration of signal of interest';K=K+1;
%             param.bandwidth='.05';     param_desc{K}='Bandwidth of detector in kHz';K=K+1;
%             param.threshold='10';  param_desc{K}='Threshold in dB to accept a detection';K=K+1;
%             param.snips_chc='1';  param_desc{K}='0 for no snips file, 1 for snips file of one channel, 2 for snips file of all channels';K=K+1;
%             param.bufferTime='0.5'; param_desc{K}='Buffer Time in seconds to store before and after each detection snip, -1 suppress snips file';K=K+1;
%             param.TolTime='1e-4';  param_desc{K}='Minimum time in seconds that must elapse for two detections to be listed as separate';K=K+1;
%             param.MinTime='0';     param_desc{K}='Minimum time in seconds a required for a detection to be logged';K=K+1;
%             param.MaxTime='3';     param_desc{K}= 'Maximum time in seconds a detection is permitted to have';K=K+1;
%             param.debug='0';       param_desc{K}= '0: do not write out debug information. 1:  SEL output.  2:  equalized background noise. 3: SNR.';K=K+1;
%
function [detect,debug]=MultipleBandEnergyDetector(x,tabs_start,params)

[~,FF,TT,B] = spectrogram(x,hanning(params.Nfft),round(params.ovlap*params.Nfft),params.Nfft,params.Fs);
B=10*log10(B); %Convert to dB
threshold=params.threshold;
Ncol=size(B,2);

%%%Set up bands

flo=params.flo_det:(0.5*params.bandwidth):params.fhi_det;
fhi=flo+params.bandwidth;
Iflo=ceil(flo*params.Nfft/params.Fs);
%Iband=ceil(params.bandwidth*params.Nfft/params.Fs);
Ifhi=ceil(fhi*params.Nfft/params.Fs);

Igood=find(Ifhi<params.fhi_det*params.Nfft/params.Fs);
Ifhi=Ifhi(Igood);Iflo=Iflo(Igood);fhi=fhi(Igood);flo=flo(Igood);
Ndet=length(Ifhi);

Ispan=Iflo(1):Ifhi(end);
[~,Iburn]=min(abs(params.burn_in_time*60-TT));
eq=mean(B(Ispan,1:Iburn),2);
Ifhi=Ifhi-Iflo(1)+1;
Iflo=Iflo-Iflo(1)+1;

dT=TT(2)-TT(1);
dT=(1-params.ovlap)*params.Nfft/params.Fs;
Imin_time=round(params.MinTime/dT);
Itol_time=round(params.TolTime/dT);
Imax_time=round(params.MaxTime/dT);

%%Sneaky fast way
% II=[];
% for I=1:length(Iflo)
%      II=[II;Iflo(I):Ifhi(I)];
% end
%  band.eq=sum(eq(II),2);
%  debug.detect=sum(B(II,:));

%%Easy-to-understand way
for I=1:length(Iflo)
    band.eq(I)=10*log10(sum(10.^(eq(Iflo(I):Ifhi(I))/10)));
    debug.detect(I,:)=10*log10(sum(10.^(B(Iflo(I):Ifhi(I),:)/10),1));
end


dn = (1-params.ovlap)*params.Nfft;
alpha = 0.01^(dn/(params.eq_time*params.Fs));  %20 dB depression

%%Initialize detector.
count_estimate=round(10*length(x)/params.Fs);  %Assume around 10 detections/second
fieldnames={'tstart','tend','duration','amplitude','fmin','fmax'};
for I=1:length(fieldnames)
    detect.(fieldnames{I})=zeros(count_estimate,1);
end
count=0;
detection_status=-2*ones(Ndet,1);  %-2: off; -1: possible off 1: possible on  2: on
active_detectors=0;
tend=zeros(Ndet,1);
tstart=zeros(Ndet,1);
magnitude=tstart;
peak_index=tstart;

tstart_total=0;tend_total=0;
global_reset=false;
write_flag=false;

%%Cycle through detector.
%%  Written as for loop to make it easy
for I=((1+Iburn):Ncol)
    global_reset=false;
    val=debug.detect(:,I);  %Slice of spectrogram
    for J=1:Ndet
        if write_flag||global_reset %Detection has offically ended, cycle thorugh rest of detectors quickly
            continue
        end
        if val(J)>=band.eq(J)+threshold  %threshold exceeded
            switch detection_status(J)
                case -2 %OFF
                    detection_status(J)=1; %POSSIBLE ON
                    tstart(J)=I;
                    tend(J)=I;
                case 1 %POSSIBLEON
                if (I-tstart(J))>=Imin_time
                    detection_status(J)=2;
                    active_detectors=active_detectors+1;
                    magnitude(J)=val(J);
                    peak_index(J)=I;
                    if active_detectors==1
                        tstart_total=I;
                    end
                end
                
                case 2 %ON
                magnitude(J)=max([magnitude(J) val(J)]);
                if magnitude(J)==val(J)
                    peak_index(J)=val(J);
                end
                
                %%%Force reset?
                if I-tstart_total>=Imax_time
                    active_detectors=0;
                    tend(J)=I;
                    tend_total=tend(J);
                    write_flag=true;
                    if ~global_reset
                        for KK=0:Ndetectors
                            detection_status(KK)=-2;
                            band.eq(KK)=(1-0.25)*val(KK)+0.25*band.eq(KK);
                            debug.eq_history(J,I)=band.eq(J);
                            tend(KK)=I;
                            %writeMe(KK)=true;
                        end
                        global_reset=true;
                    end
                    
                    
                end %Imax_time
                
                case -1 %POSSIBLE OFF: we have just dipped below threshold for less than EndTolerance time, take back
                    detection_status(J)=2;
                    tend(J)=I;
            end
        else  %threshold not obtained
            
            switch detection_status(J)
                case -2 %OFF
                    band.eq(J)=(alpha).*band.eq(J)+(1-alpha)*val(J); %Update background estimate
                case 2 %ON
                    tend(J)=I;  %%Mark end of detection
                    detection_status(J)=-1; %Possible OFF
                case 1 %POSSIBLE ON:  We have briefly crossed threshold but have dipped back below
                    detection_status(J)=-2;
                    tstart(J)=0;
                    tend(J)=0;
                case -1 %POSSIBLE OFF
                    band.eq(J)=(alpha).*band.eq(J)+(1-alpha)*val(J); %Update background estimate
                    if tend(J)+Itol_time<=I&&(tend(J)~=0)  %%%if enough time has passed since last detection
                        detection_status(J)=-2;
                        active_detectors=active_detectors-1;
                        %writeMe(J)=true;
                        if active_detectors==0
                            if Imin_time+Itol_time<=(tend(J)-tstart_total)
                                tend_total=tend(J);
                                write_flag=true;
                            else
                                reset_detect;  %detection is too short
                            end
                        end
                    end
                    
            end  %detection_status
              
        end %if threshold
                
    end %J
    debug.eq_history(:,I)=band.eq;
     debug.detection_status(:,I)=detection_status; 
    %%%We've now worked through all bands.  Create or close
    %%%  an official detection.
    
    if write_flag % detection completed
        %%%Is it too short?
            count=count+1;
            detect.tstart(count)=tstart_total;
            detect.tend(count)=tend_total;
            detect.fmin(count)=min(flo(tend~=0));
            detect.fmax(count)=max(fhi(tend~=0));
            detect.amplitude(count)=max(magnitude);
            %detection_status(:)=false;
            write_flag=false;
            reset_detect;
            
    
    end %write flag
    
    
end  %I

for I=1:length(fieldnames)
    detect.(fieldnames{I})=detect.(fieldnames{I})(1:count);
end

if params.debug
    %close all
    figure(100);hold off
    subplot(4,1,1);hold off
    imagesc([],0.5*(flo+fhi),(debug.eq_history));axis('xy'); colorbar('west')
    caxis([min(min(debug.detect)) max(max(debug.detect))]);
    subplot(4,1,2);hold off
    imagesc([],0.5*(flo+fhi),(debug.detect));axis('xy'); colorbar('west')
    caxis([min(min(debug.detect)) max(max(debug.detect))]);
    subplot(4,1,3);hold off
    imagesc([],0.5*(flo+fhi),(debug.detect-debug.eq_history));axis('xy'); colorbar('west')
    caxis([params.threshold 30]);
    subplot(4,1,4);hold off
    imagesc([],0.5*(flo+fhi),debug.detection_status);axis('xy'); colorbar('west')
   
    
    hold on
    for I=1:count
       line([detect.tstart(I) detect.tend(I)],detect.fmax(I)*[1 1],'color','w','linewidth',3); 
       line(mean([detect.tstart(I) detect.tend(I)])*[1 1],[detect.fmin(I) detect.fmax(I)],'color','w','linewidth',3); 
    end
    pause;
   
    
end

%%%Convert times to absolute times
detect.duration=detect.tend-detect.tstart;
detect.tstart=tabs_start+datenum(0,0,0,0,0,detect.tstart);
detect.tend=tabs_start+datenum(0,0,0,0,0,detect.tend);


    function reset_detect
        active_detectors=0;
        tend_total=0;
        tstart_total=0;
        write_flag=false;
        for KKK=1:Ndet
            magnitude(KKK)=0;
            peak_index(KKK)=0;
            tstart(KKK)=0;
            tend(KKK)=0;
            detection_status(KKK)=-2;
            %writeMe(KKK)=false;
            
        end
    end
end



