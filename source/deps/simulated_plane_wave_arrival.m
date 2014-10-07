%%%simulate_plane_wave_arrival.m%%%%%
% Simulate a signal arriving across a set of array elements

%function x=simulated_plane_wave_arrival(params,Iplot)
%Input:
%    params: structure with fields:
%           signal_type: 'impulse','FM'
%           c:  medium sound speed in m/s
%           Nchan: number of channels
%           angle: degrees, '0' is broadside
%           L:  array spacing in meters
%           Fs: sampling rate in Hz
%           total_window_time:  total sample duration in seconds
%           delay: relative delay of signal relative to center of signal,
%                   sec (must be positive number)
%           noise_var:  how much noise to add: variance relative to a peak
%               of '1'
%           FM: structure with FM parameters
%               fstart=100;
%               fend=300;
%               tduration=1;
%               tstart=1;
%               SNRdB=8;
%               SNRchc='PSD';
function [x,t]=simulated_plane_wave_arrival(params,Iplot)

signal_type=params.signal_type;  %'impulse','FM'
Nchan=params.Nchan;
angle=params.angle;
delay0=ceil(params.Fs*params.delay);
L=params.L;  %Array spacing in meters
%params.Fs=4000;
total_window_time=params.total_window_time;

switch signal_type
    case 'FM'
        fstart=params.FM.fstart;
        fend=params.FM.fend;
        tduration=params.FM.tduration;
        tstart=params.FM.tstart;
        SNRdB=params.FM.SNRdB;
        SNRchc=params.FM.SNRchc;
        
        disp('FM Sweep');
        [y,ysignal,ynoise,noise_var]=simulate_FMsweep_with_noise(params.Fs,fstart,fend, ...
            total_window_time,tduration,tstart,SNRdB,SNRchc);
    case 'impulse'
        total_window_time=256/params.Fs;
        
        Npt=round(total_window_time*params.Fs);
        y=zeros(Npt,1);
        y(Npt/2,1)=1;
        y(Npt/2-1,1)=0.5;
        y(Npt/2+1,1)=0.5;
        y(Npt/2-2,1)=0.25;
        y(Npt/2+2,1)=0.25;
        
        if isfield(params,'noise_var')
            noise_var=params.noise_var;
        else
            noise_var=0.001;
        end
        ysignal=y;
        
end


delay=(params.Fs*sin(angle*pi/180)*L/params.c);
delay0_offset=0;
if min(delay0)<0
   delay0_offset=abs(min(delay0))+1; 
end

%Reserve time space for delays and offsets
dn=ceil(max(delay)-min(delay));
dn0=ceil(max(delay0)-min(delay0));
x=ones(1+(Nchan+1)*dn+dn0+length(y),Nchan);
x=sqrt(noise_var)*randn(size(x)).*x;
for J=1:length(angle)
    for I=1:Nchan
        if delay(J)>0
            K=round(delay0_offset+delay0(J)+(I-1)*delay(J))+1+(1:length(y));
            
        else
            K=round(delay0_offset+delay0(J)+(Nchan-I)*abs(delay(J)))+1+(1:length(y));
            
            
            
        end
        x(K,I)=x(K,I)+ysignal;
        
    end
end

t=(1:length(x(:,1)))/params.Fs;
if Iplot>0
    figure
    
    Nplot=round(linspace(1,Nchan,10));
    for I=1:length(Nplot)
        subplot(length(Nplot),1,I)
        plot(t*1000,x(:,Nplot(I)))
    end
    xlabel('Time (msec)');
end

%
%end
end

