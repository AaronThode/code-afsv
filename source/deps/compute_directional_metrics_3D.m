%%%%compute_directional_metrics_3D.m%%%%%%%%%%%%%%%%%%%
%function [TT,FF,output_array,PdB,param, Ix,Iy,Iz]=compute_directional_metrics(x,metric_type, ...
%   Fs,Nfft, ovlap, param,filetype,reactive_flag)
% x:  array of vector data, each row is a separate channel
% metric_type{Iwant}:'Directionality','TransportVelocity','KEtoPERatio' ,'IntensityPhase
%           If a cell array output_array will be a cell array
% filetype: 'DIFAR' or 'gsi'
% reactive_flag:  if true, compute reactive intensity.  Vector same size as
%           metric_type
% Fs: sampling rate in Hz
% Nfft: FFT size used to compute frequency bin size
% ovlap: fraction ovlap when computing spectrograms
% paramaters
%     param.sec_avg=(Batch_vars.sec_avg);
%     param.climm=(Batch_vars.climm);
%     param.brefa=(Batch_vars.brefa);
%     param.alg=Batch_vars.alg;
%     param.phase_calibration='Arctic5G_2014';
%
% Output:
%  TT,FF, output_array{Nwant}:  TT vector of times and FF vector of Hz for
%               output_array grid.
%  PdB: power spectral density of pressure autospectrum
%  param:  altered parameters of input param
%  Ix, Iy:  x and y active intensity

function [TT,FF,output_array,PdB,param,Ix,Iy,Iz]=compute_directional_metrics_3D(x,metric_type, ...
    Fs,Nfft, ovlap, param,filetype,reactive_flag)

use_wavelet=any(contains(lower(metric_type),'wavelet'));
if use_wavelet
    Iend=strfind(metric_type,'wavelet')-1;
    metric_type=metric_type(1:Iend);
end

TT=[];FF=[];output_array=[];PdB=[];
if ~iscell(metric_type)
    temp=metric_type;clear metric_type
    metric_type{1}=temp;
end

if length(metric_type)~=length(reactive_flag)
    disp('Metric type not same length as reactive flag');
    return
end

Ix=[];Iy=[];
if strcmpi(filetype,'PSD')
    return
end

if ~exist('reactive_flag','var')
    reactive_flag=false;
end

if ischar(param.brefa)
    param.brefa=str2num(param.brefa);
end
%%%%Two files that can get active intensity directionality: *.gsi and
%%%%*DIFAAR*.wav files.


if strcmpi(filetype,'gsi')
    
    if length(x(1,:))<Nfft/2
        disp('Signal sample shorter than Nfft');
        return
    end
    
    %%%%Test alternate, faster version, used since July 29, 2018
    % tic
    
    if size(x,2)>size(x,1)
        x=x';
    end
    
    if use_wavelet
        M=size(x,1);
        dn=1;
        fb=cwtfilterbank('SamplingFrequency',Fs,'SignalLength',size(x,1),'FrequencyLimits',[10 0.95*Fs/2],'VoicesPerOctave',12);
        for J=1:3
            [B(J,:,:),FF,COI] = wt(fb,x(:,J));  %use cwtfilterbank for future speed
        end
        TT=(1:M)/Fs;
        Nf=length(FF);
    else
        dn=round((1-ovlap)*Nfft);
        M=floor(1+(size(x,1)-Nfft)/dn);
        B=zeros(size(x,2),Nfft/2+1,M);
        
        for J=1:4
            fprintf('Spectrogram channel %i\n',J);
            [B(J,:,:),FF,TT] = spectrogram(x(:,J),Nfft,round(ovlap*Nfft),Nfft,Fs);
        end
          Nf=Nfft/2+1;  %Should be the same as length(FF)
        clear x
    end
    
     %rho=1000;c=1500;
    
    Ix=squeeze(((B(1,:,:).*conj(B(2,:,:)))));
    Iy=squeeze(((B(1,:,:).*conj(B(3,:,:)))));
    Iz=squeeze(((B(1,:,:).*conj(B(4,:,:)))));
    
    if ~all(contains(metric_type,'Directionality'))
        pressure_autospectrum=0.5*squeeze(abs(B(1,:,:)).^2)./(param.rho*param.c0.^2);
        normalized_velocity_autospectrum=0.5*param.rho.*squeeze(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2+abs(B(4,:,:)).^2);
        energy_density=0.5*abs(normalized_velocity_autospectrum+pressure_autospectrum);
    end
    %toc
    get_newparams=false;
    clear B
    
elseif contains(filetype,'DIFAR')
    M=floor(1+(max(size(x))-Nfft)/dn);
    [Ix,Iy,TT,FF,PdB]=demultiplex_DIFAR(x,Fs,Nfft,ovlap);
    Nf=length(FF);
    get_newparams=false;
    
    pressure_autospectrum=10.^(PdB/10);
    if ~isfield(param,'brefa')
        get_newparams=true;
    end
    
end  %if gsi or DIFAR

%sec_avg=input('Enter time to average over (sec; 0 does no averaging):');
if ~isfield(param,'sec_avg')
    get_newparams=true;
    
end

if get_newparams
    Batch_vars=get_Azigram_Callback(param);
    param.sec_avg=(Batch_vars.sec_avg);
    param.climm=(Batch_vars.climm);
    if ~strcmpi(filetype,'gsi')
        param.brefa=(Batch_vars.brefa);
        %%%Don't alter hdr.brefa, which contains the correction.
    else
        %param.brefa=(Batch_vars.brefa);
        param.brefa=(param.brefa);
        
    end
    param.alg=Batch_vars.alg;
end

% Batch_vars.sec_avg	=	'0.1';	Batch_desc{1}	=	'Seconds to average PSD for long-term display, if "0" no averaging' ;
% Batch_vars.climm='[0 360]'; Batch_desc{2}='Bearing Range Color Scale';
% Batch_vars.brefa='11.7';Batch_desc{3}='Bearing bias/correction (default is 11.7 degrees)';
% Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, 'Vector Sensor Processing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Average Ix and Iy, if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(param.sec_avg)
    sec_avg=str2double(param.sec_avg);
else
    sec_avg=param.sec_avg;
end

if ~isempty(sec_avg)&&sec_avg>0
    
    Navg=min([M floor(sec_avg*Fs/dn)]);  %Spectrogram samples (columns) per avg
    %fprintf('%i averages per sample.\n',Navg);
    if Navg==0
        Navg=1;
        sec_avg=dn/Fs;
    end
    Nsnap=max([1 floor((M-Navg)/((1-ovlap)*Navg))]);  %Number of averaged samples per window.
    Ix_avg=zeros(Nf,Nsnap);
    Iy_avg=Ix_avg;
    NV_avg=Ix_avg;
    PA_avg=Ix_avg;
    TT_avg=zeros(1,Nsnap);
    for J=1:Nsnap
        index=floor((J-1)*(Navg*(1-ovlap)))+(1:Navg);
        Ix_avg(:,J)=mean(Ix(:,index),2);
        Iy_avg(:,J)=mean(Iy(:,index),2);
        NV_avg(:,J)=mean(normalized_velocity_autospectrum(:,index),2);
        PA_avg(:,J)=mean(pressure_autospectrum(:,index),2);
        TT_avg(J)=mean(TT(index));
    end
    Ix=Ix_avg;
    Iy=Iy_avg;
    
    
    TT=TT_avg;
    energy_density=0.5*abs(NV_avg+PA_avg);
    normalized_velocity_autospectrum=NV_avg;
    pressure_autospectrum=PA_avg;
end  %sec_avg

if ~all(contains(metric_type,'Directionality'))
    if ~use_wavelet
        PdB=4+10*log10(2*pressure_autospectrum./(Nfft*Fs));  %%Power spectral density output
    else
        PdB=sqrt(pressure_autospectrum);
    end
end
%[~,FF,TT,PdB1] = spectrogram(x(:,1),Nfft,round(ovlap*Nfft),Nfft,Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Correct for phase misalignment in Ix and Iy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(param,'phase_calibration')
    param.phase_calibration='none';
end
[Ix,Iy]=correct_phase(Ix,Iy,FF,param.phase_calibration);


for J=1:length(metric_type)  %%for each request
    
    switch metric_type{J}
        case 'Azimuth'
            if ~reactive_flag(J)
                mu = single(atan2d(real(Ix),real(Iy)));
            else
                mu = single(atan2d(imag(Ix),imag(Iy)));
            end
            output_array{J}=bnorm(param.brefa+mu);
            %output_array{J}=bnorm(mu);
        case 'Elevation'
            if ~reactive_flag(J)
                mu = single(atan2d(real(Iz),sqrt((real(Iy)).^2+(real(Ix)).^2)));
            else
                mu = single(atan2d(imag(Iz),sqrt((real(Iy)).^2+(real(Ix)).^2)));
            end
            output_array{J}=bnorm(param.brefa+mu);
            %output_array{J}=bnorm(mu);
        case 'TransportVelocity'
            if ~reactive_flag(J)
                intensity=sqrt((real(Ix)).^2+(real(Iy)).^2+(real(Iz)).^2);
            else
                intensity=sqrt((imag(Ix)).^2+(imag(Iy)).^2+(imag(Iz)).^2);
            end
            
            %%%Effective velocity j/Sp
            output_array{J}=intensity./energy_density;
        case 'KEtoPERatio'
            output_array{J}=normalized_velocity_autospectrum./pressure_autospectrum;
        case 'IntensityPhase'
            if ~reactive_flag(J)
                output_array{J}=atan2d(sqrt(imag(Ix).^2+imag(Iy).^2),sqrt((real(Ix)).^2+(real(Iy)).^2));
                
                %output_array{J}=atan2d((imag(Ix)),((real(Ix)))); %Ix phase
                % output_array{J}=atan2d((imag(Iy)),((real(Iy)))); %Iy phase
            else
                output_array{J}=atan2d(sqrt(real(Ix).^2+real(Iy).^2),sqrt((imag(Ix)).^2+(imag(Iy)).^2));
            end
            
        case 'IntensityPhaseZ'
            if ~reactive_flag(J)
                output_array{J}=atan2d(imag(Iz),real(Iz));
                
                %output_array{J}=atan2d((imag(Ix)),((real(Ix)))); %Ix phase
                % output_array{J}=atan2d((imag(Iy)),((real(Iy)))); %Iy phase
            else
                output_array{J}=atan2d(sqrt(real(Ix).^2+real(Iy).^2),sqrt((imag(Ix)).^2+(imag(Iy)).^2));
            end
        case 'PhaseSpeed'
            %factor=2*(1500.^2)./(1000*1500.^2);
            factor=2*param.c0.^2;
            if ~reactive_flag(J)
                output_array{J}=factor*pressure_autospectrum./sqrt((real(Ix)).^2+(real(Iy)).^2+(real(Iz)).^2);
                
            else
                output_array{J}=factor*pressure_autospectrum./sqrt((imag(Ix)).^2+(imag(Iy)).^2+(imag(Iy)).^2);
                
            end
    end
end

end

function [Ix,Iy]=correct_phase(Ix,Iy,FF,phase_calibration_chc)
phase_calibration_chc='Arctic5G_2014';
%phase_calibration_chc='none';

Np=size(Ix,2);
switch(phase_calibration_chc)
    case 'none'
        return
    case 'Arctic5G_2014'
        %slopee=4.10e-04;
        %slopee=1.25*4.1e-04;  %%%radians phase per radians frequency, original bulk run
       
        
        %%%% slopee measured at 00:06:30 local time on 17 August, 2014 5G
        %%%%23:57:28  1 seconds 10/1/2014 5G
        %%% Also checked 04-Oct-2010 02:45:47.4 DASAR 5G, within 10° below
        %%%     375 Hz.  This signal is lower received level.
        %slopee=1.1*4.1e-04;  %%%radians phase per radians frequency, for 
        slopee=1.2*4.1e-04;  %%%radians phase per radians frequency, for 1 sec averages 
       
        
        %Phasee=10*pi/180+slopee*2*pi*(FF-75);  %Phase is 10 degrees at 75 Hz, linear
        Phasee=slopee*2*pi*FF;  %Nearly identical to above (y-intercept is
        %           -1 degrees in equation above.
        
        
        Phasee=exp(-1i*Phasee*ones(1,Np));
        
        Ix=Ix.*Phasee;
        Iy=Iy.*Phasee;
        
end

end


