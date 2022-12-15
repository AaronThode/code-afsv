%%%%compute_directional_metrics.m%%%%%%%%%%%%%%%%%%%
%function [TT,FF,output_array,PdB,param, Ix,Iy]=compute_directional_metrics(x,metric_type, ...
%   Fs,Nfft, ovlap, param,filetype,reactive_flag)
% x:  array of vector data, each row is a separate channel
% metric_type{Iwant}:'Directionality','ItoERatio','KEtoPERatio',
%       'IntensityPhase', 'Polarization','AdditiveBeamforming'
%           metric_type can also be a string.
%           If a cell array output_array will be a cell array
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
%     param.instrument = 'drifterM35' or 'DIFAR'
%
% Output:
%  TT,FF, output_array{Nwant}:  TT vector of times and FF vector of Hz for
%               output_array grid.
%  PdB: power spectral density of pressure autospectrum
%  param:  altered parameters of input param
%  Ix, Iy:  x and y active intensity

function [TT,FF,output_array,PdB,param,Ix,Iy]=compute_directional_metrics(x,metric_type, ...
    Fs,Nfft, ovlap, param,reactive_flag)


use_wavelet=any(contains(lower(metric_type),'wavelet'));
if ~isfield(param,'instrument')
    param.instrument='DASAR';
end

if use_wavelet
    Iend=strfind(metric_type,'wavelet')-1;
    metric_type=metric_type(1:Iend);
end

TT=[];FF=[];output_array=[];PdB=[];Ix=[];Iy=[];
if ~iscell(metric_type)
    temp=metric_type;clear metric_type
    metric_type{1}=temp;
end

if length(metric_type)~=length(reactive_flag)
    disp('Metric type not same length as reactive flag');
    return
end

if ~exist('reactive_flag','var')
    reactive_flag=false;
end

if ischar(param.brefa)
    param.brefa=str2num(param.brefa);
end


%%%%Two files that can get active intensity directionality: *.gsi and
%%%%*DIFAR*.wav files.


if contains(param.instrument,'DIFAR')
    
    M=floor(1+(max(size(x))-Nfft)/dn);
    [Ix,Iy,TT,FF,PdB]=demultiplex_DIFAR(x,Fs,Nfft,ovlap);
    Nf=length(FF);
    get_newparams=false;
    
    pressure_autospectrum=10.^(PdB/10);
    if ~isfield(param,'brefa')
        get_newparams=true;
    end
else  %%All other data is coming on on other channels.
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
        
        for J=1:3
            fprintf('Spectrogram channel %i\n',J);
            [B(J,:,:),FF,TT] = spectrogram(x(:,J),Nfft,round(ovlap*Nfft),Nfft,Fs);
        end
          Nf=Nfft/2+1;  %Should be the same as length(FF)
        clear x
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Correct gain as needed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gains=correct_gain(FF,param.instrument);
    %temp=10*log10(abs(squeeze(B(2,:,:))));
    %myfig=gcf;
    %figure; for III=1:3,subplot(3,1,III);imagesc(TT,FF,temp);ylim([0 2000]);colorbar;end
    %figure(myfig);
    
    for J=1:3
        B(J,:,:)=squeeze(B(J,:,:)).*Gains(:,J);
    end
    
    %%% If AdditiveBeamforming, conduct here and leave...
    if any(strcmpi(metric_type,'AdditiveBeamforming'))
        [B(2,:,:),B(3,:,:)]=correct_phase(squeeze(conj(B(2,:,:))),squeeze(conj(B(3,:,:))),FF,param.instrument);
        B(2:3,:,:)=conj(B(2:3,:,:));
        
        %%%%Additive beamforming
        thta=param.thta-param.brefa;  %%Convert to local reference
        Iout=find(strcmpi(metric_type,'AdditiveBeamforming'));
        output_array{Iout}=squeeze(B(1,:,:)+sind(thta)*B(2,:,:)+cosd(thta)*B(3,:,:));
        
        
        %%%Debug options, including |v|/p
        %output_array{Iout}=squeeze((sqrt(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2))./abs(B(1,:,:)));
        %x0=x(:,1)+sind(thta)*x(:,2)+cosd(thta)*x(:,3);
        return
    end
    
   % temp2=10*log10(abs(squeeze(B(2,:,:))));
    %rho=1000;c=1500;
    
    %myfig=gcf;
    %figure; for III=1:3,subplot(3,1,III);imagesc(TT,FF,temp2-temp);ylim([0 2000]);colorbar;end
    %figure(myfig);
    
   
    Ix=squeeze(((B(1,:,:).*conj(B(2,:,:)))));
    Iy=squeeze(((B(1,:,:).*conj(B(3,:,:)))));
    
    %if ~all(contains(metric_type,'Directionality'))  %%%Attempt to save
    %time and memory, but need this if trying to use transparency
    pressure_autospectrum=squeeze(abs(B(1,:,:)).^2);
    normalized_velocity_autospectrum=squeeze(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2);
    energy_density=0.5*abs(normalized_velocity_autospectrum+pressure_autospectrum);
    polarization=(real(B(2,:,:)).*imag(B(3,:,:))-real(B(3,:,:)).*imag(B(2,:,:)));
    polarization=-2*squeeze(polarization./(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2));
    %end
    %toc
    get_newparams=false;
    clear B    
    
end  %if gsi or DIFAR

%sec_avg=input('Enter time to average over (sec; 0 does no averaging):');
if ~isfield(param,'sec_avg')
    get_newparams=true;
    
end

if get_newparams
    Batch_vars=get_Azigram_Callback(param);
    param.sec_avg=(Batch_vars.sec_avg);
    param.climm=(Batch_vars.climm);
    if ischar(Batch_vars.brefa)
        param.brefa=eval(Batch_vars.brefa);
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

%if ~all(contains(metric_type,'Directionality'))  %%An attempt to be more efficient, not worth is
    if ~use_wavelet
        PdB=4+10*log10(2*pressure_autospectrum./(Nfft*Fs));  %%Power spectral density output
    else
        PdB=sqrt(pressure_autospectrum);
    end
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Correct for phase misalignment in Ix and Iy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ix,Iy]=correct_phase(Ix,Iy,FF,param.instrument);


for J=1:length(metric_type)  %%for each request
    
    switch metric_type{J}
        case 'Directionality'
            if ~reactive_flag(J)
                mu = single(atan2d(real(Ix),real(Iy)));  %Compass convention (usually atan2d(y,x))
            else
                mu = single(atan2d(imag(Ix),imag(Iy)));
            end
            output_array{J}=bnorm((param.brefa)+mu);
            %output_array{J}=bnorm(mu);
        case 'ItoERatio'
            if ~reactive_flag(J)
                intensity=sqrt((real(Ix)).^2+(real(Iy)).^2);
            else
                intensity=sqrt((imag(Ix)).^2+(imag(Iy)).^2);
            end
            
            %%%Effective velocity j/Sp
            output_array{J}=intensity./energy_density;
        case 'KEtoPERatio'
            output_array{J}=normalized_velocity_autospectrum./pressure_autospectrum;
        case 'Polarization'
           % output_array{J}=(real(Ix).*imag(Iy)-real(Iy).*imag(Ix))./pressure_autospectrum;
            
            output_array{J}=polarization;
        case 'IntensityPhase'
            if ~reactive_flag(J)
                output_array{J}=atan2d(sqrt(imag(Ix).^2+imag(Iy).^2),sqrt((real(Ix)).^2+(real(Iy)).^2));
                
                %output_array{J}=atan2d((imag(Ix)),((real(Ix)))); %Ix phase
                % output_array{J}=atan2d((imag(Iy)),((real(Iy)))); %Iy phase
            else
                output_array{J}=atan2d(sqrt(real(Ix).^2+real(Iy).^2),sqrt((imag(Ix)).^2+(imag(Iy)).^2));
            end
            
        case 'PhaseSpeed'
            if ~reactive_flag(J)
                output_array{J}=1450*pressure_autospectrum./sqrt((real(Ix)).^2+(real(Iy)).^2);
                
            else
                output_array{J}=1450*pressure_autospectrum./sqrt((imag(Ix)).^2+(imag(Iy)).^2);
                
            end
    end
end

end

%%%%%%%SubFunctions%%%%%%%%%%%%%

function Gains=correct_gain(FF,instrument_type)
switch instrument_type
    case 'DASAR'  %%%Remember energy is arriving as normal modes so a slight vertical offset.
        Gains=ones(length(FF),3);
        %Gains(:,2:3)=Gains(:,2:3)*cosd(28)*exp(1i*2*pi*(-2)/180); %150 Hz bowhead
        %Gains(:,2:3)=Gains(:,2:3)*cosd(20)*exp(1i*2*pi*(-4)/180); %250 Hz bowhead
       
        
        %%%%Thought normal modes might explain this, but better using a
        %%%%frequency-independent factor.
       D=18; c=1500;
       k=2*pi*FF/c;
       cos_elevation=real(sqrt(k.^2-(pi/D).^2)./k);
       f_trans=150;
       cos_elevation(FF>=f_trans)=cos_elevation(FF>=f_trans).*0.8;
       
       phasse=-7*ones(size(FF));
       phasse(FF>=350)=4;
       %Gains(:,2:3)=Gains(:,2:3).*cos_elevation.*exp(1i*2*pi*(phasse)/180); %10-Apr-2020 00:02:03.000 287 (110 and -6 null) deg Arctic5G_2014 5G 100-150 Hz 17 m water depth
        Gains(:,2:3)=Gains(:,2:3).*cosd(25).*exp(1i*2*pi*(phasse)/180); %10-Apr-2020 00:02:03.000 287 (110 and -6 null) deg Arctic5G_2014 5G 100-150 Hz 17 m water depth
      
        
    case 'drifterM35'
        Gains(:,1) = getSensitivity(FF,'GTI-M35-300-omni')';
        Gains(:,2) = getSensitivity(FF,'GTI-M35-300-directional')';
        
        %%%Test of whether a true conversion to velcoity is needed:
        %%%dp/dx=rho*p/(i*k)=-i*rho*c*p/(2*pi*f)
        %Gains(:,2)=Gains(:,2)./(2*pi*FF);
        
        Gains(:,3)=Gains(:,2);
        
        % figure;
        %semilogx(FF,20*log10(Gains(:,1)),FF, 20*log10(Gains(:,2)));grid on;
    otherwise
        Gains=ones(length(FF),3);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%function [Ix,Iy]=correct_phase(Ix,Iy,FF,instrument_type,phase_calibration_chc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ix,Iy]=correct_phase(Ix,Iy,FF,instrument_type)
phase_calibration_chc='Arctic5G_2014'; %AlisonDASARCalibration Arctic5G_2014
%phase_calibration_chc='none';

Np=size(Ix,2);

switch instrument_type
    case 'DASAR'
        
        switch(phase_calibration_chc)
            case 'none'
                return
            case 'AlisonDASARCalibration'
                fe = [0 93 148 200 250 360 500]';
                pex = [0 -3 -7.8 -18.3 -25 -60 -90]';
                phaseex = interp1(fe,pex,FF);
                phaseex(isnan(phaseex)) = 0;
                Phasee=exp(1i*(pi/180)*phaseex*ones(1,Np));
                
                
                Ix=Ix.*Phasee;
                Iy=Iy.*Phasee;
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
    case 'drifterM35'
        
        %%%Oddly enough, the directional channels seem 90° out of phase...
        %%%  Perhaps the M35 measured pressure gradient (difference between
        %%%  hydrophones) and not velocity?
        Ix=Ix.*exp(1i*pi/2);
        Iy=Iy.*exp(1i*pi/2);
       
        %%%figure;
        % semilogx(FF,HH_omni,FF, HH_vel);grid on;
        return
    otherwise
        
        return
end

end


