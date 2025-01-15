%%%%compute_directional_metrics.m%%%%%%%%%%%%%%%%%%%
%function [TT,FF,output_array,PdB,param, Ix,Iy,Iz,polarization]=compute_directional_metrics(x,metric_type, ...
%   Fs,Nfft, ovlap, param,filetype,reactive_flag)
% x:  array of vector data, each row is a separate channel
% metric_type{Iwant}:'Azimuth','Elevation','ItoERatio','KEtoPERatio',
%       'IntensityPhase', 'Polarization','AdditiveBeamforming'
%           metric_type can also be a string.
%           If a cell array output_array will be a cell array
%           If 'wavelet' appended to metric_type string, will use wavelets
%          
% Fs: sampling rate in Hz
% Nfft: FFT size used to compute frequency bin size
% ovlap: fraction ovlap when computing spectrograms
% reactive_flag:  if true, compute reactive intensity.  Vector same size as
%
% param:
%     param.sec_avg=(Batch_vars.sec_avg);
%     param.climm=(Batch_vars.climm);
%     param.brefa=(Batch_vars.brefa);
%     param.alg=Batch_vars.alg;
%     param.phase_calibration='Arctic5G_2014';
%     param.instrument{1} = 'drifterM35sensor' or 'DIFARsensor' or
%     '_VS-209-omnisensor'; needs to start with '_' and end with 'sensor'
%           Each channel needs its own string.
% filetype: dead variable, can be empty
%
% Output:
%  TT,FF, output_array{Nmetric}:  TT vector of times and FF vector of Hz for
%               output_array grid, one for each value in metric_type.
%               If 'polarization' is a metric, output_array{}[freq times metric] will have
%               three dimensions, with third dimension being Stokes metric.
%
%  PdB: power spectral density of pressure autospectrum
%  param:  altered parameters of input param
%  Ix, Iy, Iz:  x and y complex intensity

function [TT,FF,output_array,PdB,param,Ix,Iy,Iz,polarization]=compute_directional_metrics(x,metric_type, ...
    Fs,Nfft, ovlap, param,reactive_flag)

TT=[];FF=[];output_array=[];PdB=[];Ix=[];Iy=[];Iz=[];polarization=[];
if ~iscell(metric_type)
    temp=metric_type;clear metric_type
    metric_type{1}=temp;
end
metric_type=lower(metric_type);

if ~exist('reactive_flag','var')
    reactive_flag=zeros(1,length(metric_type));
end
if size(x,1)<size(x,2)
    x=x.';
end
Nchan=size(x,2);

use_elevation=any(contains(metric_type,'Elevation'));
use_wavelet=any(contains(lower(metric_type),'wavelet'));
if ~isfield(param,'instrument')
    param.instrument='DASAR';
end

if use_wavelet
    Iend=strfind(metric_type,'wavelet')-1;
    metric_type=metric_type(1:Iend);
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
    if max(size(x))<Nfft/2
        disp('Signal sample shorter than Nfft');
        return
    end

    %%%%Check if four channels exist if want elevation
    if use_elevation & Nchan<4
        f = errordlg('Elevation angle requires 4-channels.', 'Elevation button error!', 'modal');
        return
    end

    if use_wavelet
        M=size(x,1);
        dn=1;
        fb=cwtfilterbank('SamplingFrequency',Fs,'SignalLength',size(x,1),'FrequencyLimits',[10 0.95*Fs/2],'VoicesPerOctave',12);
        for J=1:Nchan
            [B(J,:,:),FF,COI] = wt(fb,x(:,J));  %use cwtfilterbank for future speed
        end
        TT=(1:M)/Fs;
        Nf=length(FF);
    else
        dn=round((1-ovlap)*Nfft);
        M=floor(1+(size(x,1)-Nfft)/dn);
        B=zeros(Nchan,Nfft/2+1,M);

        for J=1:Nchan
            fprintf('Spectrogram channel %i\n',J);
            [B(J,:,:),FF,TT] = spectrogram(x(:,J),Nfft,round(ovlap*Nfft),Nfft,Fs);
        end
        Nf=Nfft/2+1;  %Should be the same as length(FF)
        clear x
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Correct gain as needed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gains=zeros(length(FF),Nchan);
    for IJ=1:Nchan
        [sensor_name]=parse_instrument_string(param.instrument{IJ});
        if isempty(sensor_name)
            disp('No sensor identified with this channel')
            return
        end
        Gains(:,IJ)=getSensitivity(FF,sensor_name,'units','uPa/V').';
    end

    % Gains=correct_gain(FF,param.instrument,Nchan);
    %temp=10*log10(abs(squeeze(B(2,:,:))));
    %myfig=gcf;
    %figure; for III=1:3,subplot(3,1,III);imagesc(TT,FF,temp);ylim([0 2000]);colorbar;end
    %figure(myfig);

    for J=1:Nchan
        B(J,:,:)=squeeze(B(J,:,:)).*Gains(:,J);
    end

    %%% If AdditiveBeamforming, conduct here and leave...
    if any(strcmpi(metric_type,'AdditiveBeamforming'))
        % [B(2,:,:),B(3,:,:)]=correct_phase(squeeze(conj(B(2,:,:))),squeeze(conj(B(3,:,:))),FF,param.instrument);
        %B(2:3,:,:)=conj(B(2:3,:,:));

        %%%%Additive beamforming
        thta=param.thta-param.brefa;  %%Convert to local reference
        Iout=find(strcmpi(metric_type,'AdditiveBeamforming'));
        output_array{Iout}=squeeze(B(1,:,:)+sind(thta)*B(2,:,:)+cosd(thta)*B(3,:,:));


        %%%Debug options, including |v|/p
        %output_array{Iout}=squeeze((sqrt(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2))./abs(B(1,:,:)));
        %x0=x(:,1)+sind(thta)*x(:,2)+cosd(thta)*x(:,3);
        return
    end

    
    %%%%Compute all cross products: intensity, velocity cross-products, and
    %%%% energy densities....

     %time and memory, but need this if trying to use transparency
   
    

    pressure_autospectrum=squeeze(abs(B(1,:,:)).^2);

    Ix=squeeze(((B(1,:,:).*conj(B(2,:,:)))));
    Iy=squeeze(((B(1,:,:).*conj(B(3,:,:)))));
    Vxx=abs(squeeze(B(2,:,:)).^2);
    Vyy=abs(squeeze(B(3,:,:)).^2);
    if Nchan>3  %%Needed to compute transport velocity
        Iz=squeeze(((B(1,:,:).*conj(B(4,:,:)))));
        Vzz=abs(squeeze(B(4,:,:)).^2);
    end

    %%Polarization is (vx.*conj(vy))./|Vx|^2+|vy|^2); 
    %%parameter
    polarization_flag=any(contains(lower(metric_type),lower('polarization')));
    if polarization_flag
        Vxy=squeeze(B(2,:,:).*conj(B(3,:,:)));
        if Nchan>3
            Vzx=squeeze(B(4,:,:).*conj(B(2,:,:)));
            Vyz=squeeze(B(3,:,:).*conj(B(4,:,:)));
        end
    end
    %end
    %toc

     %%%%Debug
    %    figure
    %    subplot(3,1,1);imagesc(TT,FF/1000,log10(abs(real(Ix)))); colormap(jet);ylim([0 12]);xlim([0 1]);caxis([10 20]);axis xy
    %    subplot(3,1,2);imagesc(TT,FF/1000,log10(abs(real(Iy)))); colormap(jet);ylim([0 12]);xlim([0 1]);caxis([10 20]);axis xy
    %    subplot(3,1,3);imagesc(TT,FF/1000,log10(abs(real(Iz)))); colormap(jet);ylim([0 12]);xlim([0 1]);caxis([10 20]);axis xy

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
    param.mask=(Batch_vars.mask);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Average cross products, if needed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    avg.PP=zeros(Nf,Nsnap);
    avg.Ix=zeros(Nf,Nsnap);
    avg.Iy=zeros(Nf,Nsnap);
    avg.Vxx=zeros(Nf,Nsnap);
    avg.Vyy=zeros(Nf,Nsnap);

    if Nchan>3
        avg.Vzz=zeros(Nf,Nsnap);
        avg.Iz=zeros(Nf,Nsnap);
    end
    avg.TT=zeros(1,Nsnap);

    for J=1:Nsnap
        index=floor((J-1)*(Navg*(1-ovlap)))+(1:Navg);
        avg.PP=mean(pressure_autospectrum(:,index),2);
        avg.Ix(:,J)=mean(Ix(:,index),2);
        avg.Iy(:,J)=mean(Iy(:,index),2);
        avg.Vxx(:,J)=mean(Vxx(:,index),2);
        avg.Vyy(:,J)=mean(Vyy(:,index),2);
        
        if Nchan>3
            avg.Iz(:,J)=mean(Iz(:,index),2);
            avg.Vzz(:,J)=mean(Vzz(:,index),2);
        end
        avg.TT(J)=mean(TT(index));
    end
    pressure_autospectrum=avg.PP;
    Ix=avg.Ix;Iy=avg.Iy;Vxx=avg.Vxx;Vyy=avg.Vyy;
    if Nchan>3
        Iz=avg.Iz;Vzz=avg.Vzz;
    end
    TT=avg.TT;

    if polarization_flag
        avg.Vxy=zeros(Nf,Nsnap);
        if Nchan>3
            avg.Vzx=zeros(Nf,Nsnap);
            avg.Vyz=zeros(Nf,Nsnap);
        end
        for J=1:Nsnap
            index=floor((J-1)*(Navg*(1-ovlap)))+(1:Navg);
            avg.Vxy(:,J)=mean(Vxy(:,index),2);
            if Nchan>3
                avg.Vzx(:,J)=mean(Vzx(:,index),2);
                avg.Vyz(:,J)=mean(Vyz(:,index),2);
            end
        end
        Vxy=avg.Vxy;
        if Nchan>3
            Vzx=avg.Vzx;Vyz=avg.Vyz;
        end
    end

    clear avg
end  %sec_avg


%%%Compute final parameters
normalized_velocity_autospectrum=Vxx+Vyy;
if Nchan>3
    normalized_velocity_autospectrum=normalized_velocity_autospectrum+Vzz;
end
energy_density=0.5*abs(normalized_velocity_autospectrum+pressure_autospectrum);

if ~use_wavelet
    PdB=4+10*log10(2*pressure_autospectrum./(Nfft*Fs));  %%Power spectral density output
else
    PdB=sqrt(pressure_autospectrum);
end
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Correct for phase misalignment in Ix and Iy and Iz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%AARON fix to allow VS-209
[Ix,Iy,Iz]=correct_phase(Ix,Iy,Iz,FF,param.instrument);


for J=1:length(metric_type)  %%for each request

    switch metric_type{J}
        case lower('Azimuth')
            if ~reactive_flag(J)
                mu = single(atan2d(real(Ix),real(Iy)));  %Compass convention (usually atan2d(y,x))
            else
                mu = single(atan2d(imag(Ix),imag(Iy)));  %Compass convention (usually atan2d(y,x))
            end
            output_array{J}=bnorm((param.brefa)+mu);
            %output_array{J}=bnorm(mu);
        case lower('Elevation')

            if ~reactive_flag(J)
                mu = single(atand(real(Iz)./sqrt((real(Iy)).^2+(real(Ix)).^2)));
            else
                mu = single(atand(imag(Iz)./sqrt((imag(Iy)).^2+(imag(Ix)).^2)));
            end
            %output_array{J}=(param.brefa+mu);
            output_array{J}=(mu);


        case lower('ItoERatio')

            if Nchan>3
                if ~reactive_flag(J)
                    intensity=sqrt((real(Ix)).^2+(real(Iy)).^2+(real(Iz)).^2);
                else
                    intensity=sqrt((imag(Ix)).^2+(imag(Iy)).^2+(imag(Iz)).^2);
                end
            else
                if ~reactive_flag(J)
                    intensity=sqrt((real(Ix)).^2+(real(Iy)).^2);
                else
                    intensity=sqrt((imag(Ix)).^2+(imag(Iy)).^2);
                end
            end

            %%%Effective velocity j/Sp
            output_array{J}=single(intensity./energy_density);
        case lower('KEtoPERatio')
            output_array{J}=single(normalized_velocity_autospectrum./pressure_autospectrum);
        case lower('Polarization')
            %%Third dimension is axis, fourth dimension is whether auto or
            %%cross-intensity
            %%% Attach description as a cell matrix to
            %%% param.polarization_desc

            S0=Vxx+Vyy;
            output_array{J}(:,:,1,1)=single(S0); param.polar_desc{1,1}='S0_xy';
            output_array{J}(:,:,1,2)=single((Vxx-Vyy)./S0);param.polar_desc{1,2}='S1_xy';
            output_array{J}(:,:,1,3)=2*single(Vxy./S0);param.polar_desc{1,3}='S{2,3}_xy';
            output_array{J}(:,:,1,4)=single(sqrt(output_array{J}(:,:,1,2).^2+output_array{J}(:,:,1,3).*conj(output_array{J}(:,:,1,3))));
            %output_array{J}(:,:,1,4)=output_array{J}(:,:,1,4)./output_array{J}(:,:,1,1);
            param.polar_desc{1,4}='DegreePolarization_xy';

            if Nchan>3
                S0=Vyy+Vzz;
                output_array{J}(:,:,2,1)=single(S0);param.polar_desc{2,1}='S0_yz';
                output_array{J}(:,:,2,2)=single((Vyy-Vzz)./S0);param.polar_desc{2,2}='S1_yz';
                output_array{J}(:,:,2,3)=2*single(Vyz./S0);param.polar_desc{2,3}='S{2,3}_yz';


                S0=Vzz+Vxx;
                output_array{J}(:,:,3,1)=single(S0);param.polar_desc{3,1}='S0_zx';
                output_array{J}(:,:,3,2)=single((Vzz-Vxx)./S0);param.polar_desc{3,2}='S1_zx';
                output_array{J}(:,:,3,3)=2*single(Vzx./S0);param.polar_desc{3,3}='S{2,3}_zx';

                for JJJ=2:3
                    output_array{J}(:,:,JJJ,4)=single(sqrt(output_array{J}(:,:,JJJ,2).^2+output_array{J}(:,:,JJJ,3).*conj(output_array{J}(:,:,JJJ,3))));
                    %output_array{J}(:,:,JJJ,4)=output_array{J}(:,:,JJJ,4)./output_array{J}(:,:,JJJ,1);
                end
                param.polar_desc{2,4}='DegreePolarization_yz';
                param.polar_desc{3,4}='DegreePolarization_zx';
            end
%             polarization=(real(B(2,:,:)).*imag(B(3,:,:))-real(B(3,:,:)).*imag(B(2,:,:)));
%             polarization=-2*squeeze(polarization./(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2));
% 
%              polarization=zeros(Nchan-1,Nfft/2+1,M);
%             polarization(1,:,:)=(real(B(3,:,:)).*imag(B(4,:,:))-real(B(4,:,:)).*imag(B(3,:,:)));
%             polarization(1,:,:)=polarization(1,:,:)./(abs(B(3,:,:)).^2+abs(B(4,:,:)).^2);
%             polarization(2,:,:)=(real(B(4,:,:)).*imag(B(2,:,:))-real(B(2,:,:)).*imag(B(4,:,:)));
%             polarization(2,:,:)=polarization(2,:,:)./(abs(B(2,:,:)).^2+abs(B(4,:,:)).^2);
%             polarization(3,:,:)=(real(B(2,:,:)).*imag(B(3,:,:))-real(B(3,:,:)).*imag(B(2,:,:)));
%             polarization(3,:,:)=polarization(3,:,:)./(abs(B(3,:,:)).^2+abs(B(2,:,:)).^2);
%             polarization=0.5*polarization;
% 
% 
%             output_array{J}=squeeze(polarization(3,:,:));
     
        case lower('IntensityPhase')
            if ~reactive_flag(J)
                output_array{J}=atan2d(sqrt(imag(Ix).^2+imag(Iy).^2),sqrt((real(Ix)).^2+(real(Iy)).^2));
            else
                output_array{J}=atan2d(sqrt(real(Ix).^2+real(Iy).^2),sqrt((imag(Ix)).^2+(imag(Iy)).^2));
            end
            output_array{J}=single(output_array{J});
        case lower('IntensityPhaseZ')
            if ~reactive_flag(J)
                output_array{J}=atan2d(imag(Iz),real(Iz));
            else
                output_array{J}=atan2d(real(Iz),imag(Iz));
            end
            output_array{J}=single(output_array{J});
        case lower('PhaseSpeed')
            if Nchan>3
                if ~reactive_flag(J)
                    output_array{J}=1450*pressure_autospectrum./sqrt((real(Ix)).^2+(real(Iy)).^2+(real(Iz)).^2);
                else
                    output_array{J}=1450*pressure_autospectrum./sqrt((imag(Ix)).^2+(imag(Iy)).^2+(imag(Iz)).^2);

                end
            else
                if ~reactive_flag(J)
                    output_array{J}=1450*pressure_autospectrum./sqrt((real(Ix)).^2+(real(Iy)).^2);
                else
                    output_array{J}=1450*pressure_autospectrum./sqrt((imag(Ix)).^2+(imag(Iy)).^2);

                end
            end
            output_array{J}=single(output_array{J});
    end
end

end

%%%%%%%SubFunctions%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%function [Ix,Iy,Iz]=correct_phase(Ix,Iy,Iz,FF,instrument_type,phase_calibration_chc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ix,Iy,Iz]=correct_phase(Ix,Iy,Iz,FF,instrument_type)
phase_calibration_chc='none'; %default choice
%phase_calibration_chc='none';

Np=size(Ix,2);

[sensor_name]=parse_instrument_string(instrument_type{1});
if isempty(sensor_name)
    disp('No sensor identified with this channel')
    return
end

switch sensor_name
    case 'DASAR'
        phase_calibration_chc='Arctic5G_2014'; %AlisonDASARCalibration Arctic5G_2014
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


