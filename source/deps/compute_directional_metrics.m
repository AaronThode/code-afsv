%%%%compute_directional_metrics
%function [TT,FF,output_array,PdB,param]=compute_directional_metrics(x,metric_type, ...
 %   Fs,Nfft, ovlap, param,filetype,reactive_flag)
% x:  array of vector data, each column is a separate channel
% metric type:'Directionality','ItoERatio','KEtoPERatio' ,'IntensityPhase
% filetype: 'DIFAR' or 'gsi'
% reactive_flag:  if true, compute reactive intensity
% Fs: sampling rate in Hz
% Nfft: FFT size used to compute frequency bin size
% ovlap: fraction ovlap when computing spectrograms
% Optional
%     param.sec_avg=(Batch_vars.sec_avg);
%     param.climm=(Batch_vars.climm);
%     param.brefa=(Batch_vars.brefa);
%     param.alg=Batch_vars.alg;

function [TT,FF,output_array,PdB,param]=compute_directional_metrics(x,metric_type, ...
    Fs,Nfft, ovlap, param,filetype,reactive_flag)

vx=[];vy=[];
if strcmpi(filetype,'PSD')
    return
end

%%%%Two files that can get active intensity directionality: *.gsi and
%%%%*DIFAAR*.wav files.
dn=round((1-ovlap)*Nfft);

if strcmpi(filetype,'gsi')
    
    if length(x(1,:))<Nfft/2
        return
    end
    
    %%%%Test alternate, faster version, used since July 29, 2018
    tic
    M=floor(1+(size(x,2)-Nfft)/dn);
    x=x';
    B=zeros(size(x,2),Nfft/2+1,M);
    %[~,FF,TT,PdB] = spectrogram(x(:,1),Nfft,round(ovlap*Nfft),Nfft,Fs);
    for J=1:size(x,2)
        [B(J,:,:),FF,TT] = spectrogram(x(:,J),Nfft,round(ovlap*Nfft),Nfft,Fs);
        
    end
    
    % vx=squeeze(real((B(1,:,:).*conj(B(2,:,:)))));
    %vy=squeeze(real((B(1,:,:).*conj(B(3,:,:)))));
    rho=1000;c=1500;
    
    vx=squeeze(((B(1,:,:).*conj(B(2,:,:)))));
    vy=squeeze(((B(1,:,:).*conj(B(3,:,:)))));
    
    pressure_autospectrum=squeeze(abs(B(1,:,:)).^2);
    normalized_velocity_autospectrum=squeeze(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2);
    energy_density=0.5*abs(normalized_velocity_autospectrum+pressure_autospectrum);
    
    Nf=Nfft/2+1;  %Should be the same as length(FF)
    toc
    get_newparams=false;
    
    
elseif ~isempty(strfind(filetype,'DIFAR'))
    M=floor(1+(max(size(x))-Nfft)/dn);
    [vx,vy,TT,FF,PdB]=demultiplex_DIFAR(x,Fs,Nfft,ovlap);
    Nf=length(FF);
    get_newparams=false;
    
    
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
sec_avg=str2num(param.sec_avg);

%%%%%Average vx and vy, if needed
if ~isempty(sec_avg)&&sec_avg>0
    
    Navg=floor(sec_avg*Fs/dn);  %Samples per avg
    fprintf('%i averages per sample.\n',Navg);
    if Navg==0
        Navg=1;
        sec_avg=dn/Fs;
    end
    Nsnap=floor((M-Navg)/((1-ovlap)*Navg));  %Number of averaged samples per window.
    vx_avg=zeros(Nf,Nsnap);
    vy_avg=vx_avg;
    NV_avg=vx_avg;
    PA_avg=vx_avg;
    TT_avg=zeros(1,Nsnap);
    for J=1:Nsnap
        index=floor((J-1)*(Navg*(1-ovlap)))+(1:Navg);
        vx_avg(:,J)=mean(vx(:,index),2);
        vy_avg(:,J)=mean(vy(:,index),2);
        NV_avg(:,J)=mean(normalized_velocity_autospectrum(:,index),2);
        PA_avg(:,J)=mean(pressure_autospectrum(:,index),2);
        TT_avg(J)=mean(TT(index));
    end
    vx=vx_avg;
    vy=vy_avg;
    
    
    TT=TT_avg;
    energy_density=0.5*abs(NV_avg+PA_avg);
    normalized_velocity_autospectrum=NV_avg;
    pressure_autospectrum=PA_avg;
end  %sec_avg

PdB=10*log10(2*pressure_autospectrum./(Nfft*Fs));  %%Power spectral density output


if ~reactive_flag
    intensity=sqrt((real(vx)).^2+(real(vy)).^2);
else
    intensity=sqrt((imag(vx)).^2+(imag(vy)).^2);
end

switch metric_type
    case 'Directionality'
        if ~reactive_flag
            mu = atan2d(real(vx),real(vy));
        else
            mu = atan2d(imag(vx),imag(vy));
        end
        output_array=bnorm(param.brefa+mu);
    case 'ItoERatio'
        %%%Effective velocity j/Sp
        output_array=intensity./energy_density;
    case 'KEtoPERatio'
        output_array=normalized_velocity_autospectrum./pressure_autospectrum;
    case 'IntensityPhase'
        output_array=atan2d(sqrt(imag(vx).^2+imag(vy).^2),sqrt((real(vx)).^2+(real(vy)).^2));
end

    
end

