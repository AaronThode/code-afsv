%%%%display_directional_diagram.m
%function [TT,FF,azi,vx,vy,VR_ratio,handles]=display_directional_diagram(handles,x,Fs,Nfft,Nfft_window, ovlap, hdr)
%
% handles.filetype
% handles.myfile
% Optional
%    handles.azigram.sec_avg=(Batch_vars.sec_avg);
%     handles.azigram.climm=(Batch_vars.climm);
%     handles.azigram.brefa=(Batch_vars.brefa);
%     handles.azigram.alg=Batch_vars.alg;
% handles.display_view,'Directionality'
% handles.checkbox_grayscale,'Value'
%  handles.edit_fmax,edit_fmin
function [TT,FF,azi,vx,vy,KV_ratio,handles]=display_directional_diagram(handles,x,Fs,Nfft,Nfft_window, ovlap, hdr)

azi=[];vx=[];vy=[];KV_ratio=[];
if strcmpi(handles.filetype,'PSD')
    return
end

%%%%Two files that can get active intensity directionality: *.gsi and
%%%%*DIFAAR*.wav files.
dn=round((1-ovlap)*Nfft);

if strcmpi(handles.filetype,'gsi')
    
    if length(x(1,:))<Nfft/2
        return
    end
    
    %axes(handles.axes1);
    %figure
    
    %%%Original code that does not use the spectrogram function (switched
    %%%to hamming instead of hanning on July 29 2018)
%     tic
%     M=floor(1+(size(x,2)-Nfft)/dn);
%     FF=linspace(0,Fs,Nfft);
%     FF=FF(1:(Nfft/2));
%     TT=(0:(M-1))*dn/Fs;
%     
%     B=zeros(3,Nfft/2,M);
%     x=x';
%     for J=1:M
%         index=(J-1)*dn+(1:Nfft);
%         xsub=x(index,:);  %xsub is [nsamp nchan];
%         Xsub=fft((hamming(Nfft)*ones(1,3)).*xsub,Nfft);
%         B(:,1:(Nfft/2),J)=Xsub(1:(Nfft/2),:).';
%     end
%     vx=squeeze(real((B(1,:,:).*conj(B(2,:,:)))));
%     vy=squeeze(real((B(1,:,:).*conj(B(3,:,:)))));
%     Nf=Nfft/2;
%     toc
    
    %%%%Test alternate, faster version, used since July 29, 2018
    tic
    M=floor(1+(size(x,2)-Nfft)/dn);
    x=x';
    B=zeros(size(x,2),Nfft/2+1,M);
    [~,FF,TT,PdB] = spectrogram(x(:,1),Nfft,round(ovlap*Nfft),Nfft,Fs);
    for J=1:size(x,2)
        [B(J,:,:),FF,TT] = spectrogram(x(:,J),Nfft,round(ovlap*Nfft),Nfft,Fs);
    end
   % vx=squeeze(real((B(1,:,:).*conj(B(2,:,:)))));
    %vy=squeeze(real((B(1,:,:).*conj(B(3,:,:)))));
    rho=1000;c=1500;
    
    vx=squeeze(((B(1,:,:).*conj(B(2,:,:)))))./(rho*c); %%Convert velocity into velocity units
    vy=squeeze(((B(1,:,:).*conj(B(3,:,:)))))./(rho*c);
    
    PdB=10*log10(squeeze(PdB));
    Nf=Nfft/2+1;  %Should be the same as length(FF)
    toc
    get_newparams=false;

    
elseif ~isempty(strfind(handles.myfile,'DIFAR'))
    M=floor(1+(max(size(x))-Nfft)/dn);
    [vx,vy,TT,FF,PdB]=demultiplex_DIFAR(x,Fs,Nfft,ovlap);
    Nf=length(FF);
    get_newparams=false;

    if exist('hdr','var')
        if ~isfield(hdr,'brefa')
        get_newparams=true;
        end
    end
end

%sec_avg=input('Enter time to average over (sec; 0 does no averaging):');
if ~isfield(handles,'azigram')
    get_newparams=true;
    
elseif ~isfield(handles.azigram,'climm')
    get_newparams=true;
    
end

if get_newparams
    Batch_vars=get_Azigram_Callback(handles);
    handles.azigram.sec_avg=(Batch_vars.sec_avg);
    handles.azigram.climm=(Batch_vars.climm);
    if strcmpi(handles.filetype,'gsi')
        handles.azigram.brefa=(Batch_vars.brefa);
        %%%Don't alter hdr.brefa, which contains the correction.
    else
        handles.azigram.brefa=(Batch_vars.brefa);
        hdr.brefa=eval(handles.azigram.brefa);

    end
    handles.azigram.alg=Batch_vars.alg;
end

% Batch_vars.sec_avg	=	'0.1';	Batch_desc{1}	=	'Seconds to average PSD for long-term display, if "0" no averaging' ;
% Batch_vars.climm='[0 360]'; Batch_desc{2}='Bearing Range Color Scale';
% Batch_vars.brefa='11.7';Batch_desc{3}='Bearing bias/correction (default is 11.7 degrees)';
% Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, 'Vector Sensor Processing');
sec_avg=str2num(handles.azigram.sec_avg);
climm=eval(handles.azigram.climm);
alg_mult=eval(handles.azigram.alg);


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
    TT_avg=zeros(1,Nsnap);
    for J=1:Nsnap
        index=floor((J-1)*(Navg*(1-ovlap)))+(1:Navg);
        vx_avg(:,J)=mean(vx(:,index),2);
        vy_avg(:,J)=mean(vy(:,index),2);
        TT_avg(J)=mean(TT(index));
    end
    vx=vx_avg;
    vy=vy_avg;
    TT=TT_avg;
end

%vx=squeeze(imag((B(1,:,:).*conj(B(2,:,:)))));
%vy=squeeze(imag((B(1,:,:).*conj(B(3,:,:)))));

if strcmpi(handles.display_view,'Directionality')
    if ~handles.checkbox_reactive.Value
        mu = atan2d(real(vx),real(vy));
    else
        mu = atan2d(imag(vx),imag(vy));
        
    end
    azi=bnorm(hdr.brefa+mu);
    
    
    hh=imagesc(TT,FF/1000,azi);
    
    
    titstr=' Azimuth';
    try
        if get(handles.checkbox_grayscale,'Value')==1
            colormap(flipud(gray));
        else
            colormap(hsv);
        end
    catch
        colormap(hsv);
    end
    caxis(climm);
    
    %%%Use alpha adjustment to display intensity as well
    if exist('PdB','var')&&(sec_avg==0)
        set(hh,'AlphaData',PdB);
        set(hh,'AlphaDataMapping','scaled')
        alim(str2double(handles.edit_mindB.String) + [0 str2double(handles.edit_dBspread.String)]);
    end
    
elseif strcmpi(handles.display_view,'EnergyRatio')||strcmpi(handles.display_view,'IntensityEnergyRatio')
    %Uncomment to show reactive intensity
    
    %pressure_autospectrum=squeeze(abs(B(1,:,:)).^2);
    %normalized_velocity_autospectrum=squeeze(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2);
   
    
    work.pressure_autospectrum=squeeze(abs(B(1,:,:)).^2);
    work.normalized_velocity_autospectrum=squeeze(abs(B(2,:,:)).^2+abs(B(3,:,:)).^2);
    
    %%These velocities have units of pressure :
    %%(rho*c*velocity)
    work.intensity_xa=squeeze(((B(1,:,:).*conj(B(2,:,:)))));
    work.intensity_ya=squeeze(((B(1,:,:).*conj(B(3,:,:)))));
    if ~handles.checkbox_reactive.Value
        work.intensity_a=sqrt((real(work.intensity_xa).^2)+(real(work.intensity_ya).^2));
    else
        work.intensity_a=sqrt((imag(work.intensity_xa).^2)+(imag(work.intensity_ya).^2));
    end
                       
     energy_density=0.5*abs(work.normalized_velocity_autospectrum+work.pressure_autospectrum);
    %energy_density=abs(work.pressure_autospectrum);
    
     if ~isempty(sec_avg)&&sec_avg>0
        Ia_avg=zeros(Nf,Nsnap);
        Energy_avg=Ia_avg;
        NV_avg=Ia_avg; PA_avg=Ia_avg;
        for J=1:Nsnap
            index=floor((J-1)*Navg*(1-ovlap))+(1:Navg);
            Iax=mean(work.intensity_xa(:,index),2);
            Iay=mean(work.intensity_ya(:,index),2);
            if ~handles.checkbox_reactive.Value
                Ia_avg(:,J)=sqrt(real(Iax).^2+real(Iay).^2);
            else
                Ia_avg(:,J)=sqrt(imag(Iax).^2+imag(Iay).^2); %reactive intensity
            end
            %Energy_avg(:,J)=mean(energy_density(:,index),2);
            
            NV_avg(:,J)=mean(work.normalized_velocity_autospectrum(:,index),2);
            PA_avg(:,J)=mean(work.pressure_autospectrum(:,index),2);
            
        end
        work.intensity_a=Ia_avg;
        energy_density=0.5*abs(NV_avg+PA_avg);
        
        work.normalized_velocity_autospectrum=NV_avg;
        work.pressure_autospectrum=PA_avg;
        
        TT=TT_avg;
    end
    
    %KV_ratio=(normalized_velocity_autospectrum./pressure_autospectrum);
    %kang=acosd(sqrt(KV_ratio));kang(imag(kang)~=0)=0;
    %%%Effective velocity j/Sp
    if strcmpi(handles.display_view,'IntensityEnergyRatio')
        
        c_eff=work.intensity_a./energy_density;
        imagesc(TT,FF/1000,c_eff);caxis([0 1])
    elseif strcmpi(handles.display_view,'EnergyRatio')
        ratioo=work.normalized_velocity_autospectrum./work.pressure_autospectrum;
        imagesc(TT,FF/1000,10*log10(ratioo));caxis([-20 20])
    end
    
    colormap(jet);axis('xy');figure(gcf);colorbar
            
    
    %%%Set colorscale
    %caxx=str2double(handles.edit_mindB.String)+[0 str2double(handles.edit_dBspread.String)];
    %caxis(caxx);
    titstr='Energy ratio (dB)';
    if get(handles.checkbox_grayscale,'Value')==1
        colormap(flipud(gray));
    else
        colormap(jet);
    end
    
end

grid on
axis('xy')

try
    fmax=str2double(get(handles.edit_fmax,'String'));
    fmin=str2double(get(handles.edit_fmin,'String'));
    if fmax==0
        ylim([0 Fs/2000]);
        set(handles.edit_fmax,'String',num2str(Fs/2000));
    else
        ylim([fmin fmax]);
    end
    %ylim([0 1]);axis('xy')
    
    
    colorbar;
    % set(gcf,'pos',[30   322  1229   426])
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (sec)');ylabel('Frequency (kHz)');
    
    set(handles.text_filename,'String',sprintf('%s, %s ... ', ...
        get(handles.text_filename,'String'),titstr));
catch
    disp('catch in display_directional_diagram');
end
