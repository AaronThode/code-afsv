%%%%display_directional_diagram.m
%function [TT,FF,azi,vx,vy,handles]=display_directional_diagram(handles,x,Fs,Nfft,Nfft_window, ovlap, hdr)
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
function [TT,FF,azi,vx,vy,handles]=display_directional_diagram(handles,x,Fs,Nfft,Nfft_window, ovlap, hdr)

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
    for J=1:size(x,2)
        [B(J,:,:),FF,TT] = spectrogram(x(:,J),Nfft,round(ovlap*Nfft),Nfft,Fs);
    end
    vx=squeeze(real((B(1,:,:).*conj(B(2,:,:)))));
    vy=squeeze(real((B(1,:,:).*conj(B(3,:,:)))));
    Nf=Nfft/2+1;  %Should be the same as length(FF)
   
    toc
    
    
elseif ~isempty(strfind(handles.myfile,'DIFAR'))
    M=floor(1+(max(size(x))-Nfft)/dn);
    [vx,vy,TT,FF]=demultiplex_DIFAR(x,Fs,Nfft,ovlap);
    Nf=length(FF);
end

%sec_avg=input('Enter time to average over (sec; 0 does no averaging):');
get_newparams=false;
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

if ~isempty(sec_avg)&&sec_avg>0
    Navg=floor(sec_avg*Fs/dn);  %Samples per avg
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
    
    mu = 180/pi*atan2(vx,vy);
    azi=bnorm(hdr.brefa+mu);
    
    
    imagesc(TT,FF/1000,azi);
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
    
elseif strcmpi(handles.display_view,'ReactiveRatio')
    %Uncomment to show reactive intensity
    active=sqrt(vx.^2+vy.^2);
    
    %vx=squeeze(imag((B(1,:,:).*conj(B(2,:,:)))));
    %vy=squeeze(imag((B(1,:,:).*conj(B(3,:,:)))));
    if ~isempty(sec_avg)&&sec_avg>0
        vx_avg=zeros(Nf,Nsnap);
        vy_avg=vx_avg;
        for J=1:Nsnap
            index=floor((J-1)*Navg*(1-ovlap))+(1:Navg);
            vx_avg(:,J)=mean(vx(:,index),2);
            vy_avg(:,J)=mean(vy(:,index),2);
        end
        vx=vx_avg;
        vy=vy_avg;
        TT=TT_avg;
    end
    
    reactive=sqrt(vx.^2+vy.^2);
    total=(reactive+active);
    imagesc(TT,FF/1000,10*log10(abs(active./total)));
    climm=[-10 0];
    caxis(climm);
    titstr='active fraction (dB)';
    if get(handles.checkbox_grayscale,'Value')==1
        colormap(flipud(gray));
    else
        colormap(jet);
    end
    % elseif strcmpi(handles.display_view,'Impedance')
    %     %Uncomment to show reactive intensity
    %
    %     vx=squeeze(B(2,:,:));
    %     vy=squeeze(B(3,:,:));
    %     v=sqrt(abs(vx).^2+abs(vy).^2);
    %     p=abs(squeeze(B(1,:,:)));
    %
    %     if ~isempty(sec_avg)&&sec_avg>0
    %         v_avg=zeros(Nfft/2,Nsnap);
    %         p_avg=v_avg;
    %         for J=1:Nsnap
    %             index=floor((J-1)*Navg*(1-ovlap))+(1:Navg);
    %             v_avg(:,J)=mean(v(:,index),2);
    %             p_avg(:,J)=mean(p(:,index),2);
    %         end
    %         v=v_avg;
    %         p=p_avg;
    %         TT=TT_avg;
    %     end
    %
    %
    %     imagesc(TT,FF/1000,20*log10(abs(v./p)));
    %     climm=[-6 10];
    %     caxis(climm);
    %     titstr='Acoustic impedance (dB)';
    %     if get(handles.checkbox_grayscale,'Value')==1,
    %         colormap(flipud(gray));
    %     else
    %         colormap(jet);
    %     end
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
