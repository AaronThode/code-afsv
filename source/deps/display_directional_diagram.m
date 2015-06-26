function display_directional_diagram(handles,x,Fs,Nfft,Nfft_window, ovlap, hdr)

if strcmpi(handles.filetype,'PSD')
    return
end

if length(x(1,:))<Nfft/2
    return
end

axes(handles.axes1);

dn=(1-ovlap)*Nfft;
M=floor(1+(size(x,2)-Nfft)/dn);
FF=linspace(0,Fs,Nfft);
FF=FF(1:(Nfft/2));
TT=(0:(M-1))*dn/Fs;

B=zeros(3,Nfft/2,M);
x=x';
for J=1:M
    index=(J-1)*dn+(1:Nfft);
    xsub=x(index,:);  %xsub is [nsamp nchan];
    Xsub=fft((hanning(Nfft)*ones(1,3)).*xsub,Nfft);
    B(:,1:(Nfft/2),J)=Xsub(1:(Nfft/2),:).';
end

vx=squeeze(real((B(1,:,:).*conj(B(2,:,:)))));
vy=squeeze(real((B(1,:,:).*conj(B(3,:,:)))));

%sec_avg=input('Enter time to average over (sec; 0 does no averaging):');
Batch_vars.sec_avg	=	'0';	Batch_desc{1}	=	'Seconds to average PSD for long-term display, if "0" no averaging' ;
Batch_vars	=	input_batchparams(Batch_vars, Batch_desc, 'Vector Sensor Processing');
sec_avg=str2num(Batch_vars.sec_avg);

if ~isempty(sec_avg)&&sec_avg>0
    Navg=floor(sec_avg*Fs/dn);  %Samples per avg
    Nsnap=floor((M-Navg)/((1-ovlap)*Navg));  %Number of averaged samples per window.
    vx_avg=zeros(Nfft/2,Nsnap);
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
    if get(handles.checkbox_grayscale,'Value')==1,
        colormap(flipud(gray));
    else
        colormap(hsv);
    end
elseif strcmpi(handles.display_view,'ReactiveRatio')
    %Uncomment to show reactive intensity
    active=sqrt(vx.^2+vy.^2);
    
    vx=squeeze(imag((B(1,:,:).*conj(B(2,:,:)))));
    vy=squeeze(imag((B(1,:,:).*conj(B(3,:,:)))));
    if ~isempty(sec_avg)&&sec_avg>0
        vx_avg=zeros(Nfft/2,Nsnap);
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
    if get(handles.checkbox_grayscale,'Value')==1,
        colormap(flipud(gray));
    else
        colormap(jet);
    end
end
grid on
axis('xy')
fmax=str2double(get(handles.edit_fmax,'String'));
fmin=str2double(get(handles.edit_fmin,'String'));
if fmax==0,
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
