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
TT=0:(M-1)*dn/Fs;

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
