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
B=zeros(3,Nfft/2,M);
x=x';
for J=1:M
    index=(J-1)*dn+(1:Nfft);
    xsub=x(index,:);  %xsub is [nsamp nchan];
    Xsub=fft(xsub,Nfft);
    B(:,1:(Nfft/2),J)=Xsub(1:(Nfft/2),:).';
end

vx=squeeze(real((B(1,:,:).*conj(B(2,:,:)))));
vy=squeeze(real((B(1,:,:).*conj(B(3,:,:)))));

mu = 180/pi*atan2(vx,vy);
FF=linspace(0,Fs,Nfft);
FF=FF(1:(Nfft/2));
TT=0:(M-1)*dn/Fs;
% handles.sgram.T		=	TT;
% handles.sgram.F		=	FF;
% handles.sgram.B		=	B;
% handles.sgram.Nfft	=	Nfft;
% handles.sgram.ovlap	=	ovlap;
% handles.sgram.Fs	=	Fs;

azi=bnorm(hdr.brefa+mu);

imagesc(TT,FF/1000,azi);%

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


 if get(handles.checkbox_grayscale,'Value')==1,
     colormap(flipud(gray));
 else
     colormap(parula);
 end

colorbar;
% set(gcf,'pos',[30   322  1229   426])
set(gca,'fontweight','bold','fontsize',14);
xlabel('Time (sec)');ylabel('Frequency (kHz)');
if ~strcmp(handles.display_view,'Spectrogram')
    title(get(handles.text_filename,'String'));
end