function [params,values,slice] = choose_instant(x_deconv,params,delay,Fs,r,c1,c2,D,Nmode)
% input: 
%  x_deconv: time series.  If beta<0 time series has already been
%  time-reversed
%  delay: sample at which to begin warping...
% slice: sum of frequency-warped spectrogram along time axis

tic
values=zeros(3,Nmode); % 3 evaluation function (assess the quality of the warping)

x_min=max(1,delay);               

%Strip signal of parts before (or after, for beta<0) signal5
%x_ok=ifft(fft(x_deconv,N).*exp(1i*2*pi*x_min/N*(1:N)'),'symmetric'); % x_deconv is delayed by x_min
x_ok=[x_deconv(x_min:end); zeros(x_min-1,1)];

%% Warping

%%%JULIEN
if isfield(params,'beta_transform')
    if ~params.beta_transform
        [s_w, Fe_w]=warp_temp_exa(x_ok,Fs,r,c1);    % s_w: warped signal, Fe_w: new warping frequency
    else
        nzero=length(x_deconv)-x_min;
        x_ok=[x_deconv(1:x_min); zeros(nzero,1)];
        %dj=(params.j-x_min)+params.jmax;
        [s_w, Fe_w]=warp_temp_exa_beta(x_ok,Fs,r,c1,params.beta,params.jmax/Fs,x_min/Fs);
        % Note that s_w is no longer time-reversed
    end
    
else
     [s_w, Fe_w]=warp_temp_exa(x_ok,Fs,r,c1);
end

%[s_w, Fe_w]=warp_temp_exa(x_ok,Fs,r,c1);     % s_w: warped signal, Fe_w: new warping frequency
% figure(50)
% subplot(211)
% plot(x_ok)
% subplot(212)
% plot(s_w)
% pause(0.1);


%clear x_ok
M=length(s_w);

if ~isfield(params,'xmax')
    params.xmax=(M-1)/Fe_w;
    params.Nww=1+2*ceil(M/8);
end


M=2^nextpow2(M);
RTF=abs(tfrstft(s_w,1:M,M,hamming(params.Nww)));
M1=floor(M/2);
RTF(M1:end,:)=[];                   % Delete negative frequencies
slice=sum(RTF,2);

t_w=(0:M-1)/Fe_w;                           % Warped time
f_w=(0:(floor(M/2)-2))*Fe_w/M;                         % Warped frequencies
params.t_w=t_w;
params.f_w=f_w;

% Adjust the size of the window and keep the same for each
% time instant
%if ~isfield(params,'fwlims')
    temp=max(RTF,[],2)./max(RTF(:));
    fw_thresh=1e-2;
    Igood=find(temp>=fw_thresh,1,'last');
    fw_max=max(Igood);
    %fw_max=M1-1;
    %fw_thresh=1e-2;
    %while(max(RTF(fw_max,:))/max(RTF(:))<=fw_thresh)
    %    fw_max=fw_max-1;
    %    if fw_max<0
     %       keyboard
    %    end
   % end
    params.fwlims=[0 fw_max]*Fe_w/M1;
    params.clims=[0 max(RTF(:))];
%end

%subplot(2,1,1)
figure(200)
subplot(2,1,1)
%,WINDOW,NOVERLAP,NFFT,Fs
Nfft=128;
spectrogram(x_ok,Nfft,round(0.9*Nfft),Nfft,Fs,'yaxis');
subplot(2,1,2)
imagescFun(t_w,f_w,20*log10(abs(RTF)),'ij')
ylim(params.fwlims)
%xlim([0 params.xmax])
fprintf('M: %i, delay: %i sec\n',M, delay/Fs);
axis('xy');
%caxis(params.clims)



%%% Pekeris cutoff frequencies
if params.beta>=0
    hold on
    pek_cutoff=c1*c2/(2*D*sqrt(c2^2-c1^2));
    colorMode=['r','w','c','y','g','m'];
    for mm=1:Nmode
        %[x_rect,y_rect,zone]=modes_rect(RTF,threshold,zone);
        %plot(x_rect/Fe_w,y_rect*Fe_w/M1,'Color','k','LineWidth',2)
        plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
        %values(:,mm)=evalWarp(RTF,x_rect,y_rect,Fe_w,M1);
    end
    hold off
end
end

