function [params,values,slice] = choose_instant(x_deconv,params,delay,Fs,r,c1,c2,D,Nmode)
% slice: sum of frequency-warped spectrogram along time axis

values=zeros(3,Nmode); % 3 evaluation function (assess the quality of the warping)

x_min=max(1,delay);               
N=length(x_deconv);
x_ok=ifft(fft(x_deconv,N).*exp(1i*2*pi*x_min/N*(1:N)'),'symmetric'); % x_deconv is delayed by x_min

%% Warping
[s_w, Fe_w]=warp_temp_exa(x_ok,Fs,r,c1);     % s_w: warped signal, Fe_w: new warping frequency
clear x_ok
M=length(s_w);

if ~isfield(params,'xmax')
    params.xmax=(M-1)/Fe_w;
    params.Nww=1+2*ceil(M/8);
end

t_w=(0:M-1)/Fe_w;                           % Warped time
f_w=(0:floor(M/2))*Fe_w/M;                         % Warped frequencies

params.t_w=t_w;
params.f_w=f_w;
RTF=abs(tfrstft(s_w,1:M,M,hamming(params.Nww)));
M1=floor(M/2);
RTF(M1:end,:)=[];                   % Delete negative frequencies

% Adjust the size of the window and keep the same for each
% time instant
if ~isfield(params,'fwlims')
    fw_max=M1-1;
    fw_thresh=1e-2;
    while(max(RTF(fw_max,:))/max(RTF(:))<=fw_thresh)
        fw_max=fw_max-1;
    end
    params.fwlims=[0 fw_max]*Fe_w/M1;
    params.clims=[0 max(RTF(:))];
end

%subplot(2,1,1)
imagescFun(t_w,f_w,RTF,'ij')
ylim(params.fwlims)
xlim([0 params.xmax])
hold on
pause(0.05);
% subplot(2,1,2);
% theta = 80:1:100;
% [R,xp] = radon(RTF,theta);
% % imshow(R,[],'Xdata',theta,'Ydata',xp,...
% %             'InitialMagnification','fit')
% imagesc(theta,xp,R);
% caxis([0 3e-3])
% xlabel('\theta (degrees)')
% ylabel('x''')
% 
% [~,I90]=min(abs(theta-90));
% slicee=R(:,I90);


%subplot(2,1,1);

%% Pekeris cutoff frequencies
%D=55;
pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
threshold=4*[1 3]; % thresholds to determine the height and the width of the modes (dB)
zone=struct('dpek',pek_cutoff*M1/Fe_w,'pek_cutoff',pek_cutoff*M1/Fe_w);

colorMode=['r','w','c','y','g','m'];
hold on

for mm=1:Nmode
    [x_rect,y_rect,zone]=modes_rect(RTF,threshold,zone);
    plot(x_rect/Fe_w,y_rect*Fe_w/M1,'Color','k','LineWidth',2)
    plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
    values(:,mm)=evalWarp(RTF,x_rect,y_rect,Fe_w,M1);
end

slice=sum(RTF,2);
values(end+1,:)=ones(size(values(end,:)))*kurtosis(slice);


caxis(params.clims)
colorbar
end

