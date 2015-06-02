function [x_deconv,source_phase,t_points,f_points] = sourceDeconv_dB(x,Fs,NFFT,Nwindow,flims,Npoint)

N=length(x);
%% plot the tf response to select the points and interpolate the 1st mode
TFR=10*log10(abs(tfrstft(x(:,1),1:N,NFFT,hamming(Nwindow))));
figure('units','normalized','outerposition',[0 0 1 1]);
imagescFun((1:N)/Fs,Fs*(1:NFFT)/NFFT,TFR,'xy');
ylim(flims)
caxis([prctile(TFR(:),70) prctile(TFR(:),100)])

%% Selecting points
disp('Beginning source deconvolution')
title('Choose the points to interpolate the 1st mode')    

% The user draw lines
points_temp(1)=imline();
if Npoint>2
    for ii=2:Npoint-1
        points_temp(ii)=imline();
        pos_av=points_temp(ii-1).getPosition;
        pos=points_temp(ii).getPosition;
        points_temp(ii).setPosition([pos_av(2,:);pos(2,:)]);
    end   
end

% Link 2 consecutives lines when a point is moved
for ii=1:Npoint-2
    points_temp(ii).addNewPositionCallback(@(pos) link1(points_temp,pos,ii+1));
    points_temp(ii+1).addNewPositionCallback(@(pos) link2(points_temp,pos,ii));
end

% Let the user adjust his selection before going on
msgbox('Press any key to continue')
pause()

% Link fonctions
function link1(h,pos,nn)
    pos1=h(nn).getPosition();
    h(nn).setPosition([pos(2,:);pos1(2,:)]);
end
function link2(h,pos,nn)
    pos1=h(nn).getPosition();
    h(nn).setPosition([pos1(1,:);pos(1,:)]);
end

% Retrieving the coordinates
points=zeros(Npoint,2);

for ii=1:Npoint-1
    temp=points_temp(ii).getPosition();
    points(ii,:)=temp(1,:);
end
points(Npoint,:)=temp(2,:);

[t_points,f_points]=deal(points(:,1),points(:,2));
[~,sortIndex]=sort(t_points);
[t_points,f_points]=deal(t_points(sortIndex),f_points(sortIndex));
  

%% Linear interpolation between the selected points
t=linspace(t_points(1),t_points(end),round(1+Fs*(t_points(end)-t_points(1))));
iflaw=interp1(t_points,f_points,t,'linear')';

%% Deconvolution (source phase cancelation)
x_f=fft(x,N);
source=fmodany(iflaw/Fs);
Fsource=fft(source,N);
source_phase=angle(Fsource);
x_deconv=ifft(x_f.*exp(-1i*source_phase),N,'symmetric');
disp('End of source deconvolution')
end