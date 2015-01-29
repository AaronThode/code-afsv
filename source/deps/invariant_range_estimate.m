%%invariant_range_estimate.m%%
function invariant_range_estimate(R,beta,ca)
%function invariant_range_estimate(R,beta,ca);
% R is in kim

%Given a beamforming plot of elevation angle (deg) vs time (sec),
% pick out a slope and translate into range using waveguide invariant.

if ~exist('beta')
    beta=-9;
end

if ~exist('ca','var')
    ca=1481; %m/s at array depth
end

if ~exist('R','var')
    R=[ 40:10:60 90 100 110];  %km
end
R=R*1000;

bv=R./(ca*beta);
s=-0.15:0.01:0.15;

set(gca,'ytick',s(1:2:end));
disp('Pick a point to draw a curve');
pts=ginput(1);

[junk,Is]=min(abs(pts(2)-s));
hold on

strr='kbrgyckbrgyckbrgyckbrgyc';
for I=1:length(R)
    t=-bv(I).*sqrt(1-s.^2);
    t=t-t(Is)+pts(1)/1000; %msec conversion
    hh(I)=plot(t*1000,s,[strr(I) '--x']);
    text(0.5*median(t)*1000+3.5*I,min(s)-0.05,num2str(R(I)/1000),'color',strr(I),'fontweight','bold','fontsize',18);
end


if 1==0 %Aaron's original
    
    for I=1:length(R)
        set(hh(I),'vis','off')
    end
    disp('Pick two points to make a slope');
    pts=ginput(2);
    t=pts(:,1)/1000;  %sec
    cosang=cos((pts(:,2))*pi/180);
    
    
    dt=diff(t);
    dcos=diff(cosang);
    
    r=-beta*dt*ca/dcos;
    
    titstr=sprintf('Your estimated range is %6.2f km ',r/1000);
    
    hold on
    hh=plot(t*1000,pts(:,2),'k-o');
    title(titstr);
    disp(titstr);
end