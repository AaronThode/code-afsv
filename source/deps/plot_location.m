function  [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,line0,Istation,VM_extra,units)
%function  [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,linel,Istation,VM_extra,units)

%DASAR_coords are [NDASAR 2] matrix
% bearings is a vector in degrees, map definition (0=north, increasing
%   clockwise to east).
% Igood, indicies of DASAR_coords to plot
% VM,A,B,ANG: all provided direcly from vmests, in units of meters....
%  linel: length of bearing line in km...
%  Istation: line to emphasize
% extra_VM: additional point to plot...
xg=[];yg=[];VMn=[];DASAR_coordsn=[];
if ~exist('units','var')
    units='km';
end
if ~exist('line0', 'var')||isempty(line0)
    line0=3500; %length of bearing lines in m
end
if nargin<10
    VM_extra=[];
end
if nargin==3
    VM=[];
    A=[];
    B=[];
    ANG=[];
elseif nargin==4
    A=[];
    B=[];
    ANG=[];
end
%Convert to km
if strcmpi(units,'km')
    VM=VM/1000;
    DASAR_coords=DASAR_coords/1000;
    A=A/1000;
    B=B/1000;
    line0=line0/1000;
end

%%%%Check that DASAR_coords do not include zero%%%
DASAR_coords=DASAR_coords(DASAR_coords(:,1)~=0,:);
%subplot(3,1,LL);
xg=mean(DASAR_coords(:,1));
yg=mean(DASAR_coords(:,2));
DASAR_coordsn(:,1)=DASAR_coords(:,1)-xg;
DASAR_coordsn(:,2)=DASAR_coords(:,2)-yg;

plot(DASAR_coords(:,1)-xg,DASAR_coords(:,2)-yg,'r^','markersize',5,'markerfacecolor',[1 0 0]);hold on
set(gca,'fontweight','bold','fontsize',14);
axis('equal');
grid on;

%Convert bearings from nautical to mathematical frame.
if ~isempty(bearings)
    bearings=(90-bearings(Igood))*pi/180;
    for I=1:length(Igood)
        XX=DASAR_coords(Igood(I),1)+[0 line0*cos(bearings(I))]-xg;
        YY=DASAR_coords(Igood(I),2)+[0 line0*sin(bearings(I))]-yg;
        if exist('Istation','var')&(Igood(I)==Istation)
            line(XX,YY,'linewidth',5);
        else
            line(XX,YY);
        end
    end
end
if ~isempty(VM)&&all(~isnan(VM))
    VMn=[VM(:,1)-xg VM(:,2)-yg];
    plot(VMn(1),VMn(2),'ks','markerfacecolor',[0 0 0],'markersize',5);
end

if ~isempty(VM_extra)
    VM_extra=VM_extra/1000;
    VMn=[VM_extra(:,1)-xg VM_extra(:,2)-yg];
    plot(VMn(1),VMn(2),'gd','markerfacecolor',[0 0 0],'markersize',5);
end

%Plot error elipps
if ~isempty(A)
    %ELLIPSE(ra,rb,ang,x0,y0)
    h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
    set(h,'linewidth',0.5);
end

if strcmpi(units,'km')
    
    xlim([-30 30]);
    ylim([-30 30]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    
elseif strcmpi(units,'m')
    xlim([-100 100]);
    ylim([-100 100]);
    xlabel('Easting (m)');
    ylabel('Northing (m)');
    
end

hold off;


if isnan(A)||isnan(B)
    title('Ellipse parameters not calculated');
end

if all(isnan(VM))
    title('No convergence');
end

end%plot_location
