function  [DASAR_coordsn,xg,yg,VMn]=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,linel,Istation)
%DASAR_coords are [NDASAR 2] matrix
% bearings is a vector in degrees, map definition (0=north, increasing
%   clockwise to east).
% Igood, indicies of DASAR_coords to plot
% VM,A,B,ANG: all provided direcly from vmests, in units of meters....
%  linel: length of bearing line in km...
%  Istation: line to emphasize


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
VM=VM/1000;
DASAR_coords=DASAR_coords/1000;
A=A/1000;
B=B/1000;

if ~exist('linel', 'var')
    linel=35; %length of bearing lines in km
end
%subplot(3,1,LL);
xg=mean(DASAR_coords(:,1));
yg=mean(DASAR_coords(:,2));
DASAR_coordsn(:,1)=DASAR_coords(:,1)-xg;
DASAR_coordsn(:,2)=DASAR_coords(:,2)-yg;

plot(DASAR_coords(:,1)-xg,DASAR_coords(:,2)-yg,'r^','markersize',5,'markerfacecolor',[1 0 0]);hold on
set(gca,'fontweight','bold','fontsize',14);
xlabel('Easting (km)');
ylabel('Northing (km)');
axis('equal');
grid on;

%Convert bearings from nautical to mathematical frame.
if ~isempty(bearings)
    bearings=(90-bearings(Igood))*pi/180;
    for I=1:length(Igood)
        XX=DASAR_coords(Igood(I),1)+[0 linel*cos(bearings(I))]-xg;
        YY=DASAR_coords(Igood(I),2)+[0 linel*sin(bearings(I))]-yg;
        if exist('Istation')&&Igood(I)==Istation
            line(XX,YY,'linewidth',5);
        else
            line(XX,YY);
        end
    end
end
if ~isempty(VM)
    VMn=[VM(:,1)-xg VM(:,2)-yg];
    plot(VMn(1),VMn(2),'ks','markerfacecolor',[0 0 0],'markersize',5);
end


%Plot error elipps
if ~isempty(A)
    %ELLIPSE(ra,rb,ang,x0,y0)
    h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
    set(h,'linewidth',0.5);
end

xlim([-20 20]);
ylim([-20 20]);

hold off;

end%plot_location
