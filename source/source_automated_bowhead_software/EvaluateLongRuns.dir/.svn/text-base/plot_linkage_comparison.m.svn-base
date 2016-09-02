function plot_linkage_comparison(movie_name,DASAR_coords,Nlongest_link,manual,auto,Ngood_all,Nlocs_linked,Isite)
%%Warning!  DASAR_coords should only include relevent DASARS...
%Imax=input('Enter maximum number of frames:');

if Isite==1
    xlimm=[-60 60];
    ylimm=[-30 30];
else
    xlimm=[-60 60];
    ylimm=[-30 30];
end
Ireview=find(Nlongest_link>=2)';
Imax=min([30 length(Ireview)]);

Icount=0;
for Icall=Ireview(1:Imax);
    figure(1);subplot(1,2,1);set(gcf,'pos',[352         296        1197         747]);
    Iwant=find(manual.individual{1}.ctime(Icall,:)>0);
    VM=[manual.localized{1}.utmx(Icall) manual.localized{1}.utmy(Icall)];
    A=manual.localized{1}.axmajor(Icall);
    B=manual.localized{1}.axminor(Icall);
    Baxis=manual.localized{1}.Ang(Icall);
    wgt=manual.individual{1}.wgt(Icall,:);

    plot_location(DASAR_coords,manual.individual{1}.bearing(Icall,:),Iwant,wgt,VM,A,B,Baxis,20) ;
    title(sprintf('N bearings: %i, Manual index %i, time %s',length(Iwant),Icall,ctime2str(manual.localized{1}.ctev(Icall))));
    disp(sprintf('number bearings: %i, Ngood_all: %i, Nlongest_link %i',length(Iwant),Ngood_all(Icall),Nlongest_link(Icall)));
    xlim(xlimm);ylim(ylimm);

    try
        Icall_auto=Nlocs_linked(Icall);
        subplot(1,2,2);
        %Iwant2=find(auto.locations_ctime{1}(Icall_auto,:)>0);
        Iwant2=find(~isnan(auto.locations{1}{Icall_auto}.bearing));
        if  isfield(auto.locations{1}{Icall_auto}.position,'location')
            VM2=[auto.locations{1}{Icall_auto}.position.location];
            A2=auto.locations{1}{Icall_auto}.position.major;
            B2=auto.locations{1}{Icall_auto}.position.minor;
            Baxis2=auto.locations{1}{Icall_auto}.position.ellipse_ang;
            wght2=auto.locations{1}{Icall_auto}.position.w;
            plot_location(DASAR_coords,auto.locations{1}{Icall_auto}.bearing,Iwant2,wght2,VM2,A2,B2,Baxis2,20) ;
            title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));

        else
            VM2=[];A2=[];B2=[];Baxis2=[]; wght2=ones(size(Iwant2));
            plot_location(DASAR_coords,auto.locations{1}{Icall_auto}.bearing,Iwant2,wght2,VM2,A2,B2,Baxis2,20) ;

            title(auto.locations{1}{Icall_auto}.position.outcome);
        end
        xlim(xlimm);ylim(ylimm);
    catch
        title('plot_linkage_comparison bug... failed')
        keyboard
    end
    Icount=Icount+1;
    FF(Icount)=getframe(gcf);
    clf


end


%movie(FF);
movie2avi(FF,movie_name,'fps',5);
end

%%plot_location(DASAR_coords,bearings,VM,A,B,ANG);
% DASAR_coords UTM, bearings in degrees, VM in UTM
function  plot_location(DASAR_coords,bearings,Igood,wght,VM,A,B,ANG,linel)

%LL=3;

if isempty(wght)
    wght=ones(size(bearings));
end
if nargin==3,
    VM=[];
    A=[];
    B=[];
    ANG=[];
elseif nargin==4,
    A=[];
    B=[];
    ANG=[];
end
%Convert to km
VM=VM/1000;
DASAR_coords=DASAR_coords/1000;
A=A/1000;
B=B/1000;

if ~exist('linel')
    linel=35; %length of bearing lines in km
end
%subplot(3,1,LL);
xg=mean(DASAR_coords(:,1));
yg=mean(DASAR_coords(:,2));
plot(DASAR_coords(:,1)-xg,DASAR_coords(:,2)-yg,'kp','markerfacecolor',[0 0 0]);hold on
set(gca,'fontweight','bold','fontsize',14);
xlabel('Easting (km)');
ylabel('Northing (km)');
grid on;

%Convert bearings from nautical to mathematical frame.
if ~isempty(bearings)
    bearings=(90-bearings(Igood))*pi/180;
    for I=1:length(Igood)
        XX=DASAR_coords(Igood(I),1)+[0 linel*cos(bearings(I))]-xg;
        YY=DASAR_coords(Igood(I),2)+[0 linel*sin(bearings(I))]-yg;
        hh=line(XX,YY);
        if wght((I))<0.2
            set(hh,'linestyle','--');
        end
    end
    
    
end
if ~isempty(VM)
    plot(VM(1)-xg,VM(2)-yg,'ro','markerfacecolor',[0 0 0]);
end


%Plot error elipps
if ~isempty(A)
    %ELLIPSE(ra,rb,ang,x0,y0)
    h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
    set(h,'linewidth',3);
end

hold off;


end
