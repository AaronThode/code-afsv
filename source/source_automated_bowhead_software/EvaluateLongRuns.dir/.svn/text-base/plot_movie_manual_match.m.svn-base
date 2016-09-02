function plot_movie(movie_name,DASAR_coords,Nlongest_link,manual,auto,Ngood_all,Nlocs_linked)

Ireview=find(Nlongest_link>=2)';
Npoints=100;

Ncol=floor(length(Ireview)/Npoints);
Ireview=reshape(Ireview(1:Ncol*Npoints),Npoints,Ncol);


for Iindex=1:size(Ireview,2)
    figure(1);subplot(1,2,1);set(gcf,'pos',[352         296        1197         747]);

    Icall=Ireview(:,Iindex);
    tim1=manual.localized{1}.ctev(Icall(1));
     tim2=manual.localized{1}.ctev(Icall(end));
   

    for I=1:length(Icall)
        Iwant=find(manual.individual{1}.ctime(Icall(I),:)>0);
        VM=[manual.localized{1}.utmx(Icall(I)) manual.localized{1}.utmy(Icall(I))];
        A=manual.localized{1}.axmajor(Icall(I));
        B=manual.localized{1}.axminor(Icall(I));
        Baxis=manual.localized{1}.Ang(Icall(I));
        plot_location(DASAR_coords,[],Iwant,VM,A,B,Baxis,20) ;
        %title(sprintf('N bearings: %i, Manual index %i, time %s',length(Iwant),Icall,ctime2str(manual.localized{1}.ctev(Icall))));
        %disp(sprintf('number bearings: %i, Ngood_all: %i, Nlongest_link %i',length(Iwant),Ngood_all(Icall),Nlongest_link(Icall)));
        xlim([-15 15]);ylim([-15 15]);
        hold on
    end
    title(sprintf('Comparison of %i manual results between %s and %s ',Npoints,ctime2str(tim1),ctime2str(tim2)));
    legend('manual');
    subplot(1,2,2);

    for I=1:length(Icall)
        try
            Icall_auto=Nlocs_linked(Icall(I));
            Iwant2=find(auto.locations_ctime{1}(Icall_auto,:)>0);
            VM2=[auto.locations{1}{Icall_auto}.position.location];
            A2=auto.locations{1}{Icall_auto}.position.major;
            B2=auto.locations{1}{Icall_auto}.position.minor;
            Baxis2=auto.locations{1}{Icall_auto}.position.ellipse_ang;


            plot_location(DASAR_coords,[],Iwant2,VM2,A2,B2,Baxis2,20) ;
           % title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));
            xlim([-15 15]);ylim([-15 15]);
            hold on;
        catch
            %title('Localization failed')
        end
        
    end
    legend('automated')
     for I=1:2
        subplot(1,2,I)
        line([6 6],[-15 15]);
        line([-6 -6],[-15 15]);
    end
    FF(Iindex)=getframe(gcf);
    clf

end
%movie(FF);
movie2avi(FF,movie_name,'fps',1)


%%plot_location(DASAR_coords,bearings,VM,A,B,ANG);
% DASAR_coords UTM, bearings in degrees, VM in UTM
function  plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,linel)

%LL=3;


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
        line(XX,YY);
    end
end
if ~isempty(VM)
    plot(VM(1)-xg,VM(2)-yg,'ro','markerfacecolor',[0 0 0],'markersize',6);
end


%Plot error elipps
% if ~isempty(A)
%     %ELLIPSE(ra,rb,ang,x0,y0)
%     h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
%     set(h,'linewidth',0.5);
% end

hold off;



