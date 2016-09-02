%function loc_index=-plot_movie_all_auto_only(DASAR_coords,auto,Isite,manual_limits,param,goodFile)
%%% Plot locations of all automated detections, without reference to manual
%%% results
%%   Includes an option select localizations from a map,  return
%%   associated indicies, and create a "false alarm" file
%%
%% Inputs:
%%%     Movie_name: strong of movie

function loc_index=plot_movie_all_auto_only(DASAR_coords,auto,Isite,manual_limits,param,goodFile,Ipass)

strr='ABCDEFGHIJKLMNOP';
review_status='write';  %%Can be 'look' or 'write'
auto_corrected=[];  %This is a vector of autocorrected data..


if isinf(manual_limits)
    disp('Changing infinite limits to 50 km')
    manual_limits=50000;
end

loc_index=[];
if Isite==1
    xlimm=[-30 30];
    ylimm=[-30 30];
else
    %xlimm=1.25*manual_limits*[-1 1]/1000;
    ylimm=[-20 20];
    xlimm=50*[-1 1];
    xlimm=[-10 30];
    
end

%%Find midnight of the day in question, based on time of first detection.. 
It=find(auto.locations_ctime{1}(1,:)>0);
ctime_min=auto.locations_ctime{1}(1,It(1));
tabs_min=(datenum(1970,1,1,0,0,ctime_min));
tvec=datevec(tabs_min);
tabs_midnt=datenum(tvec(1),tvec(2),tvec(3));

sec_elaps=datevec(tabs_min-tabs_midnt);
ctime1=ctime_min-sec_elaps(6)-60*sec_elaps(5)-3600*sec_elaps(4);

if sec_elaps(3)>0
    sec_elaps=datevec(-tabs_min+tabs_midnt);
    ctime1=ctime_min+sec_elaps(6)+60*sec_elaps(5)+3600*sec_elaps(4);
end

ctime2=ctime1+24*60*60;

if ~isfield(param,'Nframes')
    Nframes=25;
    Nframes=7;
else
    Nframes=param.Nframes;
end

if ~isfield(param,'time_inc')
    ctime_range=linspace(ctime1,ctime2,Nframes);
else
    ctime_range=ctime1:(param.time_inc*3600):ctime2;
end


for Iauto=1:size(auto.locations_ctime{1},1);
    if isfield(auto.locations{1}{Iauto}.position,'ctime')
        auto_ctime(Iauto)=auto.locations{1}{Iauto}.position.ctime;
        
    else
        disp('Warning!! Unsuccesful localizations present in data');
    end
end


if exist('Ipass')
    IJ=setdiff(1:size(auto.locations_ctime{1},1),Ipass);
    auto_ctime(IJ)=-1;
end

Icount=0;
for Iindex=1:(length(ctime_range)-1)

    figure(1);set(gcf,'pos',[352         296        1197         747]);
    Icall_auto=find(auto_ctime>=ctime_range(Iindex)&auto_ctime<=ctime_range(Iindex+1));
    Npoint=1;
       
    for I=1:length(Icall_auto)
        try

            %Icall_auto=Nlocs_all_match(Icall(I));
            %Iwant2=find(auto.locations_ctime{1}(Icall_auto(I),:)>0);
            if  isfield(auto.locations{1}{Icall_auto(I)}.position,'location')

                VM2(Npoint,1:2)=[auto.locations{1}{Icall_auto(I)}.position.location];
                Npoint=Npoint+1;
               % A2=auto.locations{1}{Icall_auto(I)}.position.major;
                %B2=auto.locations{1}{Icall_auto(I)}.position.minor;
                %Baxis2=auto.locations{1}{Icall_auto(I)}.position.ellipse_ang;
                %title(sprintf('number bearings: %i, Auto index %i, time %s',length(Iwant2),Icall_auto,ctime2str(auto.locations{1}{Icall_auto}.position.ctime)));
               
            end


        catch
            disp(auto.locations{1}{Icall_auto(I)}.position.outcome);
        end

    end
    try
    Dn=plot_location(DASAR_coords,[],[],VM2) ;
     xlim(xlimm);ylim(ylimm);
                hold on;
    end
    clear VM2

    title(sprintf('%i Automated detections %s - %s',length(Icall_auto),ctime2str(ctime_range(Iindex)),ctime2str(ctime_range(Iindex+1))));

    if Isite==1
        line(-6.8+manual_limits*[1 1],[-1 20]);
        line(-6.8-manual_limits*[1 1],[-1 20]);
        line(10.7368+manual_limits*[1 1],[-15 -5]);
        line(10.7368-manual_limits*[1 1],[-15 -5]);

    else
        line(manual_limits*[1 1]/1000,ylimm);
        line(-manual_limits*[1 1]/1000,ylimm);

    end
    Icount=Icount+1;
    FF(Icount)=getframe(gcf);

    %%Flagging a particular detection
    %Iyes=[];
    Iyes=input('Select a region?');
    
    
    while ~isempty(Iyes)
        disp('Select upper left first, bottom right next');
        tmp=ginput(2);


        %if 1==0
        xg=mean(DASAR_coords(:,1));
        yg=mean(DASAR_coords(:,2));
        tmp=tmp*1000+[1;1]*[xg yg];

        xl=tmp(1,1);
        xu=tmp(2,1);
        yu=tmp(1,2);
        yl=tmp(2,2);
        Iplot=[];
        for I=1:length(Icall_auto)
            try
                VM2=[auto.locations{1}{Icall_auto(I)}.position.location];
            catch
                VM2=[Inf Inf];
            end

            if (VM2(1)>=xl&VM2(1)<=xu&VM2(2)>=yl&&VM2(2)<=yu)
                Iplot=[Iplot Icall_auto(I)];
                h=plot((VM2(1)-xg)/1000,(VM2(2)-yg)/1000,'s','markersize',20);hold on
            end

        end

        Iyes2=input('Display associated spectrograms?');
        if ~isempty(Iyes2)
            Nplots=length(Iplot);
            Nplots=input((sprintf('Enter number of samples(%i), enter -1 to just list sorted times:',Nplots)));
            if Nplots<0
                datestr(datenum(1970,1,1,0,0,auto_ctime(Iplot)))
            else
                for II=Iplot(1:(min([length(Iplot) Nplots])))

                    %%Display spectrogram of links in question
                    param.Nfft=512;param.ovlap=.75; param.Fs=1000;
                    display_automated_crosslink_data(auto.locations{1},II,param,goodFile{1});
                    %title(sprintf('Localization num: %i',II));
                    pause
                    close(2:gcf);
                    filter_names={'Contour_global_bandwidth','Contour_fmax','Contour_fmin'};

                    for K=1:length(filter_names)
                        disp(sprintf('%s: %s',filter_names{K},mat2str(auto.locations{1}{II}.(filter_names{K}))));

                    end
                    
                end
            end
        end
        loc_index=[loc_index Iplot];
        Iyes=input('Select a region?');

    end
    
    clf

end
%movie2avi(FF,movie_name,'fps',5)


%%plot_location(DASAR_coords,bearings,VM,A,B,ANG);
% DASAR_coords UTM, bearings in degrees, VM in UTM
function  DASAR_coordsn=plot_location(DASAR_coords,bearings,Igood,VM,A,B,ANG,linel)

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
DASAR_coordsn(:,1)=DASAR_coords(:,1)-xg;
DASAR_coordsn(:,2)=DASAR_coords(:,2)-yg;

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
    plot(VM(:,1)-xg,VM(:,2)-yg,'ro','markerfacecolor',[0 0 0],'markersize',6);
end


%Plot error elipps
% if ~isempty(A)
%     %ELLIPSE(ra,rb,ang,x0,y0)
%     h=ellipse(A,B,ANG,VM(1)-xg,VM(2)-yg,'k');
%     set(h,'linewidth',0.5);
% end

hold off;



