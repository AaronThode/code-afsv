%%%master_movie_cleanup.m%%
function master_movie_cleanup
Isite=4;

param.btol=10;  %bearing tolerance in degrees.  If less than one, accept all bearings
param.min_range=10;  %minimum range from center of array in km
param.o_tol=-1; %orientation tolerance in degrees relative to zero
dirpath='../Processed.dir/Site_04/Shell12_Site4_InitialRun/morph/';
%dirpath='.';

fnames=dir(sprintf('%s/S%i*Filt*Filt*.mat',dirpath,Isite));
%fnames=dir(sprintf('%s/S%i*Filt*.mat',dirpath,Isite));

load ../DASARlocations/DASAR_locations_2012.mat
for Ifile=1:length(fnames)
    disp(fnames(Ifile).name);
    yes=menu('Process?','Yes','No');
   
    if yes==2
        continue
    end
    data=load([dirpath '/' fnames(Ifile).name]);
    [Ifinal,param]=cleanup_locations(data.locations,Site,Isite,param);
    
    %%Write TSV file, if desired...
    
    yes=menu('Write a cleaned file?','Yes','No');
    if yes==1
        [PATHSTR,NAME,EXT] = fileparts(fnames(Ifile).name);
        if isempty(PATHSTR),PATHSTR='.';end
        tsv_out=sprintf('%s_stripped',NAME);
        write_tsv(data.locations(Ifinal),data.goodName,tsv_out);
    end
    
end

end


function [Ifinal,param]=cleanup_locations(locations,Site,Isite,param)
close all

if isempty(locations)
    disp('Locations empty');
    return 
end
btol=param.btol;
otol=param.o_tol;

%%Make a crude map..

figure(1);
VM=[mean(Site{Isite}.easting) mean(Site{Isite}.northing)]';
nplot(VM(1),VM(2),'s');hold on
hh=nplot(Site{Isite}.easting,Site{Isite}.northing,'b^');
set(hh,'markersize',10,'markerfacecolor','b');

Nbear=length(locations{1}.bearing);

%Initialize more convenient data structures
feature.bearings=zeros(Nbear,length(locations));
feature.fmin=feature.bearings;
feature.fmax=feature.bearings;
feature.ornt=NaN*ones(Nbear,length(locations));  %Orientation of all components
pos=zeros(2,length(locations));

%%Course is direction of localization from center of site array.  Used for tracking airguns
course=-1*zeros(1,length(locations));
distt=course;  %distnace from array center in km
survey_flag=zeros(length(locations),1);  %Set to one if all bearings are from same direction (distant airgun array)
ship_flag=survey_flag;


%%Extract features
for I=1:length(locations)
    
    %Unpack possible 'giveaway' features
    feature.bearings(:,I)=locations{I}.bearing;
    feature.fmin(:,I)=locations{I}.Totalfmin;
    feature.fmax(:,I)=locations{I}.Totalfmax;
    feature.ornt(:,I)=[locations{I}.feature.Orientation]';
    
    %%Load position and compute course
    try
        pos(:,I)=locations{I}.position.location';
        course(I)=atan2(pos(2,I)-VM(2),pos(1,I)-VM(1));
        course(I)=pi/2-course(I);
        distt(I)=sqrt(sum((pos(:,I)-VM).^2))/1000;
    catch
        course(I)=NaN;
        distt(I)=NaN;
    end
    
    %Crude measure of entropy of background.
    Ival=find((feature.bearings(:,I))~=0);
    p=locations{I}.equalization(Ival,:)';
    p=p./(ones(size(p,1),1)*sum(p));
    feature.entropy(I)=median(sum(p.*log(p)));
    
end

%Convert course to degrees and ensure it is between 0 and 360 degrres
course=bnorm(course*(180/pi));

%Analyze all bearings from all DASARs.
btot=feature.bearings(:);
btot=btot(btot>0);

figure
subplot(2,1,1)
hist(btot,0:5:360);
xlabel('DASAR bearings (deg)');
grid on

subplot(2,1,2)
hist(course,0:5:360);grid on;
xlabel('Azimuth of detection relative to array center');

if btol>0
    
    done=false;
    param.bearing_want=[];
    while ~done
        disp('Select an airgun survey bearing:');
        tmp=ginput(1);
        
        if isempty(tmp)
            done=true;  %No more bearings...
        else
            param.bearing_want=[param.bearing_want tmp(1)];
            fprintf('Bearings selected are %6.2f \n',param.bearing_want');
        end
    end
    
    for I=1:length(survey_flag)
        Ival=find(feature.bearings(:,I)>0);
        
        %Option 1: identify survey by bearing relative to array center.
        %  Problem:  If array not symmetrical then bearing from center may not match
        %       bearing "streak"
        %survey_flag(I)=(abs(course(I)-bearing_want)<btol);
        
        %Option 2:  look for consistentcy in bearings.  This may be more robust for Site 4 (asymmetric sites)
        for J=1:length(param.bearing_want)
            survey_flag(I)=survey_flag(I)||all(abs(feature.bearings(Ival,I)-param.bearing_want(J))<=btol);
        end
        
        %%Potential restrict range to be greater than a threshold distance to array center...
        survey_flag(I)=survey_flag(I).*(distt(I)>param.min_range);
        
    end
    
end %if btol > 0

%Restrict locations to succesfull locations only
Ivalid=find(sum(pos)>0);  %Only use succesful positions
course=course(Ivalid);
pos=pos(:,Ivalid);
feature.bearings=feature.bearings(:,Ivalid);
feature.fmin=feature.fmin(:,Ivalid);
feature.fmax=feature.fmax(:,Ivalid);
feature.ornt=feature.ornt(:,Ivalid);
feature.entropy=feature.entropy(Ivalid);

survey_flag=survey_flag(Ivalid);
ship_flag=ship_flag(Ivalid);
Igun=find(survey_flag==1);
Iother=find(survey_flag==0);


%Plot flagged and other localizations from bearing spike method...
figure(1);
nplot(pos(1,Igun),pos(2,Igun),'ro');
nplot(pos(1,Iother),pos(2,Iother),'go');
grid on
hold on

xlim([-100 100]);
ylim([-100 100]);

%Translate into raw locations
%Iother=Ivalid(Iother);


%Plot properties of remaining detectinos
fmn=feature.fmin(:,Iother);
fmx=feature.fmax(:,Iother);
ornt=feature.ornt(:,Iother);  %Orientation show properties of everything, not just surveys...

fmn=fmn(:);fmx=fmx(:);ornt=ornt(:);

%fmn=fmn(Igood);fmx=fmx(Igood);

figure(2);
subplot(3,1,1);
hist(fmn(fmn>0),100);grid on
title('Minimum frequency of surviving features');

subplot(3,1,2)
hist(fmx(fmx>0),100);grid on
title('Maximum frequency of surviving features');

subplot(3,1,3);
hist(ornt,100);grid on
title('Orientation of surviving features');

if otol>0
    for I=1:length(Iother)
        
        II=Iother(I);
        Ival=find(~isnan(feature.ornt(:,II)));
        ship_flag(II)=all(abs(feature.ornt(Ival,II))<=otol);
        
        %if ship_flag(II)==1
        %   keyboard
        %end
    end
end

%%Plot survivors
figure(1);
Iremain=find(ship_flag+survey_flag==0);
[~,Xn,Yn]=nplot(pos(1,Iremain),pos(2,Iremain),'ko');
legend('Aray Center','DASAR','original','no survey','final');
%disp('Select two points to define a region:');
%tmp=ginput(2);
%Xmin=min(tmp(:,1));Xmax=max(tmp(:,1));
%Ymin=min(tmp(:,2));Ymax=max(tmp(:,2));

%Ifocus=find(Xn>Xmin&Xn<Xmax&Yn>Ymin&Yn<Ymax);
%Irest=setdiff(1:length(Xn),Ifocus);

%
% figure(10);
% subplot(2,1,1)
% hist(feature.entropy(Ifocus),100);title('Entropy selected');
% subplot(2,1,2)
% hist(feature.entropy(Irest),100);title('Entroyp rest');
% xlimm=xlim;
% subplot(2,1,1)
% xlim(xlimm);
%keyboard

Ifinal=Ivalid(Iremain);



%%plot with normalized northings and eastings
    function [hh,Xn,Yn]=nplot(X,Y,sym)
        Xn=(X-VM(1))/1000;
        Yn=(Y-VM(2))/1000;
        hh=plot(Xn,Yn,sym);
    end
end

function bearings=bnorm(bearings)
Ibig=find(bearings>=360);
bearings(Ibig)=bearings(Ibig)-360;

Ilow=find(bearings<0);
bearings(Ilow)=bearings(Ilow)+360;

end
