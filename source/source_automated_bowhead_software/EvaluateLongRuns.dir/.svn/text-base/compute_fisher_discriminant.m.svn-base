function compute_fisher_discriminant(whale,other,param,Ntrue)
%fisher discriminant
% Canonical discriminat analysis
disp('whale analysis:');
[xfeature{1},S{1},m{1},xraw{1}]=extract_feature_statistics(whale,param.feature.names);
disp('other analysis:');
[xfeature{2},S{2},m{2},xraw{2}]=extract_feature_statistics(other,param.feature.names);

projection_vector=inv(S{1}+S{2})*(m{1}-m{2});
projection_vector=projection_vector/norm(projection_vector);

for I=1:length(param.feature.names)
    disp(sprintf('projection vector onto %s is %6.4f\n',param.feature.names{I},projection_vector(I)));
end


%Histogram single feature, if desired
Ichc=menu('Historgram feature:',param.feature.names);
figure
for I=1:2,
    subplot(2,1,I)
    [N,F]=hist(xraw{I}(Ichc,:),250);
    bar(F,N,'k');grid on
    set(gca,'fontweight','bold','fontsize',14);

    xlabel(param.feature.names{Ichc});
    if I==1,
        ylabel('Whale');
        title('Distribution of durations of whale calls and false alarms');
    else
        ylabel('Other');
    end
end
[xfeature1{1},xfeature1{2},threshold]=plot_fisher_discriminant(xfeature,projection_vector,Ntrue);

keyboard
for I=1:2,
    m{I}=mean(xfeature1{I},2);
    S{I}=zeros(length(param.feature.names),length(param.feature.names));
    for II=1:size(xfeature1{I},2),
        S{I}=S{I}+(xfeature1{I}(:,II)-m{I})*((xfeature1{I}(:,II)-m{I}).');
    end
end

projection_vector=inv(S{1}+S{2})*(m{1}-m{2});
projection_vector=projection_vector/norm(projection_vector);

disp('Second round');
[xfeature1{1},xfeature1{2},threshold]=plot_fisher_discriminant(xfeature1,projection_vector,Ntrue);


for I=1:length(param.feature.names)
    disp(sprintf('Second projection vector onto %s is %6.4f\n',param.feature.names{I},projection_vector(I)));
end


end





%%%%%extract_feature_statistics%%%%%%%
function [x,S,m,xraw]=extract_feature_statistics(structure,names)
%%create vector
Nsamples=0;
Nstations=length(structure);

for JJ=1:Nstations,
    Nsamples=Nsamples+length(structure{JJ}.ctime_min);
end

x=zeros(length(names),Nsamples);
xraw=x;
for II=1:length(names),
    y=[];

    for JJ=1:Nstations,
        y=[y structure{JJ}.(names{II})];
    end

    x(II,:)=scale_feature_range(y,names{II});
    xraw(II,:)=y;
end

m=mean(x,2);
S=zeros(length(names),length(names));
for II=1:Nsamples,
    S=S+(x(:,II)-m)*((x(:,II)-m).');
end

end


%function plot_fisher_discriminant(xfeature,projection_vector,Ntrue)
% Plot projection of multidimensional data onto a single dimension
%    xfeature: two element cell matrix, 1 for whale features, 2 for other
%           Each cell contains [Nfeature x Nsample] matrix.
%    projection_vector: Fisher linear dimscrinant 1-D projection vector
%    Ntrue, number of manual detections (not localizations).
function [xtrue,xfalse,threshold]=plot_fisher_discriminant(x,projection_vector,Ntrue)

figure
Nwhale=size(x{1},2);
Nother=size(x{2},2);

for I=1:2,
    xp{I}=projection_vector'*x{I};
    subplot(3,1,I);
    if I==1,
        [NN{I},XX{I}]=hist(xp{I},200);grid on;
    else
        [NN{I},XX{I}]=hist(xp{I},XX{1});grid on;

    end
    bar(XX{I},NN{I});
    cum{I}=cumsum(NN{I})/sum(NN{I});
    set(gca,'fontweight','bold','fontsize',14);
end
xlimm=xlim;

subplot(3,1,3);
plot(XX{1},cum{1},XX{2},cum{2});
set(gca,'fontweight','bold','fontsize',14);
xlabel('Scalar threshold');
ylabel('Fraction rejected');
grid on
xlim(xlimm);

subplot(3,1,1);
xlim(xlimm);
title('Fisher Discriminant Analysis: Best axis of separation for whale and false alarms')
subplot(3,1,2);
xlim(xlimm);


figure;
plot(cum{1},cum{2},'ko');grid on
set(gca,'fontweight','bold','fontsize',14);
xlabel('missed call fraction');
ylabel('fraction of false alarms cut')
title('ROC curve for discriminant');


figure
subplot(2,1,1);
plot(XX{1},Nwhale*(1-cum{1}),'ro',XX{2},Nother*(1-cum{2}),'ko');
set(gca,'fontweight','bold','fontsize',14);
xlabel('Scalar threshold');
xlim(xlimm);
ylabel('Detectons remaining');grid on;

subplot(2,1,2);

false_alarm_ratio=(1-cum{2})*Nother./((1-cum{1})*Nwhale);
plot(cum{1},false_alarm_ratio,'ko');grid on
set(gca,'fontweight','bold','fontsize',14);

xlabel('missed call fraction');
ylabel('false alarm ratio')

figure(gcf);
disp('Click on bottom panel ... ');
tmp=ginput(1);
miss_ratio=tmp(1);
[junk,Iwant]=min(abs(cum{1}-miss_ratio));
threshold=XX{1}(Iwant);

for I=1:2,
    Ipass{I}=find(xp{I}>threshold);
    Npass(I)=length(Ipass{I});
end

disp(sprintf('%i manual calls: Before filtering, %i true and %i false alarms',Ntrue,size(x{1},2),size(x{2},2)));
disp(sprintf('%i manual calls: After filtering, %i true and %i false alarms',Ntrue,Npass(1),Npass(2)));

xtrue=x{1}(:,Ipass{1});
xfalse=x{2}(:,Ipass{2});
end


function yout=scale_feature_range(yin,feature_str)
switch lower(feature_str)

    case 'robust_fmin'
        yrange=[0 500];
    case 'fmin'
        yrange=[0 500];
    case 'fmax'
        yrange=[0 500];
    case 'robust_fmax'
        yrange=[0 500];
    case 'centroid'
        yrange=[0 500];
    case 'boundingbox'
        yrange=[0 20];
    case 'orientation'
        yrange=[-90 90];
    case'eccentricity'
        yrange=[0 1];
    case  'local_bandwidth'
        yrange=[0 500];
    case'robust_bandwidth'
        yrange=[0 500];
    case'solidity'
        yrange=[0 1];
    case 'intensity_variance',
        yrange=[1 20];
    case 'perimeter_ratio',
        yrange=[0 2];
    case 'snr'
        yrange=[-20 40];
    case 'sel'
        yrange=[ 70 140];
    case 'slope'
        yrange=[-40 40];
    case 'curv'
        yrange=[-40 40];
    case 'centroid_freq'
        yrange=[0 500];
    case 'majoraxislength'
        yrange=[0 200];
    case 'minoraxislength'
        yrange=[0 200];
    otherwise
        keyboard
        yout=yin;
        return
end

yout= (yin-yrange(1))./(yrange(2)-yrange(1));
disp(sprintf('Scaled %s min: %8.4f  max: %8.4f',feature_str,min(yout),max(yout)));
end
