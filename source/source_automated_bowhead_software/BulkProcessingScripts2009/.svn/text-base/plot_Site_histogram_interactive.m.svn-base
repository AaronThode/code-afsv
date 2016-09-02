%%%%%plot_site_histogram.m%%%%%%%%%
%%% interactive program to make histogram of calls detected per site,
%%% location, and detection stage.
%  Run from Arctic_2007/Processed folder

clear all
site_names=dir('Site*');
Ichc=menu('Select a site:',site_names(:).name);
site_name=site_names(Ichc).name;

cd(site_name);

loc_names=dir('D*');
Ichc=menu('Select a location:',loc_names(:).name);
loc_name=loc_names(Ichc).name;
cd('FinalResults');

det_options={'raw','noairgun','classified'};
Ichc=menu('Which detection option?',det_options);
det_str=det_options{Ichc};
if strcmp(det_str,'classified')

    eval_str='best_ctimes';
elseif strcmp(det_str,'raw')
    eval_str='best_ctimes_raw';
end

fnames=dir([loc_name '*' det_str '.mat']);
for I=1:length(fnames),
    load(fnames(I).name);
    ctimes=eval(eval_str);
    It=max(findstr('T',fnames(I).name)-4);
    tabs(I)=datenum([fnames(I).name(It+(0:1)) '/' fnames(I).name(It+(2:3))]);
    %temp=datevec(datenum(1970,1,1,0,0,median(ctimes)));
    %temp(:,4:6)=0;
    %tabs(I)=datenum(temp);
    ncount(I)=length(ctimes);
end


%bar(tabs,ncount);
plot(tabs,ncount,'ko');datetick('x',6);
bar(tabs,ncount)