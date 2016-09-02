
%%%%%%%%%master_rewrite_TSV.m%%%%%%%%%%%
%  Load all *crosschecked*mat files and rewrite TSV files.
%  Useful because TSV files always changing.
clear
%finaldir='/Volumes/macmussel1/Arctic_2008/Processed/Site_05/Shell08_Final_February/morph';
finaldir='.';

path(path,'../CommonScripts.dir');

mydir=pwd;
date_range = '20080821';
fnames=dir([finaldir '/*' date_range '*crosschecked_Huber.mat']);

%run_options.localization_alg='HuberCalibratedKappa';
run_options.localization_alg='Huber';
run_options.bearing_alg='sel_ratio';
run_options.filter_chc='min_weight';
run_options.plot_locations=0;
param.localize.min_weight=0.001;

for I=1:length(fnames)
    disp(sprintf('Loading %s/%s',finaldir,fnames(I).name));
    auto=load([finaldir '/' fnames(I).name]);
    auto.run_options.localization_alg= run_options.localization_alg;
    auto.run_options.bearing_alg=run_options.bearing_alg;
    auto.run_options.filter_chc=run_options.filter_chc;
    auto.param.localize.min_weight=param.localize.min_weight;
    
    if exist('date_range')
        [rawdatadir,outputdir,param,manualdir]=load_pathnames(auto.Icase,auto.param);
        [goodDASAR,goodFile,goodName]=find_DASAR_dates(date_range,auto.Isite,'*',rawdatadir,auto.Icase);
        auto.goodFile=goodFile;
    end

    %[locations,tet]=compute_bearings(auto.locations,auto.station,auto.goodFile,auto.param,run_options.bearing_alg,0);
    [locations,confidence_area{I}]=compute_position(auto.locations,auto.goodFile,auto.param,auto.Icase,auto.Isite,run_options);
    %[fname_tsv,Iwrite]=write_tsv(locations,auto.goodName,run_options.localization_alg);
    for J=1:length(auto.goodName)
        goodDASAR{J}=auto.goodName{J}(1:5);
    end
    %[ind,localized]=read_tsv(fname_tsv,0,Inf,goodDASAR);
    %Compare ind and auto.individual
    
    fname_out=sprintf('%s_%s_%s',auto.goodName{end},'crosschecked_minwght0p0',run_options.localization_alg);
    disp(sprintf('Saving %s',fname_out));
    %save(fname_out,'locations','locations_ctime','station','raw_station','param','goodName','goodFile','Icase','Isite','run_options');
    save(fname_out,'-struct', 'auto');          
end

for I=1:length(fnames),
    figure
    Igood=find(confidence_area{I}>0);
    hist(log10(confidence_area{I}(Igood)/1e6),20)
    title(sprintf('%s, %s algorithm',fnames(I).name,run_options.localization_alg));
end