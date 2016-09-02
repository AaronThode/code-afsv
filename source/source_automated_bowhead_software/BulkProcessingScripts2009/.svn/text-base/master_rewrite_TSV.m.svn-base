
%%%%%%%%%master_rewrite_TSV.m%%%%%%%%%%%
%  Load all *crosschecked*mat files and rewrite TSV files.
%  Useful because TSV files always changing.
clear
!rm *tsv
%path(path,'../CommonScripts.dir');
%path('../ComputerSpecificScripts.dir',path);
run_options.localization_alg='Huber';
run_options.bearing_alg='sel_ratio';
run_options.filter_chc='min_weight';
run_options.plot_locations=0;


mydir=pwd;
date_range = '2009';
id_str='ThresholdAdjust2_Huber_FilteredLocations';
now_str=datestr(now,30);

fid=fopen(sprintf('%s_%s_qualitycheck.txt',id_str,now_str),'w');
for Isite=1:5
finaldir=sprintf('/Volumes/macmussel1/Arctic_2009/Processed/Site_0%i/Shell09_Site%i_NeuralNetupdated/morph',Isite,Isite);


fnames=dir([finaldir '/*' date_range '*' id_str '.mat']);
%run_options.localization_alg='HuberCalibratedKappa';
for I=1:length(fnames)
    disp(sprintf('Loading %s/%s',finaldir,fnames(I).name));
    try
    auto=load([finaldir '/' fnames(I).name]);
    catch
        disp(sprintf('Failed to load %s',fnames(I).name));
        continue
    end
%     if exist('date_range')
%         [rawdatadir,outputdir,param,manualdir]=load_pathnames(auto.Icase,auto.param);
%         [goodDASAR,goodFile,goodName]=find_DASAR_dates(date_range,auto.Isite,'*',rawdatadir,auto.Icase);
%         auto.goodFile=goodFile;
%     end

    %[locations,tet]=compute_bearings(auto.locations,auto.station,auto.goodFile,auto.param,run_options.bearing_alg,0);
    %[locations,confidence_area{I}]=compute_position(auto.locations,auto.goodFile,auto.param,auto.Icase,auto.Isite,run_options);

    
    [fname_tsv,Nwrite]=write_tsv(auto.locations,auto.goodName,[id_str '_' now_str]);
    disp(sprintf('File: %s, number of detections across stations: %s',fname_tsv,mat2str(Nwrite)));
    
    if I==1
        for J=1:length(auto.goodName)
            goodDASAR{J}=auto.goodName{J}(1:5);
            fprintf(fid,'DASAR file used: %s\n',goodDASAR{J});
        end
    end
    fprintf(fid,'File: %s, number of detections across stations: %s\n',fname_tsv,mat2str(Nwrite));
    
    
    %[ind,localized]=read_tsv(fname_tsv,0,Inf,goodDASAR);
    %Compare ind and auto.individual
end

% for I=1:length(fnames),
%     figure
%     Igood=find(confidence_area{I}>0);
%     hist(log10(confidence_area{I}(Igood)/1e6),20)
%     title(sprintf('%s, %s algorithm',fnames(I).name,run_options.localization_alg));
% end

end
