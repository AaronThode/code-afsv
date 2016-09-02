%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% analyze_anchor_call_stats.m  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function [manual, SNR,dates,false,uncertain]=load_manual_results(DASAR_list,dates_list,load_false_reviews,DASAR_site_list),
%%  Aaron Thode
%% function [manual, dates]=load_manual_results(Icase,DASAR_list),
%  Output:
%      manual: cell matrix of length equal to number of days (TSV files in
%      folder)
%      dates: cell matrix of strings of format 'MMDD' corresponding with
%           TSV file
%      SNR:
%
%      false:  Optional output if a 'FalseReview' file has been loaded.
%                   Contains 'ctimes' and 'reason' fields of documented false alarms
%      each cell of manual contains fields defined in read_tsv.m
%      uncertain: Cell array with fields uncertain{Idate}.ctimes.
%               Contains times that will not count for or against automated
%               routine.

%function [manual,SNR,dates,false,uncertain]=load_manual_results(DASAR_list,dates_list,load_false_reviews,DASAR_site_list),
clear all;close all;
DASAR_site_list={'D07s3a','D07s3b','D07s3c','D07s3d','D07s3e','D07s3f','D07s3g', ...
    'D07s4a','D07s4b','D07s4c','D07s4d','D07s4e','D07s4f','D07s4g',...
    'D07s5a','D07s5b','D07s5c','D07s5d','D07s5e','D07s5f','D07s5g'};
dates_list={'0823','0824','0825','0826','0827', ...
    '0828','0829','0830','0831','0901', ....
    '0902', '0903','0904','0905','0906', ...
    '0907','0908','0909', '0910','0911', ...
    '0912','0913','0914','0915','0916', ...
    '0917','0918','0919','0920', '0921', ...
    '0922','0923', ...
    '0924','0925','0926','0927','0929', ...
    '1001','1002','1003','1004','1005', ...
    '1006','1007','1008','1009','1010'
    };
DASAR_list=DASAR_site_list(14);
DASAR_site_list=DASAR_site_list(8:14);
dates_list=dates_list(7:8);
load_false_reviews=1;

duplicate_time_tol=1;  %Remove times within this sec of previous time
minStations=2;  %Minimum Number of DASARS that have detected a call for it to be considered under 'max SNR' criteria

if strcmp(lower(computer),'mac')
    loc.base='/Users/thode/Projects/Greeneridge_bowhead_detection';
    loc.manual=[loc.base '/ManualDetectionDownload/TSV_files_Shell07'];
    path(path,[loc.base '/ManualDetectionDownload']);
elseif strcmp(lower(computer),'maci')
    loc.base='/Users/Shared/Projects/Arctic_2007';
    loc.manual=[loc.base '/TSV_files_Shell07'];
end


mdstr_i=3;  %Start of month/day string

mydir=pwd;
cd(loc.manual);

for J=1:length(dates_list),
    All_SNR=[];

    uncertain{J}.ctimes=[];

    fnames=dir(['*wc' dates_list{J} '*' DASAR_list{1}(end-1) '*.tsv']);  %Site number in DASAR_list
    %for J=1:length(fnames),
    disp(fnames(1).name);
    [manual_raw,localized]=read_tsv(fnames(1).name,0,Inf,DASAR_site_list);

    Ncalls=sum(manual_raw.stndb'>0); %Numbers of stations detecting calls
    Iloc=find(Ncalls>minStations);  %Number of DASARS needed for a location


    for I=1:length(DASAR_site_list),
        Icall=find(manual_raw.wctype(Iloc,I)>0&manual_raw.wctype(Iloc,I)<8);

        Icand=Iloc(Icall);

        %We need to log calls that need not be recognized (low SNR),
        %but shouldn't be counted as false alarms...
        [peakdB,Istation]=max(manual_raw.stndb(Icand,:)');
        Imax=find(Istation==I);

        %Iloc=find(findstr(localized.Outcome(Ilocall(Istation)

        Ifin=Icand(Imax);

        manual{I,J}.wctype=manual_raw.wctype(Ifin,I);
        manual{I,J}.ctime=manual_raw.ctime(Ifin,I);
        manual{I,J}.stndb=manual_raw.stndb(Ifin,I);

        Icall_all=find(manual_raw.wctype(:,I)>0&manual_raw.wctype(:,I)<8);

        All_SNR=[All_SNR ; manual_raw.stndb(Icall_all,I)];


        uncertain_times=manual_raw.ctime(Icall_all,I);
        Iun=[];
        for II=1:length(uncertain_times),
            if all(uncertain_times(II)~=manual{I,J}.ctime),
                Iun=[Iun II];

            end
        end
        [uncertain{J}.ctimes,Iun2]=setdiff(uncertain_times,manual{I,J}.ctime);
        [length(Iun) length(Iun2)]
        %uncertain{J}.ctimes=uncertain{J}.ctimes(Iun);
        disp(sprintf('Total bowhead calls at %s: %i, number that can be localized with local-max SNR: %i, remainder: %i', ...
            DASAR_site_list{I},length(Icall_all),length(Ifin),length(Iun)));

    end

All_SNRs{J}=All_SNR;
end

SNR_vec=2:2:40;
Nvec=length(SNR_vec);

for J=1:length(dates_list),
    for I=1:length(DASAR_site_list),
        SNRdist(:,I)=histc(manual{I,J}.stndb,SNR_vec);


    end
    meanSNR=mean(All_SNRs{J});
    stdSNR=std(All_SNRs{J});
    SNRdist_total=histc(All_SNR,SNR_vec);
    figure
    subplot(2,1,1);
    bar(SNR_vec,SNRdist,'stacked')
    title(sprintf('SNR distribution of anchor calls at %s, Site 4, mean: %6.2f std: %6.2f', dates_list{J}, meanSNR,stdSNR));grid on
    subplot(2,1,2);
    bar(SNR_vec,SNRdist_total,'stacked');
    title('SNR distribution of all calls');grid on
end