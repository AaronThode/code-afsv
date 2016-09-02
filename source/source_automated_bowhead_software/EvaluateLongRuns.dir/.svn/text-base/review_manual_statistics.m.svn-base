clear ;close all
path(path,'../BulkProcessingScripts')
path(path,'../CommonScripts.dir')
manualdir='/Users/Shared/Projects/Arctic_2007/TSV_files_Shell08';
Isite=1;
Icase='Shell08';
rawdatadir ='/Volumes/macmussel2/Arctic_2008';
detectiondir='/Volumes/macmussel1/Arctic_2008/Processed';
date_range=['20080821';'20080828';'20080906';'20080913';'20080921';'20080929'];
strr='ABCDEFG';
maxtime=Inf;
run_options.min_stations=1;
run_options.success_only=0;
run_options.plot_manual_statistics=0;


for Idd=1:size(date_range,1),
    [goodDASAR{Idd},goodFile{Idd},goodName{Idd},DASAR_str_date{Idd}]=find_DASAR_dates(date_range(Idd,:), ...
        Isite,'*',rawdatadir, Icase);
    DASAR_coords{Idd}=load_DASAR_coords(Icase,Isite,goodFile{Idd});

end

for J=1:size(date_range,1)
    [localized{J},individual{J},Istrip,localized_false{J},individual_false{J}]=load_manual_results_TSV(manualdir, ...
        goodDASAR,date_range(J,:),maxtime,run_options,mean(DASAR_coords{Idd},1));

    %%Analyze false calls.
    Nunits=size(individual_false{J}{1}.wgt,2);
    dvec=zeros(Nunits,1);
    figure
    for K=1:Nunits
        Ihit=find(individual_false{J}{1}.ctime(:,K)>0);
        dvec(K)=length(Ihit);
        subplot(Nunits,1,K)
        vecc=datevec(datenum(1970,1,1,0,0,individual_false{J}{1}.ctime(Ihit,K)));
        NN=histc(vecc(:,4),0:24);
        bar(0:24,NN,'histc');
        grid on;
        
            title(goodDASAR{J}{K});
        
        
    end
    xlabel(sprintf('local hour on Site %i, %s',Isite,date_range(J,:)));
    
    figure
    subplot(2,1,1)
    bar(1:Nunits,dvec);grid on;xlabel('DASAR number');
    
    title(date_range(J,:));
    subplot(2,1,2)
    hist(localized_false{J}{1}.wctype,1:12);grid on;xlabel('Call type')

    date_chc{J}=date_range(J,:);
end


J=menu('date?',date_chc);
Iunit=menu('DASAR?',goodDASAR{J});
Itime=input('Hour?');
 [localized{J},individual{J},Istrip,localized_false{J},individual_false{J}]=load_manual_results_TSV(manualdir, ...
        goodDASAR,date_range(J,:),maxtime,run_options,mean(DASAR_coords{Idd},1));
tabs=datenum(1970,1,1,0,0,individual_false{J}{1}.ctime(:,Iunit));
twant=datenum(str2num(date_range(J,1:4)),str2num(date_range(J,5:6)),str2num(date_range(J,7:8)),Itime,30,0);
Iwant=find(abs(tabs-twant)<=datenum(0,0,0,0,30,0));
disp(sprintf('Times that lie within 1 hour of %i  of %s on DASAR %s...',Itime,date_range(J,:),goodDASAR{J}{Iunit}));
for I=1:length(Iwant)
disp(sprintf('%s %i',datestr(tabs(Iwant(I))'), localized_false{J}{1}.wctype(Iwant(I))));
end