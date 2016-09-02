%%%review_and_trim_database.m%%%
% Add a "review" field to database to store manual corrections, but not write over original classification
% Codes:-128: not reviewed yet. -1: Eliminate from dataset.  0: reviewed, no change. Integer>0: Reclassify as
%  1:'Upsweep',2: 'Downsweep',3:'Constant',4:'U-shaped',5:'N-Shaped',6:'Undulation',7:'Complex',8:'Bearded Seal'
clear
load TrainingNValidationData2010_500samples.mat
%load ../TSV_files_Shell10/Shell10_ManualTSVFiles/wc0829S310_onlyHalf.mat

Icol=2;  %Which column to review?  1=Training, 2=Validation.
Icall=1;  %Which call type to review?

Nfft=256;
Fs=1000;
ovlap=7/8;
wctypes={'Upsweep','Downsweep','Constant','U-shaped','N-Shaped','Undulation','Complex','Bearded Seal','Accept','Reject'};

if ~isfield(database,'review')
    database.review=-128*ones(size(database.SNR),'int8');
end
Igood=find(database.wctype(:,Icol)==Icall&(database.review(:,Icol)==-128));
SNR=database.SNR(Igood,Icol);

for I=1:length(Igood)
    if rem(I,100)==0, save TrainingNValidationData2010_500samples database;end
    [S,F,T,PP2]=spectrogram(database.x{Igood(I),Icol},Nfft,round(ovlap*Nfft),Nfft,Fs,'yaxis');
    imagesc(T,F,10*log10(abs(PP2)));axis('xy');caxis([40 120])
    
    xlimm=xlim;
    titlestr=sprintf('%s number %i; SNR: %6.2f, DASAR: %i,ctime: %s, filename: %s, Site: %i', ...
        wctypes{Icall},I,database.SNR(Igood(I),Icol),database.station(Igood(I),Icol),ctime2str(database.ctime(Igood(I),Icol)), ...
        database.short_fname{Igood(I),Icol},database.site(Igood(I),Icol));
    title(titlestr);%caxis([-20 40]);
    flo=database.flo(Igood(I),Icol);
    fhi=database.fhi(Igood(I),Icol);
    duration=database.duration(Igood(I),Icol);
    ctime=database.ctime(Igood(I),Icol);
    tabs=database.tabs(Igood(I),Icol);
    tmp=datevec(tabs-datenum(1970,1,1,0,0,ctime));
    
    
    linewidth=2;
    tplot=tmp(end)+[0 duration];
    hold on
    hh=line(tplot,flo*[1 1]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
    hh=line(tplot,fhi*[1 1]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
    hh=line(tplot(1)*[1 1],[flo fhi]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
    hh=line(tplot(2)*[1 1],[flo fhi]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
    text(1,300,sprintf('dB level: %6.2f',database.dB_level(Igood(I),Icol)),'fontweight','bold','fontsize',24,'color',[1 1 1 ]);
    try
        text(1,400,sprintf('Other non-zero stations: %i',database.Nstation(Igood(I),Icol)),'fontweight','bold','fontsize',24,'color',[1 1 1 ]);
    end
    
    hold off
    set(gcf,'pos',[ 206         187        1247         883]);
    %disp(wctypes)
    
    newcode=menu('Choice:',wctypes);
    switch newcode
        case 9
            database.review(Igood(I),Icol)=0;
        case 10
            database.review(Igood(I),Icol)=-1;
        otherwise
            database.review(Igood(I),Icol)=newcode;
    end
    %newcode=input('Enter a different code to reassign, -1 to reject:');
%     if ~isempty(newcode)
%         if newcode==-1
%             database.review(Igood(I),Icol)=-1;
%         else
%             database.review(Igood(I),Icol)=newcode;
%         end
%     else
%         database.review(Igood(I),Icol)=0;  %Mark as reviwed, no change
%     end
end

save TrainingNValidationData2010_500samples database
