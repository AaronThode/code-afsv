%function display_manual_location_data(local,ind,Icall,param,filename)
%  A tool for drawing spectrograms of all manually localized calls, drawing
%  boxes around where each call should be.
%
% Inputs:
%   local:  output of 'load_manual_results_crosschecked', has fields
%           'ctev','wctype', 'Nused', 'comment' used for display purposes here.
%   ind:     output of 'load_manual_results_crosschecked' describes
%           individual arrival times of particular call at various DASARS.
%           Includes fields 'duration','anchor_station','ctime',
%           'flo','fhi','stndb (SNR)','signdb', 'wctype','DASARstr'
%   Icall: integer indexing a particular call in 'local' and 'ind'
%   param: structure of structures containing relevent parameters,
%       including:
%               param.energy.bufferTime (how much time before call start to
%               download)
%               param.Nfft,param.ovlap,param.Fs: spectrogram parameters
%
%    filename:  cell array, each cell containing complete pathname to raw
%                   data file.  e.g. 'goodFile' output of
%                   find_DASAR_dates.m

function ctimes_out=display_manual_location_data(local,ind,Icall,param,filename)

%strr='abcdefghijklmnop';
linewidth=2;
plotbox=1;

strr=ind.DASARstr;
Nstations=size(ind.duration,2);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%display_manual_location_data.m%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

disp(sprintf('Manual localized call %i at %s',Icall,datestr(datenum(1970,1,1,0,0,local.ctev(Icall)))));
disp(sprintf('Call type %6.2f, stations used %i, comment %s', local.wctype(Icall),local.Nused(Icall),local.comment{Icall}));

Istation=ind.anchor_station(Icall);
Inot=find(ind.ctime(Icall,:)==0);
ind.ctime(Icall,Inot)=-Inf;
Igood=find(ind.ctime(Icall,:)>0);
dt=-ind.ctime(Icall,Istation)+ind.ctime(Icall,Igood);

duration=(ind.duration(Icall,Igood));
flo=ind.flo(Icall,Igood);
fhi=ind.fhi(Icall,Igood);
SNR=ind.stndb(Icall,Igood);
level=ind.sigdb(Icall,Igood);
wctype=ind.wctype(Icall,Igood);
ctime=ind.ctime(Icall,Igood);
tabs=datenum(1970,1,1,0,0,ctime-param.energy.bufferTime);

disp(sprintf('Relative arrival times at %s is %s',strr(Igood),mat2str(dt,3)));
for I=1:length(Igood),
    disp(sprintf('Station %s: Flo = %6.2f, Fhi = %6.2f, duration=%6.2f, SNR=%6.2f \n',strr(Igood(I)),flo(I),fhi(I),duration(I),SNR(I)));
end


Nrows=4;
Iplot=0;
Nfigs=1;
figure;
set(gcf,'units','normalized','pos', [0.05 0.05 0.45 0.9]);
for Istation=Nstations:-1:1
    ctimes_out(Istation)=0;
    Iplot=Iplot+1;

    %%Use these parameters to split plots between two separate figure
    %%windows--easier to see.
    %     if Iplot==Nrows+1,
    %         Iplot=1;
    %         Nfigs=Nfigs+1;
    %         figure;
    %         set(gcf,'units','normalized','pos', [0.05+(Nfigs-1)*0.1 0.05 0.45 0.9]);
    %     end
    %     subplot(Nrows,1,Iplot);

    %%For presentations and plots-single figure
    if Iplot==2*Nrows+1,
        Iplot=1;
        Nfigs=Nfigs+1;
        figure;
        set(gcf,'units','normalized','pos', [0.05+(Nfigs-1)*0.1 0.05 0.45 0.9]);
    end
    subplot(Nrows,2,Iplot);


    I=find(Istation==Igood);
    if ~isempty(I),

        tlen=2*param.energy.bufferTime+duration(I);
        ctimes_out(Istation)=ctime(I)-param.energy.bufferTime;
        
        [y,t,head]=readfile(filename{Igood(I)},ctimes_out(Istation),tlen,1,'ctime','calibrate');
        [SS,FF,TT,PP]=spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
        imagesc(TT,FF,10*log10(PP));axis('xy');
        title(sprintf('Station %s, Start time of detection: %s',strr(Istation),ctime2str(ctime(I))));
        set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');
        if plotbox
        hh=line(param.energy.bufferTime + [0 duration(I)],[flo(I) flo(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
        hh=line(param.energy.bufferTime + [0 duration(I)],[fhi(I) fhi(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
        hh=line(param.energy.bufferTime + [0 0],[flo(I) fhi(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
        hh=line(param.energy.bufferTime + duration(I)*[1 1],[flo(I) fhi(I)]);set(hh,'Color',[1 1 1],'linewidth',linewidth);
        end
        %title(sprintf('Manual all %i, station %s, start time %s, type %6.2f, SNR: %6.2f, level %6.2f',Icall,strr(Igood(I)),datestr(tabs(I)),wctype(I),SNR(I),level(I)));
    else
        tlen=10;
        ctimes_out(Istation)=mean(ctime)-5;
        
        [y,t,head]=readfile(filename{Istation},ctimes_out(Istation),tlen,1,'ctime','calibrate');
        spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
        set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');

        title(sprintf('Manual station %s, %s',strr(Istation),ctime2str(mean(ctime-15))));


    end
    
    caxis([40 120])
    colorbar;
    %caxis([-20 60]);colorbar;

    hold on


end


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ');