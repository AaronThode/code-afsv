%master_plot_automated_stage_statistics.m
% function [N,Nplot,Ndiff]=plot_automated_stage_statistics(datadir,searchstr,plotlabel,Nstation)
%   Input:
%       datadir: string with location of automated data
%       search_str:  string with distinctive phrase in final output file
%       plotlabel: title to put on figure
%       Nstation: vector of station indicies...
%
%   Output:
%       N: structure of fields describing number of detections at each stage
%       Nplot: [Nstage by {Nday,Nstation}] matrix of cumulative detections
%       Ndiff: output matrix to be plotted by 'area', result should visually match Nplot

%% datadir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/Processed_2008/Site_02/Shell08_Site2_NeuralNetupdated/morph';
%%plot_automated_stage_statistics(datadir,'Filter',7);
% hh=area(1:7,Ndiff')'
% legend(fliplr(hh),{'raw','interval','image processing','neural net','link','position'})

function [N,Nplot,Ndiff,tdate]=plot_automated_stage_statistics(datadir,searchstr,plotlabel,Istation)

%date_range=expand_date_range(date_str);
fnames=dir([datadir '/*' searchstr '*']);

Ndays=length(fnames);
%Ndays=5;
Nstation=length(Istation);

N.raw=zeros(Ndays,Nstation);
N.interval=zeros(Ndays,Nstation);
N.feature_extraction=zeros(Ndays,Nstation);
N.neural_net=zeros(Ndays,Nstation);
N.links=zeros(Ndays,Nstation);
N.fixes=zeros(Ndays,Nstation);
tdate=zeros(Ndays,1);
for I=1:Ndays
    s_str='abcdefghijkl';

    try
        disp(fnames(I).name);
        results=load([datadir '/' fnames(I).name]);
        
        %Get time point
        Ijunk=find(results.locations_ctime(:)>0);
        tmp=datevec(datenum(1970,1,1,0,0,results.locations_ctime(Ijunk(1))));
        tmp(4:end)=0;
        tdate(I)=datenum(tmp);
        
        fname2=[fnames(I).name(1:22) '_morph'];
        
        %%Find succesfull localizations
        Ipass=[];
        for Il=1:length(results.locations)
            if strcmp(results.locations{Il}.position.outcome,'successful')
                Ipass=[Ipass Il];
            else
                %disp([num2str(Il) ' ' results.locations{Il}.position.outcome]);
                
            end
            
        end
        
        %%Check what stations are used in localization results, and remove
        %%if needed.
        Ibad=-ones(Nstation,1);
        for Id=1:Nstation
            II=Istation(Id);
            for J=1:length(results.goodName)
                if ~isempty(findstr(upper(s_str(II)),results.goodName{J}))
                   Ibad(Id)=1;
                end
                
            end
        end
       
        for Id=1:Nstation  %For every DASAR
            
            if Ibad(Id)<0
                disp(sprintf('Station %i is not present in final localization results, so don''t use',Id));
                continue
            end
            %%Load individual station info
            try
                 
                %%Locate station within results... be aware of missing
                %%station...
                II0=Istation(Id);
                for J=1:length(results.goodName)
                    if ~isempty(findstr(upper(s_str(II0)),results.goodName{J}))
                        II=J;
                    end
                    
                end
                
                fname2(5)=upper(s_str(II0));
                disp(sprintf('loading %s',fname2));
                ind=load([datadir '/' fname2]);
               
                
                N.raw(I,Id)=length(results.raw_station(II).raw_detections.ctime);
                if N.raw(I,Id)~=length(ind.raw_detections.ctime)
                   disp(sprintf(sprintf('!!!! N.raw is %i, ind.raw_detections.ctime is %i',N.raw(I,Id),length(ind.raw_detections.ctime))));
                end
                
                N.interval(I,Id)=length(results.raw_station(II).interval_detections.ctime);
                
                if N.interval(I,Id)~=length(ind.interval_detections.ctime)
                   disp(sprintf(sprintf('!!!! N.interval is %i, ind.interval_dections.ctime is %i',N.interval(I,Id),length(ind.interval_detections.ctime))));
            
                end
                N.original_features=length(ind.best_calls.features);
                
                N.feature_extraction(I,Id)=length(results.raw_station(II).morph_detections.ctime);
                
                %%If image processing has broken up base detection into multiples, adjust numbers of previous stages
                delN=N.feature_extraction(I,Id)-N.original_features;
                if delN>0  %image processing produced more detections than original coarse energy detections
                    N.interval(I,Id)=N.interval(I,Id)+delN;
                    N.raw(I,Id)=N.raw(I,Id)+delN;
                else
                    disp('delN is less than zero');
                end
                
                N.neural_net(I,Id)=length(results.station(II).ctime_min);
                N.links(I,Id)=length(find(results.locations_ctime(:,II)>0));
                N.fixes(I,Id)=length(find(results.locations_ctime(Ipass,II)>0));
            catch
                disp(sprintf('%s not available',fname2));
            end
            
        end
        disp('');
        
        if Ndays>1
            Nplot(1,I)=sum(N.raw(I,:));
            Nplot(2,I)=sum(N.interval(I,:));
            Nplot(3,I)=sum(N.feature_extraction(I,:));
            Nplot(4,I)=sum(N.neural_net(I,:));
            Nplot(5,I)=sum(N.links(I,:));
            Nplot(6,I)=sum(N.fixes(I,:));
        else
            Nplot(1,:)=N.raw(I,:);
            Nplot(2,:)=N.interval(I,:);
            Nplot(3,:)=N.feature_extraction(I,:);
            Nplot(4,:)=N.neural_net(I,:);
            Nplot(5,:)=N.links(I,:);
            Nplot(6,:)=N.fixes(I,:);
            
            
            
        end
    catch
        disp(sprintf('%s not loaded',fnames(I).name));
    end
end

Ndiff=flipud([-diff(Nplot); Nplot(end,:)]);

if Ndays==1
    hh=area(Istation,Ndiff');
    grid on;
    set(gca,'layer','top');   %Put grid lines on top
    legend(fliplr(hh),{'event detection stage','interval removal stage','image processing stage','neural net stage','linking stage','localization stage'})
    title(plotlabel);
    set(gca,'xticklabel',s_str(Istation)')
    xlabel('DASAR station:')
end
end
