%%create_wctype_database.m%%%
%%% From bulk database of all manually-detected sounds
%%% select a subset of samples for both training and validation.
%%% Select samples in appropriate proportions from each site and date.
clear

%%Parameters of database
Ncalls_per_cat=500;
%SNR_min=10;

fname_sum=dir('*summary.mat');  %For every site summary statistics without clips

%gather statistics on number of call types at each site/date
%  Eliminate detections of the same call at two DASARs...
%best_call_count(Idate,Istation,Iwctype)
for I=1:length(fname_sum)  %For every site
    load(fname_sum(I).name);
    
    Nstations=size(best_calls,2);
    Ndays=size(best_calls,1);
    
    if I==1
        Ncalls=zeros(Ndays,8,3);  %Ncalls and related variables arranged by (day, call type, and site)
        
    end
    
    for J=1:Ndays
        unique_indicies=[];
        call_types=[];
        for K=1:Nstations
            try
                unique_indicies=[unique_indicies best_calls{J,K}.manual_index];
                call_types=[call_types best_calls{J,K}.wctype];
            end
            
        end
        
        [unique_indicies_total{I,J},Iuq]=unique(unique_indicies);  %Ensures that a manually-localized call (row in TSV) is represented once
        call_types_total{I,J}=call_types(Iuq);
        Ncalls(J,:,I)=hist(call_types_total{I,J},1:8);
        
    end
    
    
    
end
%unique_indicies_total{Isite,Iday} are the row numbers in the original TSV file.
%call_types_total{Isite,Iday} is the corresponding call type.
%Ncalls (Idays, call type, Isite) are the number of unique
Ncalls_total=sum(sum(Ncalls,3));  %Total unique calls of each type desired for each dates and sites.
Isample=0;
for Isite=1:3
    fname=dir(sprintf('CallSnips*%i.mat',Isite+2));
    disp(sprintf('loading %s',fname.name));
    if exist('data','var') clear data; end
    data=load(fname.name);
    disp(sprintf('loaded %s',fname.name));
    
    for Icall=1:8
        tmp=Ncalls(:,Icall,Isite);
        Ncalls_want(:,Icall,Isite)=round(2*Ncalls_per_cat.*tmp./Ncalls_total(Icall));  % Total calls desired per type,day, and site
        Ncalls_want(:,Icall,Isite)=2*floor(Ncalls_want(:,Icall,Isite)/2); % Force number of calls to be even...
        
        for Iday=1:Ndays
            disp(sprintf('Site: %i, Call: %i, Day: %i',Isite+2,Icall,Iday));
            if Ncalls_want(Iday,Icall,Isite)==0
                continue
            end
            Icalls_available=find(Icall==call_types_total{Isite,Iday});  %Number of calls of that type available today
            indicies_available=unique_indicies_total{Isite,Iday}(Icalls_available);  %Which manual indicies are this type?
            
            if ~isempty(indicies_available)&&length(indicies_available)>2*Ncalls_want(Iday,Icall,Isite)
                flag=1;
                while flag
                    II=ceil(length(indicies_available)*rand(1,2*Ncalls_want(Iday,Icall,Isite)));  %Factor of three covers replications and two data sets
                    II=unique(II);  %No double counting same sample.
                    disp(sprintf('Length of indicies_available: %i, desired sample: %i',length(indicies_available),2*Ncalls_want(Iday,Icall,Isite)));
                    if length(II)>Ncalls_want(Iday,Icall,Isite)
                        flag=0;
                    end
                end
                II=II(1:(Ncalls_want(Iday,Icall,Isite)));  %Reduce number of selections to those needed for train and validation
                index_want=indicies_available(II);  %Manual indicies to select..
            elseif ~isempty(indicies_available)  %Some replication of samples unavoidable...
                disp(sprintf('Not possible: Length of indicies_available: %i, desired sample: %i',length(indicies_available),2*Ncalls_want(Iday,Icall,Isite)));
                Np=2*floor(length(indicies_available)/2);  %Force indicie numbers to be even..
                index_want=indicies_available(1:Np);
            else
                continue
            end
            
            
            %%Now march through each candidate and locate stations with highest SNR.
            Nstation=size(data.best_calls,2);
            Icol=2;
            dB_level=-1*ones(length(index_want),Nstation);
            temp_wctype=dB_level;
            for Ist=1:Nstation
                if ~isempty(data.best_calls{Iday,Ist})
                    [Idetect,II,JJ]=intersect(index_want,data.best_calls{Iday,Ist}.manual_index);
                    dB_level(II,Ist)=data.best_calls{Iday,Ist}.sigdB(JJ);  %Note
                    temp_wctype(II,Ist)=data.best_calls{Iday,Ist}.wctype(JJ); % should be Icall, no rows should be all zeros
                end
            end
            
            %temp_wctype=zeros(length(index_want),1);
            Isample_start=Isample+1;
            for Isnr=1:length(index_want)
                if Icol==1
                    Icol=2;
                else
                    Icol=1;
                    Isample=Isample+1;
                end
               
               
               [maxdB_level,Ist]=max(dB_level(Isnr,:));
               [Idetect,II,JJ]=intersect(index_want(Isnr),data.best_calls{Iday,Ist}.manual_index);
                
               database.x{Isample,Icol}=data.best_calls{Iday,Ist}.x{JJ};
               database.SNR(Isample,Icol)=data.best_calls{Iday,Ist}.stndb(JJ);
               database.bearing(Isample,Icol)=data.best_calls{Iday,Ist}.bearing(JJ);
               database.dB_level(Isample,Icol)=maxdB_level;
               database.Nstation(Isample,Icol)=sum(dB_level(Isnr,:)>0);
               database.station(Isample,Icol)=Ist;
               database.manual_index(Isample,Icol)=index_want(Isnr);
               database.main_index(Isample,Icol)=JJ;
               database.wctype(Isample,Icol)=data.best_calls{Iday,Ist}.wctype(JJ);
               database.ctime(Isample,Icol)=data.best_calls{Iday,Ist}.ctime(JJ);
               database.tabs(Isample,Icol)=data.best_calls{Iday,Ist}.tabs(JJ);
               database.flo(Isample,Icol)=data.best_calls{Iday,Ist}.flo(JJ);
               database.fhi(Isample,Icol)=data.best_calls{Iday,Ist}.fhi(JJ);
               database.duration(Isample,Icol)=data.best_calls{Iday,Ist}.duration(JJ);
               database.equalization(:,Isample,Icol)=data.best_calls{Iday,Ist}.equalization(:,JJ);
               database.nstart(Isample,Icol)=data.best_calls{Iday,Ist}.nstart(JJ);
               database.short_fname{Isample,Icol}=data.best_calls{Iday,Ist}.short_fname;
               database.site(Isample,Icol)=data.best_calls{Iday,Ist}.site;
              
              
               %Sanity checks
               if (maxdB_level~=data.best_calls{Iday,Ist}.sigdB(JJ))
                   disp('error maxsnr');
                   
               end
               
               if Ist~=data.best_calls{Iday,Ist}.station
                   disp('error Ist')
                   
               end
               %Check: data.best_calls{Iday,Ist}.wctype(index_want(Isnr))=Icall
               %temp_wctype(Isnr)=(data.best_calls{Iday,Ist}.wctype(JJ));
            end  %End of subindex
            if any(database.wctype(Isample_start:Isample)~=Icall)
                
               disp('error wctype')
               keyboard
            end
            
            %Check that number written is even
            if Icol~=2
                disp('Odd number of samples written');
                keyboard
            end
            
        end %Iday
        
    end %Icall
end %Isite
save temp database
for J=1:2
    [tmp,Isort]=sort(database.wctype(:,J));
    names=fieldnames(database);
    for I=1:length(names)
        switch(names{I})
            case 'x'
                database.x(:,J)=database.x(Isort,J);
            case 'equalization'
                database.equalization(:,:,J)=database.equalization(:,Isort,J);
            case 'short_fname'
                database.short_fname(:,J)=database.short_fname(Isort,J);
            otherwise
                database.(names{I})(:,J)=database.(names{I})(Isort,J);
        end
    end
end
database.equalization_freq=data.best_calls{1,1}.equalization_freq;

save TrainingNValidationData2010 database
