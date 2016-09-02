%function [station_out,Npass,Ntotal]=filter_with_Nnet(station,net_name,net_dir,my_threshold)
%%Input a station object, and pass features through a sequence of neural
%%networks

function [station_out,Npass,Ntotal]=filter_with_Nnet(station,net_name,net_dir,my_threshold)

%      station: station object, station(Idasar).feature has vectors of interest
%%     net_name{networks}:  name of file
%      net_dir{networks}:  path of file
%      my_threshold[networks]:  if exits, use as threshold

persistent net_name_persist netdata

%error('Must revise station structure to work here');
mydir=pwd;
if exist(net_dir,'dir')==7
    %cd(net_dir);
else
    fprintf('filter_with_Nnet: Error! %s does not exist!\n',net_dir);
    station_out=station;
    Npass=-1;
    Ntotal=-1;
end
%cd ../Neural_Networks.dir

%%%%%%Load neural networks, and let them persist...
load_networks=0;

if isempty(netdata)
    load_networks=1;
else
    for Inets=1:length(net_name)
        if ~strcmp(net_name_persist{Inets},net_name{Inets})
            load_networks=1;
        end
    end
end

if load_networks==1
    cd(net_dir);
    for Inets=1:length(net_name)
        fprintf('Loading %s\n',net_name{Inets});
        netdata{Inets}=load(net_name{Inets});
        net_name_persist{Inets}=net_name{Inets};
        disp(sprintf('Default threshold Net %i is %6.2f',Inets,netdata{Inets}.threshold));
    end
    cd(mydir);
    
    
end
for Inets=1:length(net_name)
    
    
    for Istation=1:length(station)
        if exist('pattern')
            clear pattern
        end
        if isempty(station(Istation).indicies)
            fprintf(' Station %i is empty, skip\n',Istation);
            
            continue
        end
        Nfeatures=length(netdata{Inets}.feature_names);
        %Safety check, are we looking at same parameters?
        
        %%Create input feature vector from station data...
        try
            for I=1:Nfeatures
                if isfield(station(Istation).feature,netdata{Inets}.feature_names{I})
                    pattern(I,:)=station(Istation).feature.(netdata{Inets}.feature_names{I})(1,:);
                elseif strcmp(netdata{Inets}.feature_names{I},'SNR_global')
                    pattern(I,:)=station(Istation).SNR;
                elseif strcmp(netdata{Inets}.feature_names{I},'SEL_global')
                    pattern(I,:)=station(Istation).SEL;
                    
                elseif isfield(station(Istation),netdata{Inets}.feature_names{I})
                    pattern(I,:)=station(Istation).(netdata{Inets}.feature_names{I});
                else
                    disp(sprintf('filter_with_NNet: feature %s missing',netdata{Inets}.feature_names{I}));
                    keyboard;
                    
                end
                
            end
        catch
            fprintf('Feature mismatch in %s\n',net_name{Inets});
            keyboard;
            
        end
        
        %Normalize
        %ptest=mapminmax('apply',pattern,netdata{Inets}.preprocess_info);
        ptest=mapstd('apply',pattern,netdata{Inets}.preprocess_info);
        
        if ~isempty(netdata{Inets}.preprocess_info_pca)
            ptest = processpca('apply',ptest,netdata{Inets}.preprocess_info_pca);
        end
        
        
        target_est = sim(netdata{Inets}.net,ptest);
        
        if ~exist('my_threshold')||isempty(my_threshold)
            %disp(sprintf('\nUsing stored threshold of %6.2f for network %i',netdata{Inets}.threshold,Inets));
            Igood=find(target_est(1,:)>=netdata{Inets}.threshold);
        else
            %disp(sprintf('\nUsing input threshold of %6.2f for network %i',my_threshold(Inets),Inets));
            
            Igood=find(target_est(1,:)>=my_threshold(Inets));
        end
        Npass(Inets,Istation)=length(Igood);
        Ntotal(Inets,Istation)=size(target_est,2);
        fprintf(' Station %i: %i out of %i detections pass\n',Istation,length(Igood),size(target_est,2));
        station(Istation)=trim_station(station(Istation),Igood);
    end  %Istation
end
station_out=station;
end

