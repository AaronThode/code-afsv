%function [Imiss,Ifalse,Itrue_man,Itrue_auto,Iredundant_auto]=evaluate_stations(auto_times,auto_durations,correct_times,correct_durations, ...
%min_ctime,max_ctime, auto_features,correct_features)
%% A wrapper for find_similar_elements_ovlap that restricts times between min and max_ctime...
% Input:
%       auto_times: vector of automated times
%       auto_durations: vector of automated durations
%       correct_times: manual times
%       correct_durations: manual durations
%       auto_features:  Additional features, like frequency range for automated
%       correct_features:
%       min,max_ctime: time range to chop
% Iredundant_auto: indicides of auto_times and auto_durations that satisfy
%   criteria, but are a worse fit than the values in Itrue_auto


function [Imiss,Ifalse,Itrue_man,Itrue_auto,Iredundant_auto]=evaluate_stations(auto_times,auto_durations,correct_times,correct_durations, ...
    min_ctime,max_ctime, auto_features,correct_features)
% auto and correct features have a [Nfeature x Ntime] matrix
Imiss=[];
Ifalse=[];
Itrue_man=[];
Itrue_auto=[];
Iredundant_auto=[];

%auto_times=station.ctime;
Igoodm=find(correct_times>min_ctime&correct_times<max_ctime);
correct_times=correct_times(Igoodm);
correct_durations=correct_durations(Igoodm);
if exist('correct_features')
    correct_features=correct_features(:,Igoodm);
end

if isempty(correct_times)
    disp('evaluate_stations: no manual data present');
    return
end


Igooda=find(auto_times>min_ctime&auto_times<max_ctime);
auto_times=auto_times(Igooda);
auto_durations=auto_durations(Igooda);
if exist('auto_features')
    auto_features=auto_features(:,Igooda);
end
if isempty(auto_times)
    disp('evaluate_stations: no auto data present');
    return
end

run_options.tol=0.5;
if exist('auto_features')
    [Ifalse,Imiss,Iauto_match,Imatch,miss_diff,jun2,Iredundant]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,run_options.tol,auto_features,correct_features);
    
else  %compare times only
    [Ifalse,Imiss,Iauto_match,Imatch,miss_diff,jun2,Iredundant]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,run_options.tol);
    
end
%if ~isempty(Iredundant{1})
%     keyboard;
% end


%If multiple auto candidates have satisfied criteria, include them in the
%   match as well...
%Norg_auto=length(Iauto_match);
for Ired=1:length(Iredundant)
    if ~isempty(Iredundant{Ired})
        %disp('In evalute_stations, Iredundant detected');keyboard;
        if size(Iredundant{Ired},1)>1
            Iredundant{Ired}=Iredundant{Ired}';
        end
        %Iauto_match=sort([Iauto_match Iredundant{Ired}]);
        Iredundant_auto=[Iredundant_auto Iredundant{Ired}];
    end
   
end
Iredundant_auto=sort(unique(Iredundant_auto));

%fprintf('There were %i Iauto_match and %i redundancies were then added \n',Norg_auto,length(Iredundant_auto));
    
%statss.false_fraction=length(Ifalse)/length(correct_times);
%statss.miss_fraction=length(Imiss)/length(correct_times);
%disp(sprintf('Station %s on %s, results: false_alarms: %i false alarm fraction: %6.4f matched calls %i, missed calls: %i, missed fraction: %6.4f \n', ...
%    station_name, datte,length(Ifalse),length(Ifalse)/length(correct_times), ...
%    length(Imatch),length(Imiss),length(Imiss)/length(correct_times)));

%Export corresponding indicies to automated stations...
%Ix_match=Igooda(Iauto_match);
Imiss=Igoodm(Imiss);
Itrue_man=Igoodm(Imatch);
Ifalse=Igooda(Ifalse);
Itrue_auto=Igooda(Iauto_match);

end
