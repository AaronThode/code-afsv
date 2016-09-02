%function [Imiss,Ifalse,Itrue_man,Itrue_auto,Iredundant_auto]=find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, ...
%min_ctime,max_ctime, match_tol, auto_features,correct_features,)
%% A wrapper for find_similar_elements_ovlap that restricts times between min and max_ctime...
% Input:
%       match_tol: fractional overlap in time (and optionally frequency) required to flag a match.  e.g. 0.5 is 50% ovlap, 0.75 is 75% overlap
%       auto_times: vector of automated times
%       auto_durations: vector of automated durations
%       correct_times: manual times
%       correct_durations: manual durations
%       auto_features:  Additional features, like frequency range for automated
%       correct_features:
%       min,max_ctime: time range to chop
% Iredundant_auto: indicides of auto_times and auto_durations that satisfy
%   criteria, but are a worse fit than the values in Itrue_auto


function [Imiss,Ifalse,Itrue_man,Itrue_auto,Iredundant_auto]=find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, ...
    min_ctime,max_ctime, match_tol,auto_features,correct_features)

if ~exist('match_tol')
   match_tol=0.5; 
end
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



Igooda=find(auto_times>min_ctime&auto_times<max_ctime);
auto_times=auto_times(Igooda);
auto_durations=auto_durations(Igooda);
if exist('auto_features')
    auto_features=auto_features(:,Igooda);
end


if isempty(correct_times)
  %  disp('find_similar_elements_ovlap_timelimit: no manual data present');
    Ifalse=1:length(auto_times);
    return
end

if isempty(auto_times)
   % disp('find_similar_elements_ovlap_timelimit: no auto data present');
    Imiss=1:length(correct_times);
    return
end

%run_options.tol=0.5;
if exist('auto_features')
    [Ifalse,Imiss,Iauto_match,Imatch,~,~,Iredundant]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,match_tol,auto_features,correct_features);
    
else  %compare times only
    [Ifalse,Imiss,Iauto_match,Imatch,~,~,Iredundant]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,match_tol);
    
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
