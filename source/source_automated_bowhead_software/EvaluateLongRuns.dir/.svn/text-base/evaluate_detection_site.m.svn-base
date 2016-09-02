%%%%evaluate_detection_sequence.m%%%%
% function [false_alarms,missed_calls]=evaluate_detection_site(auto_times,correct_times, datte, station_name, run_options),

function [false_alarms,missed_calls,matched_calls]=evaluate_detection_site(auto_times,correct_times, datte, station_name, run_options),

false_alarms=-1*ones(size(auto_times,1),1);
missed_calls=-1*ones(size(correct_times,1),1);
matched_calls=missed_calls;
matched_auto=false_alarms;

for Icall=1:length(correct_times),
    Nhits=0;
    Igood=1:size(auto_times,1);
    for Id=1:length(station_name),
        if correct_times(Icall,Id)==0,
            continue;
        end

        Icand=find(abs(correct_times(Icall,Id)-auto_times(Igood,Id))<=run_options.tol);

        if ~isempty(Icand),
            Nhits=Nhits+1;
        end
        Igood=Igood(Icand);
        if ~isempty(Igood)&Nhits>=run_options.minStations,
            matched_calls(Icall)=Icall;
            matched_auto(Icall)=min(Igood);
        end

    end


end

missed_calls=find(matched_calls<0);
matched_calls=find(matched_calls>0);

false_alarms=find(matched_auto<0);
matched_auto=find(matched_auto>0);

disp(sprintf('Date: %s  ',datte));
disp(sprintf('Station results: false_alarms: %i false alarm fraction: %6.4f missed calls: %i, missed fraction: %6.4f \n', ...
    length(false_alarms),length(false_alarms)/length(correct_times), ...
    length(missed_calls),length(missed_calls)/length(correct_times)));