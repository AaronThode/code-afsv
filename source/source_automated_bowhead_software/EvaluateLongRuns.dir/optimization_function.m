%function fit=optimization_function(params),
%  Given a parameter vector, run call detection multiple times...
%% DELPHINE, this may seem confusing and large, but the actual
%% optimization commands are very simple...

%% April 30, major overhaul, adding morphological processing and
%% filtering..
%% September 14--add noise-only components as row 15 in TEST DATA, change
%%      scoring so that penalities are applied, and frequencies of calls
%%      used in evaluation.   Time duration of call also considered, but
%%      since call duration not accurate manually, weighted less.

function fit=optimization_function(opt_vec,manual,date_str_new,goodFile,goodName,goodDASAR,run_options,param_in)
%persistent data_all head Interval best_calls best_auto_ctime_org

%%DELPHINE, note that this same 'global' command is in the
%%master_optimizate_script.  This allows us to access large variables
%% without having to pass them through the function.  For technical
%% reasons, this saves RAM memory on the computer.
global  final Iruns

Iruns=Iruns+1;
%fit = sum(opt_vec.^2);
%return
datestr(now)
%try,

%%DELPHINE, it is very important that an optimization function never
%%crash-it must always return a value.  So use 'try' , 'catch' statements.
try
    %if isempty(data_all),
    %    data_all{1}=[];
    %end
    disp(sprintf('Normalized input vector: %s ', mat2str(opt_vec,3)));
    tic

    %%DELPHINE, I find it convenient to define parameters for constructing
    % optimization function in a separate variable 'final'
    %%Define optimization criteria

    tol=final.criteria.tol;
    fit_weight=final.criteria.fit_weight;
    fit_goal=final.criteria.fit_goal;

    %     fit_weight=[1000 1];
    %     fit_goal=[0.05 0.25];

    total.missed=0;
    total.false=0;
    total.correct=0;

    %% DELPHINE, here the input vector is turned into fields of 'param' for
    %% convenience.
    %%Extract trial parameters into parameter vector...
    [param,opt_vec_scaled]=extract_optimization_parameters(param_in,opt_vec,'assign',run_options.opt_param_chc);
    %keyboard;
    param.morph
   
    for Idate=1:length(goodFile)
        
        for Istation=1:length(goodFile{Idate})
            

            %Start loop through files
            strr='abcdefghijkl';
            Iletter=findstr(goodDASAR{Idate}{Istation}(5),upper(strr));
            if isempty(Iletter)
                continue
            end
    

            [best_calls,debug_ctimes,debug_durations,raw_detections,interval_detections]=process_one_unit(goodFile{Idate}{Istation},goodName{Idate}{Istation},param,run_options);

            if isempty(debug_durations)
                disp(sprintf('Date: %s Station: %s had process_one_unit failure',date_str_new{Idate},strr(Iletter)));
               continue 
            end
            [Imiss,Ifalse,miss_diff,Itrue,SNR_miss,SNR_true]=evaluate_stations(best_calls.ctime,best_calls.duration,manual.individual{Idate}.ctime(:,Iletter),...
                manual.individual{Idate}.duration(:,Iletter),manual.individual{Idate}.stndb(:,Iletter),run_options);

        %Weight missed results by SNR--lower the SNR, the less the penalty
        wght=SNR_miss/param.morph.SNRmin;
        wght(wght>1)=1;
        
        wght2=SNR_true/param.morph.SNRmin;
        wght2(wght2>1)=1;
       
%         total.missed=total.missed+sum(wght);
%         total.correct=total.correct+sum(wght2);
%         
         total.missed=total.missed+length(wght);
         total.correct=total.correct+length(wght2);
%         

        total.false=total.false+length(Ifalse);

        end
    end


    final.missed_call_fraction=total.missed./(total.missed+total.correct);
    final.false_alarm_ratio=total.false./(total.missed+total.correct);
    if total.false==0| ~isnan(final.missed_call_fraction)
        %%%%%%%%%
        %%DELPHINE, an example of an optimization function.  I compute the
        %%missed fraction and false alarm fraction and convert into a
        %%single number 'fit'.  The lower the value of 'fit' the better the
        %%result.  I have to use the 'max' command and 'fit_goal' variables
        %% to ensure that the optimization program does not force the
        %% missed fraction to zero while making the false alarm rate huge.
        %% The 'fit_weight' vector assigns relative importance to the
        %% missed call fraction and false alarm fraction

        fit=fit_weight(1)*(max([final.missed_call_fraction-fit_goal(1) 0]).^2)+fit_weight(2)*(max([final.false_alarm_ratio-fit_goal(2) 0]).^2);
        disp(sprintf('Final Input vector: %s \n missed frac: %6.4f false alarm ratio: %6.4f, current fit: %8.6f \n\n', ...
            mat2str(opt_vec_scaled,3),final.missed_call_fraction,final.false_alarm_ratio,fit));
    else
        fit=fit_weight(1)*(max([1 0]).^2)+fit_weight(2)*(max([100-fit_goal(2) 0]).^2);

    end
catch
    disp('Program crashed');
    %keyboard;
    myerror=lasterror
    fclose('all');
    !rm *.detsum *.snips

    savestr=sprintf('Dump%s',datestr(now,30));
    save(savestr);
    fit=fit_weight(1)*(max([1 0]).^2)+fit_weight(2)*(max([100-fit_goal(2) 0]).^2)

end
end

function [Imiss,Ifalse,miss_diff,Itrue,SNR_miss,SNR_true]=evaluate_stations(auto_times,auto_durations,correct_times,correct_durations, correct_SNR,run_options)

%auto_times=station.ctime;
Igoodm=find(correct_times>0&correct_times<=run_options.max_ctime);
correct_times=correct_times(Igoodm);
correct_durations=correct_durations(Igoodm);


Igooda=find(auto_times>0&auto_times<=run_options.max_ctime);
auto_times=auto_times(Igooda);
auto_durations=auto_durations(Igooda);

run_options.tol=0.25;
[Ifalse,Imiss,Iauto_match,Imatch,miss_diff]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,run_options.tol);
%         else
%            [Ifalse,Imiss,Iauto_match,Imatch,miss_diff]=find_similar_elements(auto_times, correct_times,run_options.tol,auto_fmean,correct_fmean);
%
%         end
%false_alarms=Igooda(Ix_nomatch);
%missed_calls=Igoodm(Iy_nomatch);
%matched_calls=Igoodm(Iy_match);
%   disp(sprintf('\nDate: %s  Station: %s number of true calls: %i',datte ,station_name,length(correct_times)));
disp(sprintf(' number of true calls: %i',length(correct_times)));
disp(sprintf('Station results: false_alarms: %i false alarm fraction: %6.4f matched calls %i, missed calls: %i, missed fraction: %6.4f \n', ...
    length(Ifalse),length(Ifalse)/length(correct_times), ...
    length(Imatch),length(Imiss),length(Imiss)/length(correct_times)));
% disp(sprintf('Number of succesful manual calls detected: %i %i\n', ...
%    length(Ix_match), length(Iy_match)));


%Export corresponding indicies to automated stations...
%Ix_match=Igooda(Iauto_match);
Imiss=Igoodm(Imiss);
SNR_miss=correct_SNR(Imiss);
Ifalse=Igooda(Ifalse);
Itrue=Igoodm(Imatch);
SNR_true=correct_SNR(Itrue);

end
