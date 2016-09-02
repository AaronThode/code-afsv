%%%%evaluate_linking_false.m%%%%
% function [whale,other,Ntrue]=evaluate_linking(param,manual_locations,auto_locations,auto_ctimes,date_str,goodDASAR,run_options)
%    Compare manual detections with automated, in terms of detection times
%    at stations.
%     Inputs:
%           manual_locations:  manual manual_locations results
%   Outputs (all same length as number of manual detections):
%           Nstations_match: The number of time/frequency matches present in
%           locations.
%           Nmanual_match:  a {Nauto}{Ndasr} cell array that contains
%               indicies of stations that match given manual call
%           (before linkage)
%           Nstation_extra: The maximum number of linked DASARS
%           Nsplit_links: The number of automated locations that contain
%               parts of the manual detection
%
function [false_alarm]=evaluate_linking_false(param,manual_locations,auto_locations,auto_ctimes,date_str,goodDASAR,run_options)

%Set maximum time to evaluate--are we evaluating only a partial day?
false_alarm.Nmanual_match=[];
false_alarm.manual_match_index=[];
false_alarm.Nstation_miss=[];
false_alarm.Nstation_extra=[];
false_alarm.Nauto_stations=[];
false_alarm.Nsplit_links=[];

if isempty(manual_locations.ctime)|isempty(auto_ctimes)
    disp('empty data')
    return
end
first_time=auto_ctimes(1,manual_locations.ctime(1,:)>0);

if ~isfield(param.energy,'nsamples')||param.energy.nsamples==0,

    max_ctime=Inf;
else
    max_ctime=min(first_time)+param.energy.nsamples/param.Fs;
end


Nmanual=length(manual_locations.maxSNR);
Nauto=size(auto_ctimes,1);

Nstations=length(auto_locations{1}.ctime_min);
if size(manual_locations.ctime,2)~=Nstations
    disp('manual and auto dasar numbers don''t match');

end

%Unpack automated information..
auto_durations=zeros(Nauto,Nstations);
auto_flo=auto_durations;
auto_fhi=auto_durations;
for I=1:Nauto
    auto_durations(I,:)=auto_locations{I}.Totalduration';
    auto_flo(I,:)=auto_locations{I}.Totalfmin';
    auto_fhi(I,:)=auto_locations{I}.Totalfmax';

end

false_alarm.Nmanual_match=zeros(Nauto,1);
false_alarm.Nstation_extra=false_alarm.Nmanual_match;
false_alarm.Nstation_miss=Inf+false_alarm.Nmanual_match;
false_alarm.Nauto_stations=false_alarm.Nmanual_match;
false_alarm.Nsplit_links=false_alarm.Nmanual_match;
false_alarm.manual_match_index=false_alarm.Nmanual_match;
for I=1:Nauto,
    clear match_index

    Ilinks=[];
    Ivalid=find(auto_ctimes(I,:)>0);
    false_alarm.Nauto_stations(I)=length(Ivalid);
    Icount=0;
    if size(Ivalid,1)>1
        Ivalid=Ivalid';
    end
    for J=Ivalid
        %disp('Final detection comparison with frequency feature matching')
        [Imiss,Ifalse,Iman_pass,Itrue]=evaluate_stations( manual_locations.ctime(:,J),  manual_locations.duration(:,J), ...
            auto_ctimes(I,J),auto_durations(I,J),max_ctime,date_str, goodDASAR{J}, ...
            [manual_locations.flo(:,J)'; manual_locations.fhi(:,J)'],[auto_flo(I,J); auto_fhi(I,J)]);

        if isempty(Itrue)
            false_alarm.Nstation_extra(I)=false_alarm.Nstation_extra(I)+1;
        end

        Ilinks=union(Ilinks,Itrue);
        Icount=Icount+1;
        match_index{Icount}=Itrue;

    end
    
    false_alarm.Nsplit_links(I)=length(Ilinks);
    if size(Ilinks,1)>1
        Ilinks=Ilinks';
    end
    
    for J=Ilinks  %For every manual detection associated with this automated call
        Ntemp=0;
        for K=1:length(match_index)
            Ntemp=Ntemp+length(find(match_index{K}==J));
        end
        if Ntemp>=false_alarm.Nmanual_match(I)
            false_alarm.Nmanual_match(I)=Ntemp;
            false_alarm.manual_match_index(I)=J;
        end
    end
    if ~isempty(Ilinks)
        Nmanual_loc=length(find(manual_locations.ctime(false_alarm.manual_match_index(I),:)>0));
        false_alarm.Nstation_miss(I)=Nmanual_loc-false_alarm.Nmanual_match(I);
    end

end




    function [Imiss,Ifalse,Itrue_man,Itrue_auto]=evaluate_stations(auto_times,auto_durations,correct_times,correct_durations, max_ctime,datte, station_name, auto_features,correct_features)
        % auto and correct features have a [Nfeature x Ntime] matrix

        %auto_times=station.ctime;
        Igoodm=find(correct_times>0&correct_times<=max_ctime);
        correct_times=correct_times(Igoodm);
        correct_durations=correct_durations(Igoodm);
        if exist('correct_features')
            correct_features=correct_features(:,Igoodm);
        end


        Igooda=find(auto_times>0&auto_times<=max_ctime);
        auto_times=auto_times(Igooda);
        auto_durations=auto_durations(Igooda);
        if exist('auto_features')
            auto_features=auto_features(:,Igooda);
        end

        run_options.tol=0.5;
        if exist('auto_features')
            [Ifalse,Imiss,Iauto_match,Imatch,miss_diff,jun2,Iauto_match2]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,run_options.tol,auto_features,correct_features);

        else
            [Ifalse,Imiss,Iauto_match,Imatch,miss_diff,jun2,Iauto_match2]=find_similar_elements_ovlap(auto_times, auto_durations, correct_times,correct_durations,run_options.tol);

        end
        %if ~isempty(Iauto_match2{1})
        %     keyboard;
        % end
        if ~isempty(Iauto_match2)
            if size(Iauto_match2{1},1)>1
                Iauto_match2{1}=Iauto_match2{1}';
            end

            Iauto_match=sort([Iauto_match Iauto_match2{1}]);
        end
        %disp(sprintf('Station %s on %s, results: false_alarms: %i false alarm fraction: %6.4f matched calls %i, missed calls: %i, missed fraction: %6.4f \n', ...
        %     station_name, datte,length(Ifalse),length(Ifalse)/length(correct_times), ...
        %     length(Imatch),length(Imiss),length(Imiss)/length(correct_times)));

        %Export corresponding indicies to automated stations...
        %Ix_match=Igooda(Iauto_match);
        Imiss=Igoodm(Imiss);
        Itrue_man=Igoodm(Imatch);
        Ifalse=Igooda(Ifalse);
        Itrue_auto=Igooda(Iauto_match);

    end

end

%%%%%%%%%%find_similar_elements_ovlap.m%%%%%%%%
%function [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match,Iy_match_diff,Iy_nomatch_diff]=find_similar_elements_ovlap(x,dx,y,dy,ovlap_tol,xfeature,yfeature);
% x, y:  two vectors with times in same units
% tol: how close x(I) and y(J) have to be to be considered a "match".
% xfeature,yfeature: Optional features of x and y matricies, used to select
%   best x if multiple values within tol of y. [Nfeature x Nsamples]
% Output:
%  Ix_nomatch:  Indicies of x that have no counterpart in y
%  Ix_match:  Indicies of x that match y
%  Iy_match:  Indicies of y that match x, thus
%       x(Ix_match(k))=y(Iy_match(k));
%  Iy_match:  Indicies of y that have no counterpart in x
%  Ix_nomatch_diff:  how far this element of x was to nearest y element
%  Ix_redundant=other elements of x that match this element of y.
%
%  NOTE: the sum of lengths of Ix_match and Ix_nomatch might not equal
%  length of x, if more than one x element lies within tolerance of a given
%  y.

function [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match,Iy_nomatch_diff,Iy_match_diff,Ix_redundant]=find_similar_elements_ovlap(x,dx,y,dy,ovlap_tol,xfeature,yfeature)
tol=5;

Nel=10000;
Ix_match=[];
Ix_nomatch=[];
Ix_nomatch_diff=[];
Iy_match=[];
Iy_nomatch=[];
Iy_nomatch_diff=[];
Iy_match_diff=[];
Ix_redundant=[];

Nx=[];
if exist('xfeature')
    Nfeatures=size(xfeature,1);
end

if isempty(x)
    Iy_nomatch=1:length(y);
    return
end

if isempty(y)
    Ix_nomatch=1:length(x);
    return;
end
Igood=1:length(x);
if length(x)>1 & length(y)>1
    Igood=find((x<=max(y))&(x>=min(y)));
    disp(sprintf('x trimmed for time range: x was %i long, now %i points remain',length(x),length(Igood)));

    %keyboard;
    x=x(Igood);
    dx=dx(Igood);
    if exist('xfeature')
        xfeature=xfeature(:,Igood);
    end

    if isempty(x)
        disp('No x points remain');
        return;
    end
end

xref=x(1);
Nx=length(x);
Ny=length(y);

%March through reference times
for I=1:Ny,  %these are the "true" manual values, for example
    Ix_redundant{I}=[];

    fail_flag=0;
    %compute overlap
    Ixp=find((abs(y(I)-x)<=2*tol)&(abs(y(I)+dy(I)-x-dx)<=2*tol));  %Automated signals that that lie within time tolerance
    if isempty(Ixp)
        fail_flag=1;
        ovlap=-1;
        Imin=[];
    else  %potential matches exist
        ovlap=-1*ones(size(Ixp));
        if exist('xfeature')
            feature_ovlap=-1*ones(size(Ixp));
        end

        for J=1:length(Ixp),
            index=Ixp(J);
            tovlap=min([y(I)+dy(I) x(index)+dx(index)])-max([y(I) x(index)]);  %Time of overlap in sec
            ovlap(J)=(tovlap./min([dx(index) dy(I)]));  %Negative if no overlap
            if exist('xfeature')  %If extra features provided
                fovlap=min([yfeature(2,I) xfeature(2,index)])-max([yfeature(1,I) xfeature(1,index)]);
                feature_ovlap(J)=(fovlap./min([diff(yfeature(:,I)) diff(xfeature(:,index))]));  %Negative if no overlap

            end
        end

        if exist('xfeature')
            Imin=find((ovlap>=ovlap_tol)&(feature_ovlap>=ovlap_tol));
        else
            Imin=find(ovlap>=ovlap_tol);
        end
    end

    if isempty(Imin)|fail_flag==1,
        Iy_nomatch=[Iy_nomatch I];
        Iy_nomatch_diff=[Iy_nomatch_diff max(ovlap)];
        continue;
    end

    Iy_match=[Iy_match I];

    if length(Imin)==1,
        Ix_match=[Ix_match (Ixp(Imin))];
        Iy_match_diff=[Iy_match_diff ovlap(Imin)];

    else  %multiple x candidates

        [junk,Ibest]=max(ovlap(Imin));

        Ix_match=[Ix_match Ixp(Imin(Ibest))];
        %Should we remove other candidates from "false" category, if we have feature information?

        if exist('xfeature')
            Ix_redundant{I}=Ixp(Imin(setdiff(1:length(Imin),Ibest)));
        else
            Ix_redundant{I}=Ixp(Imin(setdiff(1:length(Imin),Ibest)));
        end

    end

end


Ix_match=sort(Ix_match);
if isempty(Ix_match),
    Ix_nomatch=1:length(x);
    return;
end


%Now look for false alarms in automated data...
Ix_nomatch=setdiff(1:length(x),Ix_match);

if Ny==1
    return
end
%Sanity check for false detections-do some actually qualify?
Nfunny=length(Ix_nomatch);

Ix_nomatch_ovlap=zeros(1,Nfunny);
for I=1:Nfunny,  %For each false detection
    indexx=Ix_nomatch(I);
    Iyp=find((abs(y-x(indexx))<=2*tol)&abs(y+dy-x(indexx)-dx(indexx))<=2*tol);
    if isempty(Iyp)
        Ix_nomatch_ovlap(I)=-1;
    else

        for J=1:length(Iyp),
            indexy=Iyp(J);
            tovlap=min([y(indexy)+dy(indexy) x(indexx)+dx(indexx)])-max([y(indexy) x(indexx)]);
            Ix_nomatch_ovlap(I)=max([Ix_nomatch_ovlap(I) (tovlap./min([dx(indexx) dy(indexy)]))]);  %Negative if no overlap
            if Ix_nomatch_ovlap(I)>ovlap_tol
                if exist('xfeature')  %If extra features provided
                    fovlap=min([yfeature(2,indexy) xfeature(2,indexx)])-max([yfeature(1,indexy) xfeature(1,indexx)]);
                    feature_ovlap=(fovlap./min([diff(yfeature(:,indexy)) diff(xfeature(:,indexx))]));  %Negative if no overlap
                    if feature_ovlap<ovlap_tol
                        Ix_nomatch_ovlap(I)=-1;
                    end

                end

            end

        end
    end
end
Itruly_bad=find(Ix_nomatch_ovlap<ovlap_tol);
Ix_nomatch=Ix_nomatch(Itruly_bad);
Ix_nomatch_ovlap=Ix_nomatch_ovlap(Itruly_bad);
%disp(sprintf('Out of %i false matches, %i/ %i were actually overlapped more than %6.2f  with a true value and not included in final false result', Nfunny,length(Ix_redundant), ...
%    Nfunny-length(Itruly_bad),ovlap_tol));

Ix_match=Igood(Ix_match);
Ix_nomatch=Igood(Ix_nomatch);

%disp(sprintf('Detected calls: %i  plus missed calls: %i = %i should be %i',length(Iy_match),length(Iy_nomatch),length(Iy_match)+length(Iy_nomatch),length(y)));
disp('');
% Iskip=find(diff(Iy_match)>1);
% for I=1:length(Iskip),
%     Istart=Iy_match(Iskip(I))+1;
%     Iend=Iy_match(Iskip(I)+1)-1;
%     Iy_nomatch=[Iy_nomatch Istart:Iend];
% end

end