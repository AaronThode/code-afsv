%%%%%%%%%evaluate_miss_fraction.m%%%%%%%%%%%%%%%%%%%%%%
%function [fraction_pass,Ipass_out,Imiss,Iloc_match]=evaluate_miss_fraction(manind,autolocal,auto_ctimes,tol, ...
%  min_station, max_auto_ctime,Idebug);
%
% Input:
%       manind: structure of individual manual detections--
%       autolocal: localization structure of automated detections
%       auto_ctimes:  ctimes that correspond to autolocal... [Ncalls x Nstations]
%       tol: three element vector with time tolerance, fmin tolerance, fmax
%           tolerance
%       compare_func_str: "absolute_error" or "relative_error"
%       frac_cutoff: percent of stations that must be present to pass
%       max_auto_ctime: Maximum ctime in auto_ctimes
%       excluded_stations:  string of characters representing stations not
%           to use in evaluations (e.g. Station  'e' of Site 4 has bad time
%           info)
%    
%       Idebug: if not zero, print debug output
% Output:
%       fraction_pass.ctime: For a given manual locations, the fraction of
%           manually-detected stations that match arrival times in automated
%           locations.
%       Npass.ctime: For a given manual location, max number of automated
%           stations that match..
%       fraction_pass.feature:  fraction of manual-detected station that
%           also pass feature tolerance criteria.
%       Ipass_out.ctime:  indicies of manual detections that correspond to
%           locations with greater than frac_cutoff fraction time match
%       Ipass_out.feature: indicies of manual detection that correspond to
%           locations with greater than frac_cutoff time and frequency match
%       Imiss: Complementary indicies to Ipass_out
%       Iloc_match: indicies of locations that correspond to
%               Ipass_out.feature...
function [Ipass_out,Imiss,Iloc_match]=evaluate_miss_fraction(manind,autolocal,auto_ctimes,tol, ...
    min_station,max_auto_ctime,Idebug)

disp(' ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%% evaluate_miss_fraction.m%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%strr_org='abcdefghijklmnop';
strr_org=manind.DASARstr;


Nlocs=length(autolocal);
disp(sprintf(' %i Manual results , %i automated localizations.',size(manind.ctime,1),Nlocs))
compare_func_str='absolute_error';
compare=str2func(compare_func_str);
%Nstations=size(manind.duration,2);

%Trim manual results if a complete day has not been processed in automatic data.
Ipass=find(max(manind.ctime')<=max_auto_ctime);

names={'ctime','flo','fhi'};
for I=1:length(names),
    manual.(names{I})=manind.(names{I})(Ipass,:);
    %manual.(names{I})(:,Ireject)=-1;
end

Iloc_match=zeros(1,length(Ipass));

Ncalls_manual=length(Ipass);
Ncalls_auto=length(auto_ctimes);

disp(sprintf(' %i Manual results less than %s, %i automated localizations.',length(Ipass),ctime2str(max_auto_ctime),Ncalls_auto));


%Force non-detected station times to be infinite).
Ibad=find(auto_ctimes==0);
auto_ctimes(Ibad)=Inf;

for I=1:Ncalls_manual,
    
    cdiff=abs(ones(Ncalls_auto,1)*manual.ctime(I,:)-auto_ctimes);
    % check=min(cdiff');
    Npass=sum((cdiff<=tol(1))');

    [maxpass(I),Iloc_match(I)]=max(Npass);
    if maxpass<min_station,
        Iloc_match(I)=-1;
    end

    if Idebug>0,
        %disp(' ');disp(sprintf('Call %i manually detected on stations %s',I,strr_org(Iman_pass)));
        disp(sprintf('Best time match for Manual call %i at %ith location, %i matches found ', ...
            I,Iloc_match(I),maxpass(I)));
    end

    %%Time for feature comparison!
    %Ipass_ctime=find((Npass./Ngood_manual_stations(I))>frac_cutoff);
    Ngood(I)=length(find(Npass>=min_station));

  
    

end



Ipass_out.count=find(maxpass>=min_station);
Imiss.count=setdiff(1:Ncalls_manual,Ipass_out.count);
Npass_struct.count=length(Ipass_out.count);

disp(sprintf('%i out of %i, or %6.2f percent of manual locations have more than %i matches with an automated localization', ...
    Npass_struct.count,Ncalls_manual,round(100*Npass_struct.count/Ncalls_manual),min_station));

disp(sprintf('%i locations remain, or false/true ratio of %6.2f ', ...
    Ncalls_auto-Npass_struct.count,((Ncalls_auto-Npass_struct.count)/Npass_struct.count)));


disp('%%%%%%%%%%%%%%%%');

