
function stats=find_corresponding_matches(base_indiv,compare_indiv,run_options,bufferTol,min_ctime,max_ctime,stations)
%function stats=find_corresponding_matches(base_indiv,compare_indiv,run_options,bufferTol,min_ctime,max_ctime,stations)
%          Input:
%          run_options.match_tol:  fractional overlap required for a match.  e.g. 0.25 means 25% ovlap in time
%          run_options.strict_eval_tolerance:  If one, manual and auto detections must
%               overlap 50% in time AND frequency.  If two, compare manual
%               and auto detections via time tolerances only, with no overlap
%               requirement.  If three, try to match
%               by location...
%           bufferTol:  estimated uncertainty in manual time selection in seconds...only used for 'time_only'
%           criteria
% Compare manual vs. automated linkages...
% stats:  structure, with all fields (except one) same length as number of base_indiv detections:
%           Nhits: Number of stations(detections) used in the base_indiv location.
%           Ncompare_matches: The number of time/frequency matches from a base linkage
%               that are present in
%               any compare linkage.  Anwsers whether a given base
%               station detection  has been linked into any compare detection. 
%           Nstation_match: The number of base detections present in any
%               compare  station.  Answers whether any detections were available for
%               linking.
%           Nstations_index_match:  a [Nmanual Ndasr]  array that contains
%               indicies of stations that match given manual call at that
%               DASAR, thus determining whether automated detections were
%               available to link. If the number has a fraction component,
%               then multiple station matches occurred, and 100*fraction is
%               the number of matches.  If '-1' manual detection not present at station.
%               If '0' manual detection exists at station, but not automated match found.
%           Nsplit_links: The number of compare locations that contain
%               portions of a given base location.  Should be one,
%               anything higher indicates a linkage problem.  Different
%               than Ncompare_matches in that the former counts total
%               detections matching the bsase, while Nsplit_links lists
%               the number of compare links that have contributed these
%               detections.
%           Ilocs_best_match:  The index of the automated location that
%               contains the largest number of matches with a given
%               manual result.
%           Nlongest_link: The maximum number of automated detections present
%               in the manual location given by 'Ilocs_best_match'
%
%           match_flag:  Contains '1' if associated manual detection is
%               matched, otherwise zero
%


%%Set fields for output stats information

if ~exist('stations')
    stations=[];
end

Nbase=size(base_indiv.ctime,1);  %Number of base detections
Nstations=size(base_indiv.ctime,2);

names={'Nstations_index_match','Nstation_match','Ncompare_matches','Nlongest_link', ...
    'Nsplit_links','Ilocs_best_match','match_flag','Nhits'};
names=names(3:end);
if ~isempty(stations)

    stats.Nstations_index_match=-1*ones(Nbase,Nstations);
end


%Initialize output statistics
for I=1:length(names)
    stats.(names{I})=zeros(Nbase,1);
end



%Process stations efficiently...

if ~isempty(stations)
    for J=1:Nstations
        Igood=find(base_indiv.ctime(:,J)>0);
        switch run_options.strict_eval_tolerance
            case 1  %strict, match frequecy, time, and duration
                %find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, ...
                %min_ctime,max_ctime, match_tol,auto_features,correct_features)
                [Imiss,Ifalse,Iman_pass,Itrue]=find_similar_elements_ovlap_timelimit(stations(J).ctime_min,stations(J).Totalduration, base_indiv.ctime(Igood,J), ...
                    base_indiv.duration(Igood,J),0,max_ctime,run_options.match_tol,...
                    [stations(J).Totalfmin; stations(J).Totalfmax],[base_indiv.flo(Igood,J) base_indiv.fhi(Igood,J)]');
            case 2  %match time only
                %bufferTol=3;
                [Imiss,Ifalse,Iman_pass,Itrue]=find_similar_elements_ovlap_timelimit(stations(J).ctime_min,stations(J).Totalduration, base_indiv.ctime(Igood,J)-bufferTol, ...
                    base_indiv.duration(Igood,J)+2*bufferTol,0,max_ctime,run_options.match_tol);
        end
        %length(Iman_pass)+length(Imiss)=length(Igood)
        stats.Nstation_manual(J)=length(Igood);
        stats.Nstation_match(J)=length(Iman_pass);
        stats.Nstation_miss(J)=length(Imiss);
        stats.Nstation_false(J)=length(Ifalse);
        stats.Nstations_index_match(Igood(Iman_pass),J)=Itrue';
        
    end
end

for I=1:Nbase  %For every base localization
    clear match_index
    if rem(I,100)==0, disp(sprintf('%6.1f percent done',100*I/Nbase));end
    
    Ilinks=[];
    Ivalid=find(base_indiv.ctime(I,:)>0);
    stats.Nhits(I)=length(Ivalid);
    stats.Nlongest_link(I)=-1;
    
    Icount=0;
    if size(Ivalid,1)>1
        Ivalid=Ivalid';
    end
     %stats.Nstations_index_match{I}=cell(1,7);
    
   
    
     %For every station used in base location
     for J=Ivalid
         
%          if ~isempty(stations)
%              switch run_options.strict_eval_tolerance
%                  case 1  %strict, match frequecy, time, and duration
%                      %find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, ...
%     %min_ctime,max_ctime, match_tol,auto_features,correct_features)
%                      [Imiss,Ifalse,Iman_pass,Itrue]=find_similar_elements_ovlap_timelimit(stations(J).ctime_min,stations(J).Totalduration, base_indiv.ctime(I,J), ...
%                          base_indiv.duration(I,J),0,max_ctime,run_options.match_tol,...
%                          [stations(J).Totalfmin; stations(J).Totalfmax],[base_indiv.flo(I,J) base_indiv.fhi(I,J)]');
%                  case 2  %match time only
%                      %bufferTol=3;
%                      [Imiss,Ifalse,Iman_pass,Itrue]=find_similar_elements_ovlap_timelimit(stations(J).ctime_min,stations(J).Totalduration, base_indiv.ctime(I,J)-bufferTol, ...
%                          base_indiv.duration(I,J)+2*bufferTol,0,max_ctime,run_options.match_tol);
%              end
%              
%              if ~isempty(Itrue)
%                  stats.Nstation_match(I)=stats.Nstation_match(I)+1;
%                 
%                  stats.Nstations_index_match(I,J)=min(Itrue)+length(Itrue)/100;
%              else
%                  stats.Nstations_index_match(I,J)=0;
%             
%              end
%              
%              
%          end
         
         
         
         %disp('Final detection comparison with frequency feature matching')
         switch run_options.strict_eval_tolerance
             case 1  %strict, match frequecy, time, and duration
                 %[Imiss,Ifalse,Itrue_man,Itrue_auto]=find_similar_elements_ovlap_timelimit(auto_times,auto_durations,correct_times,correct_durations, max_ctime, auto_features,correct_features)
                 
                 [Imiss,Ifalse,Iman_pass,Itrue]=find_similar_elements_ovlap_timelimit(compare_indiv.ctime(:,J),compare_indiv.duration(:,J), base_indiv.ctime(I,J), ...
                     base_indiv.duration(I,J),min_ctime,max_ctime,run_options.match_tol,[compare_indiv.flo(:,J) compare_indiv.fhi(:,J)]',[base_indiv.flo(I,J) base_indiv.fhi(I,J)]');
             case 2  %match time only
                 %bufferTol=3;
                 [Imiss,Ifalse,Iman_pass,Itrue]=find_similar_elements_ovlap_timelimit(compare_indiv.ctime(:,J),compare_indiv.duration(:,J), base_indiv.ctime(I,J)-bufferTol, ...
                     base_indiv.duration(I,J)+2*bufferTol,min_ctime,max_ctime,run_options.match_tol);
                 
         end
         
         
         if ~isempty(Itrue)
             stats.Ncompare_matches(I)=stats.Ncompare_matches(I)+1;
             
         end
         
         %Ilinks contain indicies of locations that match base result on
         %at least one DASAR.  The union operation removes double-counting
         
         Ilinks=union(Ilinks,Itrue);
         Icount=Icount+1;
         match_index{Icount}=Itrue;  %match_index{Istation} lists matching stations for a given base detction
         
     end
    
    
    %Track all locations that match base detection at at least one
    %   point...
    %Istation_match=union(Istation_match,Ilinks);
    
    stats.Nsplit_links(I)=length(Ilinks);
    
    
    
    if size(Ilinks,1)>1
        Ilinks=Ilinks';
    end
    
    %Go through all compare locations that have a partial match
    for J=Ilinks
        Ntemp=0;
        for K=1:length(match_index)
            Ntemp=Ntemp+length(find(match_index{K}==J));
        end
        stats.Nlongest_link(I)=max([stats.Nlongest_link(I) Ntemp]);
        if Ntemp==stats.Nlongest_link(I)
            stats.Ilocs_best_match(I)=J;
        end
    end
    
    %If manual and automated share two DASAR detections, flag as a match
    if stats.Nlongest_link(I)>1
        stats.match_flag(I)=1;
    end
    
   
end

%stats.auto_match_flag(Istation_match)=1;

end

   
