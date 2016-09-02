

%%%load_automated_results.m%%%
%[locations,locations_ctime,stations,raw_stations,automated_dir,auto_param,error_ellipse]=load_automated_results_crosschecked(keyword, ...
%       detectiondir, goodName, date_range, max_ctime, run_options, d_coords)
%   Input:
%        keyword eg. 'Site_4_short.morph.crosschecked'
%               Format is 'scenario.algorithm.process_stage', divided by periods.
%               algorithm:  %'contour','morph','both','none':
%               process_stage:  String of concatenated phrases like
%                   'RawCrosscheckedRandomClassified'.
%        detectiondir: string with location of 'Processed' directory
%        goodName: cell matrix [Ndays x 1] each cell itself a cell
%                   array containing name strings like 'D07s4e', e.g., stations
%                   that are recording on that particular day.
%        date_range: [Ndays x1] matrix of strings of dates that cover desired
%                       time period.  Each row is a different string for a distinct
%                       day
%        max_ctime:  Maximum ctime allowed for an automated detection
%        run_options:
%              .success_only:  If 1, then keep only locations that were succesfully tracked and lie within
%                           run_options.center_dist_limit km of the centerline.
%              .localization_alg:  'keep_stored_result'--don't recompute...
%              .auto_min_stations[min max]:  minimum and maximum number of stations in localization
%                           permitted
%              .center_dist_limit: used if .success_only activited--see above
%              .Ikeep_only:  If 1, adjust locations_ctime so that only DASARS represented in Ikeep are used.
%                       Only matters if locations were 'pruned' by setting a minimum Huber weight > 0.
%        d_coords:  UTM coords [east west] of DASARS, used only if run_options.success_only=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Output:
%       locations: see cross_channel_match.m
%       locations_ctime:  ctimes of located signals in locations
%       stations: individual results of each DASAR--see
%               extract_feature_vector_noshape.m
%       automated_dir: string with location of automated data...
%       auto_param: structure of parameters used to produce results
%       auto_goodName: cell arry of names of DASAR array


function [locations,locations_ctime,stations,raw_stations,automated_dir,auto_param,auto_goodName,error_ellipse]=load_automated_results(keyword, ...
    detectiondir, goodName, date_range,max_ctime,run_options,d_coords)

error_ellipse=[];
locations{1}=[];
locations_ctime{1}=[];
stations=[];
raw_stations=[];
automated_dir=[];
auto_param=[];
auto_goodName=[];

if isempty(goodName{1})
    return
end

% fnames_auto=dir([automated_dir '/*' detect_str '.mat']);
for I=1:size(date_range,1),
    
    if isempty(date_range(I,:)),
        continue;
    end
    
    %%Get appropriate search string for files, given keyword and date_range
    %%desired
    %[site_number,search_str]=lookup_files(I,date_range,goodName,keyword);
    site_number=goodName{I}{1}(2);  %Site number
    search_str=['*S' site_number '*' date_range(I,:) '*' keyword.stage '*.mat'] ;
    
    
    automated_dir=[detectiondir '/Site_0' site_number '/' keyword.scenario '/' keyword.algorithm];
    disp(sprintf('loading automated detection results from %s',automated_dir));
    if isempty(automated_dir)
        disp('automated directory results not found in load_automated_results.m');
    end
    %Find file name that matches dates and site desired..
    fname_auto=dir([automated_dir '/' search_str ]);
    if length(fname_auto)>1
        disp(sprintf('More than one automated detection file for %s',search_str));
        for JJJ=1:length(fname_auto)
            disp(fname_auto(JJJ).name);
        end
        JJJ=input('Which file?');
        fname_auto=fname_auto(JJJ);
    elseif length(fname_auto)==0
        disp('No files found...');
        continue
    end
    
    %%%Load automated date from file..
    fname_auto=[automated_dir '/' fname_auto.name];
    disp(['loading ' fname_auto]);
    auto=load(fname_auto);
    
    %%Safety check that DASARS in package match desired DASARS from manual
    %%data...
    Imatch=zeros(1,length(goodName{I}));
    for J=1:length(goodName{I}),
        if isfield(auto,'goodName')&&isempty(strfind(goodName{I}{J},auto.goodName{J}))
            fprintf('Warning: requested DASAR %s  doesn''t match automated DASAR %s, rearranging.\n',goodName{I}{J},auto.goodName{J});
            
            %%Rearrange automated structure to reflect input DASAR structure.
            for Idasar=1:length(auto.goodName)
                if ~isempty(strfind(goodName{I}{J},auto.goodName{Idasar}))
                    Imatch(J)=Idasar;
                    break
                end
                
            end
            
            
        else
            Imatch(J)=J;
        end
        
    end
    
    
    %%If any stations need to be rearranged, rearrange all columns...
    Ishould=1:length(goodName{I});
    
    if sum(abs(Ishould-Imatch))>0
        disp(sprintf('New DASAR order is %s',mat2str(Imatch)));
        
        auto.goodName=auto.goodName(Imatch);
        auto.goodFile=auto.goodFile(Imatch);
        auto.station=auto.station(Imatch);
        auto.raw_station=auto.raw_station(Imatch);
        auto.locations_ctime=auto.locations_ctime(:,Imatch);
        
        for Iloc=1:length(auto.locations)
            fields=fieldnames(auto.locations{Iloc});
            for Ifld=1:length(fields)
                switch fields{Ifld}
                    
                    case {'position','equalization','credibility'}
                    otherwise
                        auto.locations{Iloc}.(fields{Ifld})=auto.locations{Iloc}.(fields{Ifld})(Imatch);
                end
            end
            
        end
    end
    
    
    %%Trim stations and locations to maximum desired time
    trim_max_ctime;
    stations{I}=auto.station;
    raw_stations{I}=auto.raw_station;
    auto_goodName{I}=auto.goodName;
    auto_param=auto.param;
    
    
    
    %%%%If no locations (training set), continue next loop
    if ~isfield(auto,'locations')||isempty(auto.locations)
        continue
    end
    
    
    
    %First, option to recompute localization...
    %run_options.localization_alg={'keep stored result','Huber','medianHuber','HuberCalibratedKappa','Andrews'};
    %Iloc_chc=menu('Select localization method:',run_options.localization_alg);
    if ~strcmp(run_options.localization_alg,'keep stored result')
        disp('Recomputing automated positions...');
        auto.locations=compute_position(auto.locations,[],auto.param,auto.Icase,auto.Isite,run_options);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Ensure only locations that have minimum number of stations and less than another sum are
    %tracked...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Igood=[];
    for KK=1:length(auto.locations),
        Ngood(KK)=length(find(auto.locations{KK}.station_indicies>0));
        if Ngood(KK)>=run_options.auto_min_stations(1)&&Ngood(KK)<=run_options.auto_min_stations(2)
            Igood=[Igood KK];
        end
    end
    
    disp(sprintf('Out of %i automated locations, %i linked between %i and %i DASARS', ...
        length(auto.locations),length(Igood),run_options.auto_min_stations(1),run_options.auto_min_stations(2)));
    locations{I}=auto.locations(Igood);
    locations_ctime{I}=auto.locations_ctime(Igood,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%If desired, use only succesfully located calls%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if run_options.auto_success_only==1
        if isempty(run_options.center_dist_limit)
            center_dist_limit=21000;  %Default limit
        else
            center_dist_limit=run_options.center_dist_limit;
            %center_dist_limit=Inf;
        end
        disp(sprintf('Succesful localizations only, ranges  less than %i km from centerline used.',center_dist_limit/1000));
        Igood=[];
        error_ellipse=[];
        abratio=[];
        for KK=1:length(locations{I}),
            if ~isempty(findstr(lower(locations{I}{KK}.position.outcome),'successful'))
                
                %p = [localized{J}.utmx localized{J}.utmy];
                p=locations{I}{KK}.position.location;% Whale call locations
                %p1 = mean(d_coords(1:2,:));  % Get endpoints of a line segment through
                %p2 = mean(d_coords(6:7,:));  % Get endpoints of a line segment through
                %   the approximate center of the array.
                %center_distance = segment_points_dist_2d(p1,p2,p);
                
                if locations{I}{KK}.position.center_dist<=center_dist_limit
                    Igood=[Igood KK];
                    
                    [error_ellipse,abratio]=get_error_params(error_ellipse,abratio);
                    if run_options.Ikeep_only==1
                        tmp=locations_ctime{I}(KK,:);
                        locations_ctime{I}(KK,:)=0;
                        Ikeep=locations{I}{KK}.position.Ikeep;
                        if sum(tmp>0)>length(Ikeep)
                            %     keyboard;
                        end
                        locations_ctime{I}(KK,Ikeep)=tmp(Ikeep);
                    end
                    
                end
            end
        end
        fprintf('Out of %i automated surviving locations, %i successfully localized and have distance less than %i km \n', ...
            length(locations{I}),length(Igood),center_dist_limit/1000);
        
        locations{I}=locations{I}(Igood);
        locations_ctime{I}=locations_ctime{I}(Igood,:);
    else
        disp('Automated localizations do not have to be successful');
    end
    
    
end  %date looop

%%%%%%%%%%%Subfunctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [error_ellipse,abratio]=get_error_params(error_ellipse,abratio)
        if isfield(locations{I}{KK}.position,'Area')
            error_ellipse=[error_ellipse locations{I}{KK}.position.Area];
            abratio=[abratio locations{I}{KK}.position.major/locations{I}{KK}.position.minor];
        else
            error_ellipse=[error_ellipse -1];
            abratio=[abratio -1];
        end
    end

%%%%%%%%%%%%%function trim_max_ctime%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function trim_max_ctime
        %%Trim fields in auto object to be less than max_ctime...
        
        raw_names=fieldnames(auto.raw_station(1));
        for II=1:length(goodName{1})
            if ~isfield(auto.station(II),'ctime_min')
                continue
            end
            Itrim=find(auto.station(II).ctime_min<=max_ctime);
            auto.station(II)=trim_station(auto.station(II),Itrim);
            disp(sprintf('%i automated detections in station occur before %s',length(Itrim),ctime2str(max_ctime)));
            
            for JJ=1:length(raw_names),
                if ~isempty(findstr(raw_names{JJ},'detctions'))||isempty(auto.raw_station(II).raw_detections)||isempty(auto.raw_station(II).interval_detections)
                    continue
                end
                Itrim=find(auto.raw_station(II).(raw_names{JJ}).ctime<=max_ctime);
                auto.raw_station(II).(raw_names{JJ}).ctime=auto.raw_station(II).(raw_names{JJ}).ctime(Itrim);
                
                auto.raw_station(II).(raw_names{JJ}).duration=auto.raw_station(II).(raw_names{JJ}).duration(Itrim);
            end
            
            disp(sprintf('Maximum time in automated station %s is %s', goodName{I}{II},ctime2str(max(auto.station(II).ctime_min))));
        end
        
        
        if ~isfield(auto,'locations')
            locations=[];
            locations_ctime=[];
            return
        end
        
        Itrim=[];
        for KK=1:length(auto.locations)
            Ngood=(find(auto.locations{KK}.station_indicies>0));
            if mean(auto.locations{KK}.ctime_min(Ngood))<=max_ctime
                Itrim=[Itrim KK];
            end
        end
        fprintf('Out of %i original locations, %i occur before %s \n',length(auto.locations),length(Itrim),ctime2str(max_ctime));
        auto.locations=auto.locations(Itrim);
        auto.locations_ctime=auto.locations_ctime(Itrim,:);
        
        
    end
fprintf('Finished loading automated results\n\n');
end


function d = segment_points_dist_2d ( p1, p2, p )

%% SEGMENT_POINTS_DIST_2D: distance ( line segment, point ) in 2D.
%
%  Discussion:
%
%    A line segment is the finite portion of a line that lies between
%    two points.
%
%    The nearest point will satisfy the condition
%
%      PN = (1-T) * P1 + T * P2.
%
%    T will always be between 0 and 1.
%
%  Modified:
%
%    03 May 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real P1(2), P2(2), the endpoints of the line segment.
%
%    Input, real P(2), the point whose nearest neighbor on the line
%    segment is to be determined.
%
%    Output, real D, the distance from the point to the line segment.
%
%  Modified 19 February 2009 from Burkardt's SEGMENT_POINT_DIST_2D
%    to handle multiple points in an m-by-2 matrix (each row is a point).
%    Still assumes a single line segment.         Chris Nations

p1 = p1(:)';
p2 = p2(:)';
[m,n] = size(p);
if n~=2,
    error('p must be m-by-2')
end
one = ones(m,1);
p1 = one*p1;
p2 = one*p2;

%  If the line segment is actually a point, then the answer is easy.

if p1==p2,
    t = zeros(m,2);
else
    bot = sum((p2 - p1).^2,2);
    t = sum((p - p1).*(p2 - p1),2) ./ bot;
    t = max(t,0);
    t = min(t,1);
    t = t*[1 1];
end

pn = p1 + t.*(p2 - p1);
d = sqrt(sum((pn - p).^2,2));
end





