%%%function [locations, locations_ctime]=cross_channel_match(Isite,station,feature_params,goodFile,run_options,debug_params),
%%% Given a collection of detections across stations, match features across
%%% stations to link calls
%%% Aaron Thode, Sept 18, 2008
%%% Revised August 4. 2009
%
% Input:
%
%   Isite:  Integer representing a site number--used to acquire distances
%       between stations.
%
%   station:  a structure with following fields, most of which have 'Ncall'
%           columns, where a call is a collection of segments that belong to a single
%               vocalization event...
%           .SNR:  [1 x Ncall] vector of dB signal to noise ratios for entire call (all segments included)
%           .Totalduration [1 x Ncall]: Total duration of bounding box
%               surrounding all call segments, in seconds.
%           .Totalfmax [1 x Ncall]: Maximum frequency (Hz) of bounding box
%               surrounding all call segments.
%           .Totalfmin [1 x Ncall]: same as above, but minimum frequency in
%               Hz.
%           .ctime: [Nsegments xNcall] matrix of c-times for each call
%           segment
%           .ctime_debug: [1 xNcall] vector of c-times that reconstruct
%               original time series subject to morphological processing.
%               Intended for display purposes.
%           .ctime_min: [1 xNcall] vector of c-times that contain the start
%               time of the earliest call segment.  Intended for
%               matched-filtering.
%           .duration_debug: [1 x Ncall] vector of durations (sec) used
%               with .ctime_debug to reconstruct original analyzed time series.
%           .equalization: [Nfreq x Ncall] matrix of power spectral
%               densities (dB re 1uPa^2/Hz), associated with each final
%               call
%           .param.equalization_freq: [1 x Ncall] vector of frequencies (Hz)
%               associated with equalization.  Identical to epononymous input.
%           .feature:  a structure of many features.  See below for
%               details.
%           .param.feature_index: [1 x Nfeature] vector of integers listing which
%               index of a multi-element
%               feature to extract (e.g. "2" to access frequency  of Centroid, which is a 2-element vector)
%           .param.feature_name: {1 x Nfeature} cell matrix describing the
%               desired subset of features extracted from 'features' into
%               'station'.
%           .indicies: [2x Ncall] list of reference indicies used to link
%               output to original 'features' structure.  Row 1,2 access
%               features{1}(2)--that is, row 1 access cell, row 2 access
%               segment in cell
%
%
%           MORE DETAILS ON .feature structure..
%               .each string in 'feature_name' above will have a [Nsegments
%               x Ncall] matrix containing data on each segment (row) in each call (column).
%               The segments are arranged by minimum frequency (row 1 has
%               lowest frequency segment)

%   feature_params:
%       .weight=[1 0 0];
% 	    .tol=[25 0 0 ];
% 	    .names: %Same as station(I).feature_name;
%       .Fs: sampling rate in Hz
%
%   goodFile:   cell array of files names associated with station, e.g. goodFile{1}='D07s4aT20070829T000000'
%   run_options: structure with field
%       .min_stations: mininum number of DASARS required to constitute a
%           'location'
%   debug_params:
%       J_anchor:  Station index to begin for anchor station search, usually 1
%           unless debugging.
%       J_anchor_start:  index of call in anchor station to be checked
%       cross_channel:  if 1, some debug printing, if 2, more debug
%               printing..
%
%  Output:
%
%  locations:  cell array of location objects, organized as
%             locations{J}, which refers to location number J
%
%             All fields are 7 element vectors, first element
%               referring to station 'a', etc.
%
%             dt: [7x1 double]  arrival times relative to anchor station.
%               If station does not have location, value is NaN
%             indicies: [7x1 double]  indicies of original station objects.
%                e.g. locations{I}{J}.indicies(K) is the index of station
%                (K).  A station object can then be used to refer to
%                original 'best_calls' variable.  value of -1 is station
%                does not have.
%             ctime: [7x1 double].  Value 0 if station does not have.
%             SNR: [7x1 double]
%             mismatch: [0 Inf Inf Inf 0.0608 0.0034 0.0080]  %Mismatch of
%               best-fit.  Value of Inf if station does not have.
%             feature(K).(feature_name), this field uses name contained in
%                   station(J).feature_name, and vector has feature values
%               across all arrays.  feature(K).(feature_name) acesses the
%               value of 'feature_name' at station K.  If station K doesn't
%               have location then feature values are empty matricies.
%
%   locations_ctime: a [ Ncall x Nstation] matrix of ctimes associated with
%               locations.  Useful for sorting.
%   loc_index: convenient matrix of linked station indicies...


function [locations,locations_ctime,loc_index]=cross_channel_match(Isite,station,feature_params,goodFile,run_options,debug_params)

locations_ctime=[];
J_anchor=1;
Nanchors=length(goodFile);

if exist('debug_params')&&debug_params.cross_channel>1,
    J_anchor=debug_params.J_anchor;
    Nanchors=J_anchor;

end


feature_params.SNR_tol=0;
feature_params.min_stations=run_options.min_stations;  %Minimum number of matches between stations...
feature_names=feature_params.names;


locations_all=[];
loc_index.I=zeros(Nanchors,500);
loc_index.dt=NaN*ones(Nanchors,500);
loc_index.dt_xcorr=NaN*ones(Nanchors,500);
loc_index.credibility_score=zeros(1,500);
%loc_index.mismatch=loc_index.dt;
Nlocs=0;
%Itotal=0;

for Ianchor=J_anchor:Nanchors,  %Start cycling through DASARS
    if isempty(station(Ianchor))||isempty(station(Ianchor).feature.ctime),
        continue;
    end

    station_anchor=station(Ianchor);
    if exist('debug_params')&&debug_params.cross_channel>1,
        station_anchor=trim_station(station_anchor,debug_params.J_anchor_start);  %Only keep desired indicies..
    end

    %%Initialize a location object for this anchor station
    %locations=initialize_location_object(station_anchor,Ianchor,length(goodFile),'basic');

    for Iother=Ianchor:length(goodFile)
        %keyboard;
        if Iother==Ianchor||isempty(station(Iother))||isempty(station(Iother).feature.ctime)||isempty(station_anchor.feature.ctime),
            continue,
        end

        %if exist('debug_params')&&debug_params.cross_channel>1,

        %    debug_params.J_other=Iother;
        %end
        distance=get_DASAR_separation(Isite,goodFile{Ianchor},goodFile{Iother});
        feature_params.time_tol=distance/1400+run_options.dt_slop;

        disp('%%%%%%cross_channel_match.m%%%%%%%%%')
        disp(sprintf('Processing %i calls in %s, compare with %i calls in %s, time tolerance %6.2f',length(station_anchor.ctime_min),goodFile{Ianchor}, ...
            length(station(Iother).ctime_min),goodFile{Iother},feature_params.time_tol));

        tic
        [match_index,dt,dt_xcorr,best_mismatch]=match_stations(station_anchor,station(Iother),feature_params,debug_params,goodFile{Ianchor},goodFile{Iother});
        toc
    
        
        [loc_index,Nlocs,problem_index{Ianchor,Iother}]=update_loc_index(loc_index,Nlocs,match_index,dt,dt_xcorr,Ianchor,Iother);
    end  %Iother cycle
end

Igood=find(sum(loc_index.I)>0);
loc_index.I=loc_index.I(:,Igood);
loc_index.dt=loc_index.dt(:,Igood);
loc_index.credibility_score=loc_index.credibility_score(Igood);

Nobjects=length(Igood);
%options='basics';
locations=cell(1,Nobjects);

for I=1:size(loc_index.I,2)
    for J=1:Nanchors
        locations{I}.dt(J)=loc_index.dt(J,I);
        locations{I}.dt_xcorr(J)=loc_index.dt_xcorr(J,I);

        locations{I}.station_indicies(J)=loc_index.I(J,I);
        %locations{I}.mismatch(J)=loc_index.mismatch(J,I);
        locations{I}.credibility(I)=loc_index.credibility_score(I);
    end
end

%%%Limit locations to situations where call is present on minimum
%%%number of DASARS.
Ifinal=[];
if ~isempty(locations),
    for J=1:length(locations),
        Nstations=find(locations{J}.station_indicies>0);
        if length(Nstations)>=feature_params.min_stations
            Ifinal=[Ifinal J];

        end
    end
end
locations=locations(Ifinal);

%%Finally, copy all information from station to location for convenience..
station_features=fieldnames(station(1));
feature_names=fieldnames(station(1).feature);
if ~isempty(locations),
    for J=1:length(locations),
        for Istation=1:length(goodFile),
            Iwwant=locations{J}.station_indicies(Istation);
            %%Safety check,
            %if Istation==Ianchor,
            %    disp(sprintf('Stored SNR: %10.4f, new SNR: %10.4f',locations{J}.SNR(Istation),station(Istation).SNR(JJ)));
            %end

            for KK=1:length(station_features),
                switch station_features{KK},
                    case {'param','indicies','Image'}
                        continue
                    case 'feature'
                        for Iname=1:length(feature_names),
                            locations{J}.feature(Istation).(feature_names{Iname})=NaN;
                            if Iwwant>0,
                                locations{J}.feature(Istation).(feature_names{Iname})=station(Istation).feature.(feature_names{Iname})(:,Iwwant)';
                            end
                        end
                    otherwise
                        if Iwwant>0,
                            value=station(Istation).(station_features{KK})(:,Iwwant);
                            % if length(value)>1&Istation==2,
                            %     keyboard;
                            % end
                            locations{J}.(station_features{KK})(Istation,1:length(value))=value';
                        else
                            locations{J}.(station_features{KK})(Istation,:)=0;

                        end
                end
            end
        end
    end

    %locations_all=[locations_all locations];

    %Finally, rearrange locations by mean time, instead of anchor DASAR.
    for I=1:length(locations),
        ctime=locations{I}.ctime_min;
        Igood=find(ctime>0);
        mean_ctime(I)=mean(ctime(Igood));
        locations_ctime(I,:)=ctime';

    end

    [mean_ctime,Isort]=sort(mean_ctime);
    locations=locations(Isort);
    locations_ctime=locations_ctime(Isort,:);
end  %if isempty(locations)
end

function [loc_index,Nlocs,problem_index]=update_loc_index(loc_index,Nlocs,match_index,dt,dt_xcorr,Ianchor,Iother)
%Repeated indicies are allowed in loc_index, to provide for opportunity for
% discrepency resolution in the future

problem_index=[];
if isempty(match_index)
    return
end

for I=1:size(match_index,2)
    Imatch1=find(match_index(1,I)==loc_index.I(Ianchor,:));
    Imatch2=find(match_index(2,I)==loc_index.I(Iother,:));

    %If neither index present in loc_index.I, we have a new linkage
    if isempty(Imatch1)&isempty(Imatch2)
        Nlocs=Nlocs+1;
        if Nlocs>size(loc_index.I,2)
            names=fieldnames(loc_index);
            for K=1:length(names)
                Nrows=size(loc_index.(names{K}),1);
                if ~strcmp(names{K},'I')
                    loc_index.(names{K})=[loc_index.(names{K}) NaN*ones(Nrows,100)];
                else
                    loc_index.(names{K})=[loc_index.(names{K}) zeros(Nrows,100)];
                end
            end
        end
        loc_index.I(Ianchor,Nlocs)=match_index(1,I);
        loc_index.I(Iother,Nlocs)=match_index(2,I);
        loc_index.dt(Ianchor,Nlocs)=0;
        loc_index.dt_xcorr(Ianchor,Nlocs)=0;
        loc_index.dt(Iother,Nlocs)=dt(I);
        loc_index.dt_xcorr(Iother,Nlocs)=dt_xcorr(I);
        %loc_index.mismatch(Ianchor,Nlocs)
        continue;
        %if anchor detection used, append location
    elseif ~isempty(Imatch1)&isempty(Imatch2)
        if length(Imatch1)>1
            disp(sprintf('Index %i is repeated %i times in loc_index, row %i',match_index(1,I),length(Imatch1),Ianchor));
        end
        loc_index.I(Iother,Imatch1)=match_index(2,I)*ones(1,length(Imatch1));
        loc_index.dt(Iother,Imatch1)=loc_index.dt(Ianchor,Imatch1)+dt(I);  %t0-tanchor + tanchor-tother=t0-tother
        loc_index.dt_xcorr(Iother,Imatch1)=loc_index.dt_xcorr(Ianchor,Imatch1)+dt_xcorr(I);  %t0-tanchor + tanchor-tother=t0-tother

        %if other detection used, append location with 'Ianchor'
    elseif ~isempty(Imatch2)&isempty(Imatch1)
        if length(Imatch2)>1
            disp(sprintf('Index %i is repeated %i times in loc_index, row %i',match_index(2,I),length(Imatch2),Iother));
        end
        loc_index.I(Ianchor,Imatch2)=match_index(1,I)*ones(1,length(Imatch2));
        loc_index.dt(Ianchor,Imatch2)=loc_index.dt(Iother,Imatch2)-dt(I);  %t0-tother + (tother-tanchor)   =t0-tanchor
        loc_index.dt_xcorr(Ianchor,Imatch2)=loc_index.dt_xcorr(Iother,Imatch2)-dt_xcorr(I);  %t0-tother + (tother-tanchor)   =t0-tanchor

        %%What if both indicies have been used before?  Due to our
        %%restrictions the indicies should be only visible once
    else
        %%Good news situation: both indicies are in same location.
        %%Increases our credibility
        if Imatch1==Imatch2
            loc_index.credibility_score(Imatch1)=loc_index.credibility_score(Imatch1)+1;
        else
            %%Now a problem has arisen--these indicies have been assigned
            %%to separate localizations.  For now, just flag the issue
            problem_index=[problem_index I];
            %loc_index.I(Ianchor,Imatch2)=match_index(1,I)*ones(1,length(Imatch2));


        end
    end

end

end

