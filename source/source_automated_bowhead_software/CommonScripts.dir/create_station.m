%function station=create_station(equalization_freq,equalization,features,debug_ctimes,debug_durations,feature_name,global_feature_name,feature_index,Nsegments),
% Extract and sort feature detections, merging 'best_call' and 'shape'
% indicies..
% Input:
% equalization_freq: vector of frequencies (Hz) that will define equalization
%   curve
% equalization: [Nfreq x Ncall] matrix of equalizations in dB PSD.
% features--features{I}(J) is a structure array of features of detection J from
%   energy detection I.  Output of 'process_one_unit': best_calls.features
% final_image: cell array of spare matricies containing labeled images.
% debug_ctime:  ctime of beginning of time series used for morph
%       processing: same as 'cstart' in extract_image_features.
% debug_durations: duration of call used to reproduce extract_image_features
%   output.
% feature_name:  cell array of strings (feature names)
% global_feature_name: cell array of strings that reference features that
%   have a single value per call (i.e. Contour and Total features)
% feature_index: [1 x Nfeature] vector of integers listing which
%               index of a multi-element
%               feature to extract (e.g. "2" to access frequency  of Centroid, which is a 2-element vector)
% Nsegments:  maximum number of segments permitted per call.

% Output...
% station:  a structure with following fields, most of which have 'Ncall'
%   columns, where a call is a collection of segments that belong to a single
%   vocalization event...
%           .Image; {1x Ncall} cell array of 2-D image, with east segment I
%               labeled in image with integer 'I'
%           .Nsegments:  number of segments in call (max value Nsegments);
%           .SNR:  [1 x Ncall] vector of dB signal to noise ratios for entire call (all segments included)
%           .SEL: [1 x Ncall] of total SEL (sum of .feature.SEL)
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
%               The segments are arranged by maximum time-bandwidth product
%                 


function station=create_station(equalization_freq,equalization,features,final_image,debug_ctimes,debug_durations,feature_name, ...
    global_feature_name,feature_index,Nsegments)

station.param.equalization_freq=equalization_freq;

%Initialize...

station.feature=[];
station.indicies=[];
station.feature.ctime=[];
station.ctime_debug=[];
station.duration_debug=[];
station.SNR=[];
station.SEL=[];
station.param.feature_name=feature_name;
station.param.feature_index=feature_index;
station.equalization=[];
station.ctime_min=[];
station.Nsegments=[];
for JJ=1:length(global_feature_name),
    station.(global_feature_name{JJ})=[];
end

if length(features)==0,
    return
end
dT=features{1}.dT;
dF=features{1}.dF;

Ncalls_raw=length(features);
Ncalls=Ncalls_raw;
for J=1:length(features)
    Ncalls=Ncalls+length(features{J})-1;
end
Icount=0;

%Preallocate memory
for JJ=1:length(feature_name),
    if strcmp(feature_name{JJ},'Image')&Nsegments==1,
        station.feature.Image{1}=0;
    else
        station.feature.(feature_name{JJ})=-Inf*ones(Nsegments,Ncalls);
    end

end

station=orderfields(station);
station.feature=orderfields(station.feature);

station.feature.ctime=-Inf*ones(Nsegments,Ncalls);  %ctime of each segment of call
station.ctime_debug=-Inf*ones(1,Ncalls);  %ctime of raw detection, used to derive morph calculations, permits easy display of call
station.ctime_min=station.ctime_debug;
station.duration_debug=station.ctime_debug;
station.SNR=station.ctime_debug;
station.SEL=station.ctime_debug;
station.indicies=-Inf*ones(2,Ncalls);
station.equalization=zeros(length(equalization_freq),Ncalls);

        
for I=1:Ncalls_raw %For each energy detection
    if rem(I,1000)==0,disp(sprintf('create_station: %i percent done',round(100*I/Ncalls_raw)));end

     for J=1:length(features{I}),  %Review each detection within energy detection
        Nshape_max=min([Nsegments length(features{I}(J).ctime)]);
        Icount=Icount+1;

        %Items that do not differ over J loop
        %tmp=[features{I}(J).ctime];
        station.Nsegments(Icount)=Nshape_max;
        station.feature.ctime(1:Nshape_max,Icount)=features{I}(J).ctime(1:Nshape_max)';
        
         
        station.ctime_min(Icount)=min(features{I}(J).ctime(1:Nshape_max));
        station.ctime_debug(Icount)=debug_ctimes(I);
        station.duration_debug(Icount)=debug_durations(I);
        station.SNR(Icount)=features{I}(J).SNR;
        station.SEL(Icount)=(features{I}(J).SEL);
        station.Image{Icount}=final_image{I}{J};
        
        station.indicies(:,Icount)=[I; J];
        station.equalization(:,Icount)=equalization{I};

        %%Features that only have one value for the entire detection-no
        %%segments.
        for JJ=1:length(global_feature_name),
           station.(global_feature_name{JJ})(Icount)=features{I}(J).(global_feature_name{JJ})(1); 
        end
       
        %%New addition to add image...%%
        for Ishape=1:Nshape_max,  %Number of 

            for JJ=1:length(feature_name),
                
                %%Go through remaining features, making sure no NaN's exist
                value_raw=features{I}(J).(feature_name{JJ});
                if size(value_raw,2)>1,
                    value=transform(value_raw(feature_index(JJ),Ishape),feature_name{JJ},dT,dF);
                    station.feature.(feature_name{JJ})(Ishape,Icount)=value;  %Feature index used for features that are vectors, e.g. Bounding Box
                else  %A value that is the same for all components of a call (e.g., SNR or SEL)
                    station.feature.(feature_name{JJ})(1,Icount)=transform(value_raw(feature_index(JJ)),feature_name{JJ},dT,dF);

                end
            end
        end
    end
end

disp(sprintf('create_station: Final fraction calls processed: %i percent ',round(100*I/Ncalls_raw)));
% if strcmp(feature_name{JJ},'Image')&Nsegments==1,
%     return;
% end

%%Clean out detections that don't have results for features...
Igood=1:length(station.SNR);
for JJ=1:length(feature_name),
    if Nsegments>1
        Igood2=find(~any(isnan(station.feature.(feature_name{JJ})(:,Igood))));
    else
        Igood2=find(~(isnan(station.feature.(feature_name{JJ})(:,Igood))));

    end
    Igood=Igood(Igood2);
end

 for JJ=1:length(global_feature_name),
     if Nsegments>1
        Igood2=find(~any(isnan(station.(global_feature_name{JJ})(:,Igood))));
    else
        Igood2=find(~(isnan(station.(global_feature_name{JJ})(:,Igood))));

    end
    Igood=Igood(Igood2);
 end
         
station=trim_station(station,Igood);

%Finally, sort by time...
[junk,Isort]=sort(station.ctime_min);
station=trim_station(station,Isort);


% for JJ=1:length(names),
%     if strcmp(names{JJ},'ctime'),
%         continue;
%     end
%     station.feature.(names{JJ})=station.feature.(names{JJ})(:,Isort);
% 
% end
% %station.ctime=station.ctime(Isort);
% station.ctime_debug=station.ctime_debug(Isort);
% station.duration_debug=station.duration_debug(Isort);
% station.SNR=station.SNR(Isort);
% station.indicies=station.indicies(:,Isort);
% station.equalization=station.equalization(:,Isort);


end


function y=transform(x,name,dT,dF)

y=[];
if isempty(x), return;end
switch name,
    case 'Eccentricity'
        %y=min([-log10(x) 0.5]);
        %y=min([1-x 0.2]);
        %y=x;
        y=x;
        %disp('taking the negative log of eccentricity')
    %case {'Area','Area2'}
     %   y=dT*dF*x;
        %disp('Area')
    case 'Centroid'
        %y(1)=x(1)*dT;
        %y(2)=x(2)*dF;
        y=x;
    otherwise,
        y=x;

end

end
