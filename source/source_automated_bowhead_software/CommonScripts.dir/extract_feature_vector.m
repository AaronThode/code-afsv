%function feature_vector=extract_feature_vector(features,feature_name,feature_index,Nshapes),
% features--features{I}(J) is a structure array of features of shape J from detection I.
% feature_name:  cell array of strings (feature names)
% feature_index:  vector of integers listing which index of multi-element
% features to extract (e.g.Centroid is a 2-element vector)
% Nshapes:  maximum number of shapes to extract per detections.
%
% Output:
%   feature_vector:  feature_vector{I}(J,K) contains value of feature I,
%       contained in detection J, shape K.
function feature_vector=extract_feature_vector(features,feature_name,feature_index,Nshapes),

if length(features)==0,
    feature_vector=[];
    return
end
dT=features{1}.dT;
dF=features{1}.dF;

%keyboard;

for JJ=1:length(feature_name),
    feature_vector{JJ}=NaN*ones(length(features),Nshapes);
    for Ifeature=1:length(features) ,
        list=features{Ifeature};
        for Ishape=1:min([Nshapes length(list)]),
            %value=transform(getfield(features(Ifeature),feature_name{JJ}),feature_name{JJ});
            value=transform(list(Ishape).(feature_name{JJ}),feature_name{JJ},dT,dF);
            if Ishape>1&&(strcmp(feature_name{JJ},'Area2')|strcmp(feature_name{JJ},'vertical_percent')),  %This feature is only assigned to the largest detection
                value=feature_vector{JJ}(Ifeature,1);
                feature_index(JJ)=1;
            end
            if ~isempty(value),
                feature_vector{JJ}(Ifeature,Ishape)=value(feature_index(JJ));
            end
        end

    end

end


%             p.transform(JJ)=1;
%             if strcmp(feature_name{JJ},'Eccentricity')
%                 p.transform(JJ)=2;
%                 %feature_name{JJ}='-Log Eccentricity';
%             end
%
%             p.false{JJ}=NaN*ones(size(I_false_alarms ));
%             p.match{JJ}=NaN*ones(size(I_auto_match ));
%             p.missed{JJ}=NaN*ones(size(I_missed_calls ));

end

function y=transform(x,name,dT,dF);

y=[];
if isempty(x), return;end
switch name,
    case 'Eccentricity'
        %y=min([-log10(x) 0.5]);
        y=min([1-x 0.2]);
        %y=x;
        %disp('taking the negative log of eccentricity')
    case {'Area','Area2'}
        y=dT*dF*x;
        %disp('Area')
    case 'Centroid'
        %y(1)=x(1)*dT;
        %y(2)=x(2)*dF;
        y=x;
    otherwise,
        y=x;

end

end