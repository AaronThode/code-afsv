function [features_out,ctimes_out,Istrip,Ifail]=filter_feature_vector(feature_vector,feature_name,ctimes,optvec,filt_str,operator,debug),
%function [features_out,ctimes_out,Istrip,Ifail]=filter_feature_vector(feature_vector,feature_name,ctimes,optvec,filt_str,operator,debug),

% Input:
%   feature_vector:  feature_vector{I}(J,K) contains value of feature I,
%       contained in detection J, shape K.
%   feature_name:  cell array of strings (feature names)
%   ctimes:  c-times associated with each detection J.
%   optvec:  2-D vector of acceptable ranges for each feature.
%       [feature_1_min feature_1_max; feature_2_min feature_2_max ...]
%   filt_str: 'all' or 'first' --how many shapes per detection to evaluate?  If
%           'all' accept detection if at least one shape passes all criteria..
%   operator:  if 'and', shape has to pass all rules.  If 'or' needs to
%   pass only one
% Output:
%   features_out: filtered feature_vector
%   ctimes_out: filtered ctimes
%   Istrip:  Indicies of feature_vetor that pass filtering

if isempty(optvec)
    optvec=manual_select_thresholds(feature_name);
    optvec=reshape(optvec,2,[])';
end


total_features=size(feature_vector{1},1);
Nfeatures=length(feature_vector);
Nshapes=size(feature_vector{1},2);
if strcmp(filt_str,'first'),
    Nshapes=1;
end

if length(feature_name)~=Nfeatures,
    disp('feature_name does not list same number of features as feature_vector');
    keyboard;
    return;
end
disp(sprintf('At start of discrimination %i detections',total_features));

%%Plot first two features against each other.
if debug>0,
    try,
        figure(1);
        plot(feature_vector{1}(:,1),feature_vector{2}(:,1),'x');
        %keyboard;
    end

end

Istrip=1:total_features;
Ipass=[];

for I=1:Nfeatures,
    Ifail{I}=[];

    for J=1:Nshapes,
        Itest=(find(feature_vector{I}(Istrip,J)>=optvec(I,1)&feature_vector{I}(Istrip,J)<=optvec(I,2)));
        Ipass=[Ipass Istrip(Itest)];

        if debug>0,
            Itest=(find(feature_vector{I}(Istrip,J)<optvec(I,1)|feature_vector{I}(Istrip,J)>optvec(I,2)));
            Ifail{I}=[Ifail{I} Istrip(Itest)];
        end
    end

    %Remove duplicates
    %if Nshapes>1,
    %Ipass=sort(Ipass,2,'ascend');
    %Idiff=[diff(Ipass) 1];
    %Ipass=Ipass(Idiff>0);
    Ipass=unique(Ipass);

    if debug>0&~isempty(Ifail{I}),
        Ifail{I}=unique(Ifail{I});
        %Ifail{I}=sort(Ifail{I},2,'ascend');
        %Idiff=[diff(Ifail{I}) 1];
        %Ifail{I}=Ifail{I}(Idiff>0);
    end
    %end

    if strcmp(operator,'and'),
        Istrip=Ipass;
        Ipass=[];
        disp(sprintf('After feature %s discrimination, %i detections pass all rules',feature_name{I},length(Istrip)));

    else
        %Ipass=sort(Ipass,2,'ascend');
        %Idiff=[diff(Ipass) 1];
        %Ipass=Ipass(Idiff>0);
        
        disp(sprintf('After feature %s discrimination, %i detections pass at least one rule',feature_name{I},length(Ipass)));

    end
    % if exist('debug'),
    %   figure(gcf+1)
    %end
end
if strcmp(operator,'or'),
    Istrip=Ipass;
end

ctimes_out=ctimes(Istrip);
for I=1:Nfeatures,
    features_out{I}=feature_vector{I}(Istrip,:);
end



