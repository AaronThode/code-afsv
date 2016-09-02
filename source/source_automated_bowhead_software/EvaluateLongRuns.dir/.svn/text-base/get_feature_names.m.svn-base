%%%Obtain feature names and secondary index if needed..
function [feature_name, feature_index]=get_feature_names(features);

prompt='Number of features:';
answer=inputdlg({'Number of features'},'Feature Name Selector',1,{'2'});
Nfeatures=str2num(answer{1});
names=fieldnames(features{1});
%Nfeatures=4;
for JJ=1:Nfeatures,
    Ichc=menu(sprintf('Select Feature %i',JJ),names);
    feature_name{JJ}=names{Ichc};
    test_value=getfield(features{1},feature_name{JJ});
    if sum(size(test_value))==2,
        feature_index(JJ)=1;
    else
        feature_index(JJ)=input(sprintf('Select index between 1 and %i: ',length(test_value)));
    end


end