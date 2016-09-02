function [Ipass,feature_matrix]=crude_filter_feature_stations(station,filter_names,filter_criteria)

%filter criteria is is a [2 Ncrit] matrix
% filter_names is a cell array with Ncrit cells

%filter_names={'Contour_global_bandwidth','Contour_fmax'};
%Nstations=length(locations{1}.dt);


Ipass=1:length(station.ctime_min);
feature_matrix=zeros(length(filter_names),length(Ipass));
for J=1:length(filter_names)
    feature_matrix(J,:)=station.(filter_names{J});

    values=station.(filter_names{J})(Ipass);

    Ipass1=find(values>=filter_criteria(1,J)&values<=filter_criteria(2,J));
    if isempty(Ipass1),
        Ipass=[];
        
    else
        Ipass=Ipass(Ipass1);
    end


end

