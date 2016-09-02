function [Ipass,feature_matrix]=crude_filter_feature_locations(locations,filter_names,filter_criteria,min_stations,Idebug)

%filter criteria is is a [2 Ncrit] matrix
% filter_names is a cell array with Ncrit cells

%filter_names={'Contour_global_bandwidth','Contour_fmax'};
%Nstations=length(locations{1}.dt);
Ipass=1:length(locations);
feature_matrix=zeros(length(filter_names),length(locations));
for I=1:length(locations)
    Ikeep=find(locations{I}.station_indicies>0);  
    Igood=Ikeep;
    for J=1:length(filter_names) 
        values=locations{I}.(filter_names{J})(Igood);   
        feature_matrix(J,I)=max(locations{I}.(filter_names{J}));
        
        Ipass1=find(values>=filter_criteria(1,J)&values<=filter_criteria(2,J));
        if isempty(Ipass1),
            Igood=[];
        else
           Igood=Igood(Ipass1);
        end
       
       
    end
    
    if length(Igood)<min_stations
        Ipass(I)=-1;
    end

end

Ipass=Ipass(Ipass>0);
end