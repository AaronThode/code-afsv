function boat_detections=merge_structure(boat_detections,new_struct,Maxbands,Npeaks)


if isempty(boat_detections)
    varnames=fieldnames(new_struct);
    
    for I=1:length(varnames)
        boat_detections.(varnames{I})=-1*ones(Maxbands,Npeaks,'single');
        
    end
    boat_detections.index=1;
    boat_detections.maxbands=1;
end

%Transfer new struct to boat_detections
varnames=fieldnames(new_struct);
II=boat_detections.index;
for I=1:length(varnames)
    tmp=new_struct.(varnames{I});
    if isempty(tmp)
        return  %No change to structure if no new information
    end
    if size(tmp,2)>1
        tmp=tmp';
    end
    boat_detections.(varnames{I})(1:length(tmp),II)=tmp;
    boat_detections.maxbands=max([length(tmp) boat_detections.maxbands]);

end
    

%%If extra input, then trim the results
if nargin>4
    for I=1:length(varnames)
        boat_detections.(varnames{I})=boat_detections.(varnames{I})(1:boat_detections.maxbands,1:boat_detections.index);
    end
end

boat_detections.index=II+1;


end