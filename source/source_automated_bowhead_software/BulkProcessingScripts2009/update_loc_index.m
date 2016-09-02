function [loc_index,Nlocs,problem_index]=update_loc_index(loc_index,Nlocs,match_index,dt,Ianchor,Iother)
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
                loc_index.(names{K})=[loc_index.(names{K}) zeros(Nrows,100)];
            end
        end
        loc_index.I(Ianchor,Nlocs)=match_index(1,I);
        loc_index.I(Iother,Nlocs)=match_index(2,I);
        
        continue;
    %if anchor detection used, append location
    elseif ~isempty(Imatch1)&isempty(Imatch2)
        if length(Imatch1)>1
            disp(sprintf('Index %i is repeated %i times in loc_index, row %i',match_index(1,I),length(Imatch1),Ianchor));
        end
        loc_index.I(Iother,Imatch1)=match_index(2,I)*ones(1,length(Imatch1));
    %if other detection used, append location
    elseif ~isempty(Imatch2)&isempty(Imatch1)
        if length(Imatch2)>1
            disp(sprintf('Index %i is repeated %i times in loc_index, row %i',match_index(2,I),length(Imatch2),Iother));
        end
        loc_index.I(Ianchor,Imatch2)=match_index(1,I)*ones(1,length(Imatch2));

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
        end
    end

end

end