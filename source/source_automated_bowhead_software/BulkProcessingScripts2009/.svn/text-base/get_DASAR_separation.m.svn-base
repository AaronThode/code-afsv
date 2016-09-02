%%%%%%%%get_DASAR_separation(Isite,name1,name2);
%function distance=get_DASAR_separation(Isite,name1,name2);
%  Isite: integer representing site number
%  name1: string with single letter. Can be a complete filename with '/'.
%       In that case, file extension '.gsi' is assumed.
%  name2; string with single letter...
% Output:
%   distance: separation in kilometers
function distance=get_DASAR_separation(Isite,name1,name2)

%%If Isite==0, actually load Site 6 (this is Shell Triplet)


[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);
miss_flag=1;
for Iyear=7:20
    
    if ~isempty(findstr(name1,sprintf('S%i%02i',Isite,Iyear)))
        load(sprintf('../DASARlocations/DASAR_locations_20%02i.mat',Iyear));
        miss_flag=0;
        
    end
    
end

if Isite==0
    Isite=6;
end

if miss_flag==1
    error(sprintf('Missing location data for %s',name1));
end

%Trim input names if they are complete file paths

name{1}=lower(name1(2:end));
name{2}=lower(name2(2:end));

for I=1:2,
    Islash=max(findstr(name{I},'/'))+1;
    if ~isempty(Islash)
        name{I}=name{I}(Islash:(end-4));
    end
end

for I=1:length(Site{Isite}.stations)
    if ~isempty(strfind(lower(Site{Isite}.stations{I}),'t'))  %Avoid the 'T' in the file name
        %disp('skip t')
        continue
    end
    
    %Use name(2:end) to skip the 's'
    if ~isempty(strfind(name{1}(2:end),lower(Site{Isite}.stations{I})))
        loc1=[Site{Isite}.northing(I) Site{Isite}.easting(I)];
    elseif ~isempty(strfind(name{2}(2:end),lower(Site{Isite}.stations{I})))
        loc2=[Site{Isite}.northing(I) Site{Isite}.easting(I)];
        
    end
    
end

distance=sqrt(sum((loc2-loc1).^2));
