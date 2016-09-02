%%%%%%%load_DASAR_coords.m
% function DASAR_coords=load_DASAR_coords(Icase,Isite,goodFile,wantLL)
%  if wantLL exists, return lat long
%   Icase: string that must contain either 'Shell07' or 'Shell08' or
%       'Arctic_2008'
%   Isite: integer of site in question
%   goodFile: file locations, used to determine what DASARs to keep or throw away...
function DASAR_coords=load_DASAR_coords(Icase,Isite,goodFile,wantLL)
%persistent Site
%if isempty(Site)

%%If Isite==0, actually load Site 6 (this is Shell Triplet)



[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);
miss_flag=1;
for Iyear=7:20
    if  (~isempty(findstr(Icase,sprintf('Shell%02i',Iyear))))||~isempty(findstr(Icase,sprintf('Arctic_20%02i',Iyear)))
        try
            load(sprintf('../DASARlocations/DASAR_locations_20%02i.mat',Iyear));
        catch
            load(sprintf('DASAR_locations_20%02i.mat',Iyear));
            
        end
        miss_flag=0;
    end
    
end

%%If site desired is zero, actually use index 6
if Isite==0
    Isite=6;
end

if miss_flag==1
    error(sprintf('Missing location data for %s',Icase));
end

if ~exist('wantLL')
    DASAR_coords=[Site{Isite}.easting Site{Isite}.northing];
elseif exist('wantLL')&isempty(goodFile);
    DASAR_coords=[Site{Isite}.lon' Site{Isite}.lat'];
    
end

%Strip away bad DASARs...Dasars that we don't use
if exist('goodFile')&&~isempty(goodFile)
    goodletter=[];
    for I=1:length(goodFile),
        Islash=max(findstr(goodFile{I},'/'))+5;
        if isempty(Islash) %eg, entered goodName
            Islash=5;
        end
        goodletter=[goodletter goodFile{I}(Islash)];
    end
    
    Ikeep=[];
    for I=1:length(Site{Isite}.stations)
        if ~isempty(strfind(lower(Site{Isite}.stations{I}),'t'))  %Avoid the 'T' in the file name
            continue
        end
        if ~isempty(strfind(goodletter,upper(Site{Isite}.stations{I})))
            Ikeep=[Ikeep I];
        end
        
        
    end
    
    DASAR_coords=DASAR_coords(Ikeep,:);
end

