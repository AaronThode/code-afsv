function Event=update_GSI_localization(Event,Istation, bearing,kappa)
%Recalculate and update fields for GSI localization, when a bearing at a
%   station is altered, added, or removed.
%   When is station is removed, used 'NaN' for bearing
% Istation, bearing, and kappa must be numbers, not strings

%First check that station number is consistent
%if Istation~=str2num(Event.Istation)
%   errormsg('update_GSI_localization: mismatched Istation');
%   return
%end

%Update bearing information if needed


if exist('Istation','var')
    if Istation==str2num(Event.Istation)
        Event.bearing=num2str(bearing);
    end
    
    Event.localization.bearings_all(Istation)=bearing;
    Event.localization.kappa(Istation)=kappa;
    
    %If bearing is a NaN, then kill the hashtag too
    if isnan(bearing)
        nblanks=size(Event.link_hashtags,2);
        Event.link_hashtags(Istation,:)=blanks(nblanks);
        Event.link_hashtags(Istation,(end-1):end)='-1';
    end
    
    %If hashtags are NaN, then kill bearings
    hashtags=str2num(Event.link_hashtags);
    Inan=find(hashtags<0);
    Event.localization.bearings_all(Inan)=NaN;
    Event.localization,kappa(Inan)=NaN;
    
else
    Istation=str2num(Event.Istation);  %Use this to recompute range using my own information.
end

%Recompute localization
theta=Event.localization.bearings_all;
kappa=Event.localization.kappa;
Ikeep=find(~isnan(theta));

station_position=Event.localization.station_position;
DASAR_coords=[station_position.easting station_position.northing];

[Event.localization.location,Event.localization.Qhat,~,Event.localization.outcome] = vmmle_r(theta(Ikeep),DASAR_coords(Ikeep,:),'h',kappa(Ikeep));
%if ~strcmpi(Event.localization.outcome,'Successful')
%    keyboard
%end
    

Event.localization.Nused=length(Ikeep);

%Update range and bearing
Event=update_localization_fields(Event);
%Event.localization.range=sqrt((station_position.easting(Istation)-Event.localization.location(1)).^2+ ...
%    (station_position.northing(Istation)-Event.localization.location(2)).^2);
%Event.range=num2str(Event.localization.range/1000);


end