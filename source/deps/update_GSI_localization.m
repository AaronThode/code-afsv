function Event=update_GSI_localization(Event,Istation, bearing,kappa)
%Recalculate and update fields for GSI localization, when a bearing at a
%   station is altered, added, or removed.
%   When is station is removed, used 'NaN' for bearing

%First check that station number is consistent
%if Istation~=str2num(Event.Istation)
%   errormsg('update_GSI_localization: mismatched Istation');
%   return
%end

%Next, recompute position
station_position=Event.localization.station_position;
DASAR_coords=[station_position.easting station_position.northing];

if Istation==str2num(Event.Istation)
    Event.bearing=num2str(bearing);
end
Event.localization.bearings_all(Istation)=bearing;
%Event.localization.kappa(Istation)=kappa;

theta=Event.localization.bearings_all;
kappa=Event.localization.kappa;
Ikeep=find(~isnan(theta));

[Event.localization.location,Event.localization.Qhat,~,Event.localization.outcome] = vmmle_r(theta(Ikeep),DASAR_coords(Ikeep,:),'h',kappa(Ikeep));

%Update rest of fields
Event.position=num2str(Event.localization.location);
Event.localization.range=sqrt((station_position.easting(Istation)-Event.localization.location(1)).^2+ ...
    (station_position.northing(Istation)-Event.localization.location(2)).^2);
Event.range=num2str(Event.localization.range/1000);
Event.localization.Nused=length(Ikeep);
end