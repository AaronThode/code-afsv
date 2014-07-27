function Event=update_localization_fields(Event)
% Update the fields 'bearing','range', and 'position' for an Event, using
%    data stored in localization field
%end

%Update bearing information if needed

station_position=Event.localization.station_position;
DASAR_coords=[station_position.easting station_position.northing];
Istation=str2num(Event.Istation);

Event.localization.range=sqrt((station_position.easting(Istation)-Event.localization.location(1)).^2+ ...
    (station_position.northing(Istation)-Event.localization.location(2)).^2);
Event.range=num2str(Event.localization.range/1000);
%Event.localization.Nused=length(Ikeep);
Event.position=num2str(Event.localization.location);

%Check that bearing is the same.
bearing=Event.localization.bearings_all(Istation);
if bearing~=str2num(Event.bearing)
    disp('bearings do not match');
    Event.bearing=num2str(Event.localization.bearings_all(Istation));
end


end