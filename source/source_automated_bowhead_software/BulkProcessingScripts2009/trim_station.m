function station=trim_station(station,Igood)
%function station=trim_station(station,Igood)

%%%Trim station to contain only indicies present in 'Igood'
%function station=trim_station(station,Igood)

topnames=fieldnames(station);
names = fieldnames(station.feature);

for I=1:length(topnames),
    switch topnames{I}
        case 'feature'
            for JJ=1:length(names),
                station.feature.(names{JJ})=station.feature.(names{JJ})(:,Igood);
            end
        case 'param'
            continue;
        case 'Image'
            station.Image=station.Image(Igood);
        otherwise
            station.(topnames{I})=station.(topnames{I})(:,Igood);

    end
end




end
