function [sensor_name,bits,Vmax]=parse_instrument_string(namee)
%DASAR_5Volt_24bit_DASARsensor

Islash=strfind(namee,'_');
Iend=strfind(namee,'sensor')-1;
sensor_name=namee((Islash(end)+1):Iend);

try
    Iend=strfind(namee,'Volt')-1;
    Vmax=str2double(namee((Islash(1)+1):(Iend)));

    Iend=strfind(namee,'bit')-1;
    bits=str2double(namee((Islash(2)+1):(Iend)));
catch
    bits=[];
    Vmax=[];
end

end

