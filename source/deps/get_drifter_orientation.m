%%%%get_drifter_orientation.m%%%%
%  [brefa,prefa]=get_drifter_orientation(filename);
% Input:
%    filename: full filename including pathname (This allows drifter ID to
%    be identified from folder).
% Output:
%   brefa: azimuthal bearing correction.
%   prefa: tilt correction
function [brefa,prefa]=get_drifter_orientation(filename)
%keyboard
brefa=0;prefa=0;

ind=strfind(filename,'drifter-');
if length(ind)>1
    error('get_drifter_orientation:  filename contains drifter- more than once');
end
%%datestr
ind=ind+length('drifter-5V')+1;
tabs_txt=filename(ind:(ind+14));
Year=str2double(tabs_txt(1:4));
Day=str2double(tabs_txt(6:8));
hr=str2double(tabs_txt(10:11));
minn=str2double(tabs_txt(12:13));
sec=str2double(tabs_txt(14:15));
tabs=datetime(Year,1,Day,hr,minn,sec);

ind=strfind(filename,'Drifter')+length('Drifter');
DrifterID=str2double(filename(ind));

%%Add magnetic declination
if tabs>datetime(2024,6,1,0,0,0)& tabs<datetime(2024,8,1,0,0,0)  %June 2024 seamounts
        brefa=brefa-15-16/60;
elseif tabs>datetime(2023,10,1,0,0,0)& tabs<datetime(2023,11,1,0,0,0)  %Oct 2023 Kelvin
    brefa=brefa-15-14/60;

else %San Diego
    brefa=brefa+11+7/60;
end

switch DrifterID
    case 6


end