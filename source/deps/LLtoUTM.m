function [UTMNorthing,UTMEasting,ZoneNumber,ZoneLetter]=LLtoUTM(Iellipsoid,  Lat,  Long)
%function [Northing,Easting,UTMZone,UTMLetter]=LLtoUTM(Iellipsoid,  Lat,  Long)
% Iellipsoid should be 23 to use WGS-84
% //LatLong- UTM conversion.cpp
% //Lat Long - UTM, UTM - Lat Long conversions
% %#include "constants.h"
% #include "LatLong-UTMconversion.h"
% /.*Reference ellipsoids derived from Peter H. Dana's website- 
% http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
% Department of Geography, University of Texas at Austin
% Internet: pdana@mail.utexas.edu
% 3/22/95
% 
% Source
% Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
% 1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency
% //converts lat/long to UTM coords.  Equations from USGS Bulletin 1532 
% //East Longitudes are positive, West longitudes are negative. 
% //North latitudes are positive, South latitudes are negative
% //Lat and Long are in decimal degrees
% 	//Written by Chuck Gantz- chuck.gantz@globalstar.com

pi = 3.14159265;
FOURTHpi = pi / 4;
deg2rad = pi / 180;
rad2deg = 180.0 / pi;

Ellipsoid(23).name='WGS-84';
Ellipsoid(23).EquatorialRadius=6378137;
Ellipsoid(23).eccentricitySquared=0.00669438;

a = Ellipsoid(Iellipsoid).EquatorialRadius;
eccSquared = Ellipsoid(Iellipsoid).eccentricitySquared;
k0 = 0.9996;


%//Make sure the longitude is between -180.00 .. 179.9
LongTemp = (Long+180)-floor((Long+180)/360).*360-180; %// -180.00 .. 179.9;

LatRad = Lat.*deg2rad;
LongRad = LongTemp.*deg2rad;

ZoneNumber = floor((LongTemp + 180)/6) + 1;

if( Lat >= 56.0 & Lat < 64.0 & LongTemp >= 3.0 & LongTemp < 12.0 )
    ZoneNumber = 32;
end


strr='CDEFGHJKLMNPQRSTUVWX';
JJJ=floor(min([ 84+80 Lat+80])/8);
ZoneLetter=strr(JJJ);

%Each zone is divided into horizontal bands spanning 8 degrees of latitude. 
% These bands are lettered, south to north, beginning at 80° S with the 
%letter C and ending with the letter X at 84° N. 
% The letters I and O are skipped to avoid confusion with the numbers one and zero. 
% The band lettered X spans 12° of latitude.



%// Special zones for Svalbard
if( Lat >= 72.0 & Lat < 84.0 )  
    if(      LongTemp >= 0.0  & LongTemp <  9.0 )
        ZoneNumber = 31;
    elseif( LongTemp >= 9.0  & LongTemp < 21.0 )
        ZoneNumber = 33;
    elseif( LongTemp >= 21.0 & LongTemp < 33.0 )
        ZoneNumber = 35;
    elseif( LongTemp >= 33.0 & LongTemp < 42.0 )
        ZoneNumber = 37;
    end
end

LongOrigin = (ZoneNumber - 1).*6 - 180 + 3;  %//+3 puts origin in middle of zone
LongOriginRad = LongOrigin .* deg2rad;

%//compute the UTM Zone from the latitude and longitude
%sprintf(UTMZone, '%d%c', ZoneNumber, UTMLetterDesignator(Lat));

eccPrimeSquared = (eccSquared)./(1-eccSquared);

N = a./sqrt(1-eccSquared.*sin(LatRad).*sin(LatRad));
T = tan(LatRad).*tan(LatRad);
C = eccPrimeSquared.*cos(LatRad).*cos(LatRad);
A = cos(LatRad).*(LongRad-LongOriginRad);

M = a.*((1	- eccSquared./4		- 3.*eccSquared.*eccSquared./64	- 5.*eccSquared.*eccSquared.*eccSquared./256).*LatRad - ...
    (3.*eccSquared./8	+ 3.*eccSquared.*eccSquared./32	+ 45.*eccSquared.*eccSquared.*eccSquared./1024).*sin(2.*LatRad) + ...
    (15.*eccSquared.*eccSquared./256 + 45.*eccSquared.*eccSquared.*eccSquared./1024).*sin(4.*LatRad) - ...
    (35.*eccSquared.*eccSquared.*eccSquared./3072).*sin(6.*LatRad));

UTMEasting = (k0.*N.*(A+(1-T+C).*A.*A.*A./6 + (5-18.*T+T.*T+72.*C-58.*eccPrimeSquared).*A.*A.*A.*A.*A./120)+ 500000.0);

UTMNorthing = (k0.*(M+N.*tan(LatRad).*(A.*A./2+(5-T+9.*C+4.*C.*C).*A.*A.*A.*A./24 + (61-58.*T+T.*T+600.*C-330.*eccPrimeSquared).*A.*A.*A.*A.*A.*A./720)));
if(Lat < 0)
    UTMNorthing = UTMNorthing+10000000.0; %././10000000 meter offset for southern hemisphere
end

