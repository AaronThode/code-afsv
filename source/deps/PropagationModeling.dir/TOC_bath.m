%%%%%%%%%TOC_bath.m%%%%%%
% Input: 
%   D: water depth under source in meters
%   case_bath: string with scenario.  Some common cases:
%        'slopeXkmYm':  creates a slope such that X km away water depth is X-Y m (shallows Y meters)
%        'wedgeXkmYkmZm': creates a wedge such that X km away wedge starts, so that by Y km depth is D-Z
% Note: ranges expressed in kilometers in first column; depths are positive downward expressed in meters in second column.
% the first row is depth under source.
% For reciprocity calc: manipulate so a range of 0 m is depth under buoy
%  bath dimensions: [ranges 2];  
%
% function bath=TOC_bath(case_bath,D);
function bath=TOC_bath(case_bath,D)
if strcmp(case_bath,'select_from_menu'),
    menu_chc={'flat','wedge','Ewing_shallow_slope','test_slope','wedge_demo'};
    Ichc=menu('Select a bathymetry profile:',menu_chc)
    case_bath=menu_chc{Ichc};
end

%Check if numbers for slope embedded in case_bath
I1=findstr(case_bath,'slope')+length('slope'); %'slopeXkmYm'
I2=findstr(case_bath,'km')-length('km')+1;
I3=findstr(case_bath,'km')+length('km');
if ~isempty(I1)&&~isempty(I2)
    bath=[0 D;
        str2double(case_bath(I1:I2)) D-str2num(case_bath(I3:(end-1)))];
    
end
switch case_bath
    case 'flat'
        bath=[0 D];
    case {'wedge','wedge_demo'}
        bath=[0 200;
            5 200;
            12.5 0];
    
    case 'slope'
        
        prompt={'vector of ranges in kilometers relative to source','vector of water depths (m)'};
        dlgTitle='Parameters for bathymetry selection';
        def={'linspace(0,18,10)','D*linspace(1,45/55,10)'};
        answer=inputdlg(prompt,dlgTitle,1,def);
        bath=[eval(answer{1})' eval(answer{2})']
        
    case 'Ewing_shallow_slope'
        mydir=pwd;
        cd '/Users/thode/Projects/Airgun_modeling/Calibration_geography/ShallowSiteJune02/Bathymetry_profiles'
        load shallow-range-depth-profile-12km ;
        bath = Darray{1} ;    %  manipulate so a range of 0 m is depth under buoy
        figure;
        plot(bath(:,1),bath(:,2),'x');
        title(sprintf('Source bottom depth %6.2f, Buoy bottom depth %6.2f',source_bottom_depth(1),rcvr_bottom_depth(1)));
        pause;
        cd(mydir);
end