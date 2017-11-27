function [list_names,filter_params,station_locations]=load_bowhead_detector_params(outputdir,GSI_location_dir_template)
%function [list_names,filter_params]=load_bowhead_detector_params(outputdir)
% list_names is a cell array of full pathnames to automated detection files
% filter_params is a structure with two-element arrays for filtering
%   automated results.  'easting' and 'northing' are built in
% station_locations:  structure with 'northing' and 'easting' fields for
%       each DASAR, from 'A'(1)  to 'G'

list_names=[];
filter_params=[];
station_locations=[];

% Location folder with automated results to convert...
dirname = uigetdir(outputdir, 'Select a folder with automated bowhead detector results...');

if dirname==0
    return
end

mydir=pwd;
cd(dirname);

% Define a template for selecting files

%%Attempt to guess site and year from directory name; assumme a substring
%%like 'Shell14_Site5' exists near end of string...
Ishell=max(strfind(dirname,'Shell'))+length('Shell');
if isempty(Ishell)
    uiwait(msgbox(sprintf('Cannot find string ''Shell'' in %s',dirname)));
    
    prompt={'Year (YYYY)','Site'};
    name='Year and Site';
    numlines=1;
    options.Resize='on';
    options.WindowStyle='normal';
    defaultanswer={'2010','5'};
    answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    year=eval(answer{1});
    Site=eval(answer{2});
    
    % temp=input('Enter year and Site [year Site]:');
    % year=temp(1);
    %Site=temp(2);
else
    

    year=num2str(2000+str2num(dirname(Ishell+(0:1))));
    Site=str2num(dirname(Ishell+7));
end


GSI_location_dir=sprintf(GSI_location_dir_template,year,year);

prompt={'Enter a template string for a location output file:', ...
    'Enter geographical restriction (Range, West,Center,East), 2x2 UTC matrix, or ''All'':', ...
    'UTM location to use if ''Range'' selected above (''2010 array'',''2012 array'',''2014 array'',or two-element vector):', ...
    'Range (km) [used only if ''Range'' selected above]:','Min Frequency range (Hz)'};
name='File Template';
numlines=1;
options.Resize='on';
options.WindowStyle='normal';
%defaultanswer={sprintf('S%i*_Huber_FilteredLocations.mat',Site),'All',sprintf('%s array',year),'10','[30 200]'};
defaultanswer={sprintf('S%i*_Huber_FilteredLocations.mat',Site),'West','Site5','3','[30 200]'};

answer=inputdlg(prompt,name,numlines,defaultanswer,options);
fnames=dir(answer{1});
loc_keyword=answer{2};
UTM_keyword=answer{3};
filter_params.range=eval(answer{4})*1000;
filter_params.freq_range=eval(answer{5});

filter_params.keyword=loc_keyword;
filter_params.UTM_keyword=UTM_keyword;

% Look for possible matches in folder, and select all relevent files
for I=1:length(fnames)
    list_names{I}=fnames(I).name;
end

if isempty(fnames)
    cd(mydir);
    uiwait(msgbox('No files found!','modal'));
    return
end


[Sel, OK]	=	listdlg('ListString', list_names,...
    'Name', 'Available files',...
    'PromptString', 'Select file(s) to load:',...
    'OKString', 'Load',...
    'CancelString', 'New', ...
    'ListSize',[320 300]);

list_names=list_names(Sel);
%%Get site number

for I=1:length(list_names)
    Site_number(I)=str2num(list_names{I}(2));
    list_names{I}=fullfile(dirname,list_names{I});
    
end
cd(mydir);

station_locations=get_Site_data(station_locations,GSI_location_dir);


%If range restriction, determine it...
[filter_params,isrange]=get_range_restriction(loc_keyword,filter_params);
if isrange
    station_locations.extra=filter_params.UTM_center;
    return
end

%What DASARS used to compute box?
IDASAR=logical([1 1 1 1 1 1]);  %DASARS CDEG

station_locations.northing=station_locations.northing(IDASAR);
station_locations.easting=station_locations.easting(IDASAR);
%Site5_northing=Site5_northing(IDASAR);
%Site5_easting=Site5_easting(IDASAR);

xbuffer=100; %A little extra shrinking of the center site width, in meters (makes sure box of Center is inside of all DASARS
filter_params.northing=[min(station_locations.northing)+xbuffer max(station_locations.northing)-xbuffer];  %Boxes at same latitude
dl=abs(min(station_locations.easting)-max(station_locations.easting));  %How wide should a box be?  For now, make width the same as the width of the 'center' box..
%dl=10000; %10 km box width

filter_params.easting=[min(station_locations.easting)+xbuffer max(station_locations.easting)-xbuffer];

switch(lower(loc_keyword))
    case 'center'
        %Keep the same
    case 'west'
        filter_params.easting=filter_params.easting-(2*xbuffer+dl);
    case 'east'
        filter_params.easting=filter_params.easting+(2*xbuffer+dl);
    case 'all'
        filter_params.easting=[];
    otherwise
        return
        
end
%%%%inner function get_Site_data
    function station_locations=get_Site_data(station_locations,GSI_location_dir)
        %%% try to load %%DASAR_locations_YYYY.mat
        %%%% format will be Site{5}.easting/northing
        
        
        %find year and Site from list_names{1}
        try
            if Site==0
                iSite=6;
            else
                iSite=Site;
            end
            locs=load(GSI_location_dir);
            
            mySite=locs.Site{iSite};
            Igood=find(mySite.northing>0);
            station_locations.northing=mySite.northing(Igood);
            station_locations.easting=mySite.easting(Igood);
        catch
            uiwait(warndlg('load_bowhead_detector_params: Can''t access Site location information, using default locations for Site 5'));
            % Determine UTC boundaries
            Site5_easting=[
                412659
                418898
                412949
                419139
                413253
                413492];
            
            Site5_northing=[
                7795035
                7798311
                7802045
                7805333
                7809065
                7816110];
            
            Site4_easting=[
                560332
                554220
                560195
                554068
                560066
                553948
                559900
                541824
                541960
                554007
                560054
                554185
                560126];
            
            Site4_northing=[
                7794371
                7798004
                7801398
                7805012
                7808391
                7812033
                7815436
                7811774
                7804792
                7819070
                7822482
                7826006
                7829492];
            
            %Save data to station_locations output variable
            switch iSite(1)
                case 5
                    station_locations.northing=Site5_northing;
                    station_locations.easting=Site5_easting;
                    
                case 4
                    station_locations.northing=Site4_northing;
                    station_locations.easting=Site4_easting;
                    
                otherwise
                    uiwait(errordlg('Can''t access Site location information'));
            end
        end
        
    end

%%%%%%%%%%
%inner function-get_range_restriction
    function [filter_params,isrange]=get_range_restriction(loc_keyword,filter_params)
        isrange=false;
        switch(lower(loc_keyword))
            case 'range'
                isrange=true;
                switch(lower(UTM_keyword))
                    case '2010 array'
                        filter_params.UTM_center=eval('[4.174860699660919e+05 7.817274204098196e+06]');
                    case '2012 array'
                        filter_params.UTM_center=eval('[4.213046781419147e+05 7.809012407783064e+06]');
                        %                 VA deployment - August 20th 2012
                        % water depth : 53.1 meters
                        %
                        % horizontal separation: 210 feet (between anchor and additional recorder)
                        %
                        % 17:26 - Floats in water - WP 136 (Susanna GPS)
                        % N 70 22.486	- W 143 06.000
                        %
                        % 17:27 - Unit 2 of VA in water - WP 137
                        % N 70 22.480 - W 143 05.984
                        %
                        % **17:28 - 1st Anchor (bottom of vertical array)  in water - WP 138
                        % N 70 22.479 - W 143 05.984
                        %
                        % 17:29 - Autonomous Recorder in water (Unit 2) (210 feet away from anchor) - WP 139
                        % N 70 22.473 - W 143 05.998
                        %
                        % 17:32 - 2nd Anchor in water  - WP 140
                        % N 70 22.415 - W 143 05.832
                        %
                        % Calibrations made on August 20th 2012
                        %
                        % C24 - 17:53:05
                        % 70 24.36014,N,143 02.72971,W
                        %
                        % C14 - 18:29:25
                        % 70 22.50987,N,143 12.44420,W
                        %
                        % C22 - 19:27
                        % 70 20.60442,N,143 02.68435,W
                        %
                        % CTD cast done at DASAR D. Alex didn't extract the data from the SBE yet.
                        
                        %WP 138
                        
                        
                        
                    case '2014 array'
                        %                 Date: August 16, 2014
                        %
                        % Time: 20:30
                        %
                        % Three GPS readings:
                        %
                        %             Buoys in water:  N_70°20.532 W_143°15.233 time  20:35:41 waypoint_450
                        %
                        %             Bruce released:  N_70°20.504 W_143°15.532 time 20:38:56  waypoint_451
                        %
                        %             Danforth released: N70°20.493 W143°15.756 time 20:41:53 waypoint 452
                        %
                        % Depth from Tiger: 47.4 m  Wind speed (avg): 12 m/s  Water temp: 4.6 C
                        filter_params.UTM_center=eval('[4.152065807469279e+05 7.805557380281989e+06]');
                        
                    case 'site0'
                        filter_params.UTM_center=[5.3164e+05 7.8128e+06];
                    otherwise
                        filter_params.UTM_center=eval(UTM_keyword);
                end
                
        end
        
    end
end

