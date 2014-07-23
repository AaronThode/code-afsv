function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params)
success_flag=0;
%list_names='S510G0T20100831T000000_BearingInterval_Huber_FilteredLocations';
[list_names,filter_params,station_pos]=load_bowhead_detector_params('.');

data=load(list_names);

%%Create annotation file names for output
[Defaults.Description, Defaults.Template, edit_fields]	=	load_default_annotation_template();
Defaults.Events		=	Defaults.Template;

for I=1:length(data.goodName)
    annotation_names{I}	=	[data.goodName{I} '-notes-BowheadAutomated.mat'];    
    GUI_params	=	[];
    Data{I}				=	Defaults;
    
end

%%Set flags for filtering
filter_position=0;
if isfield(filter_params,'northing')&&isfield(filter_params,'easting')
    filter_position=1;  %Filter locations by position
end

Igood= false(1,length(data.locations));

%%Restrict annotations to certain UTM positions bounded by filt_param
%%'northing' and 'easting' restrictions.
if filter_position
    poss=zeros(length(data.locations),2);

    for I=1:length(data.locations)
        if ~strcmp(data.locations{I}.position.outcome,'successful')
            continue
        end
        
        %Does succesfull localization fit?
        pos=data.locations{I}.position.location;
        poss(I,:)=pos;
        Igood(I)= pos(1)>=filter_params.easting(1) && pos(1)<=filter_params.easting(2);
        Igood(I)= Igood(I) && pos(2)>=filter_params.northing(1) && pos(2)<=filter_params.northing(2);
        
    end
    
    % Check our filtering
    %plot(poss(:,1),poss(:,2),'kx');hold on;grid on; axis('equal');
    %plot(poss(Igood,1),poss(Igood,2),'ro');hold on;grid on; axis('equal');
    
end

%%Finish filtering locations
data.locations=data.locations(Igood);
data.locations_ctime=data.locations_ctime(Igood,:);

%Now convert each location into a series of annotations
Ns=length(data.goodName);

for J=1:Ns
    Icount=1;
    for I=1:length(data.locations)
        newEvent=createEvent(data.locations{I},I,J,Defaults.Template,annotation_names);
        if isempty(newEvent)
            continue
        end
        Data{J}.Events(Icount)=newEvent;
        Icount=Icount+1;
    end
    
end

end

function newEvent=createEvent(location,Ihash,Istation,Template)

%start_time	=	handles.tdate_start...
%    +	datenum(0,0,0,0,0,min(Times));
%min_freq	=	(1000*min(Freq));
%max_freq	=	(1000*max(Freq));
%duration	=	(abs(Times(2) - Times(1)));
%sig_type='FM';
newEvent=[];

start_time=location.ctime_min(Istation);
if start_time==0
    return
end
min_freq=location.Totalfmin(Istation);
max_freq=location.Totalfmax(Istation);
duration=location.Totalduration(Istation);


newEvent=Template;
newEvent.start_time	=	datenum(1970,1,1,0,0,start_time,annotation_names);
newEvent.sig_type		=	sig_type;
newEvent.min_freq		=	min_freq;
newEvent.max_freq		=	max_freq;
newEvent.duration		=	duration;
newEvent.hash_tag       =   now+datenum(0,0,0,0,0,Ihash);  %Just need a unique number
newEvent.link_names    =    annotation_names;

end
