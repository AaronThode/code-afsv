function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params)
%function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params)
%  Creates an annotation event from a location structure
%  Note that the Event will have 'automated' and 'localization' fields,
%   which are not present in a standard manual annotation.
success_flag=1;



%debug statements below...
%list_names='S510G0T20100831T000000_BearingInterval_Huber_FilteredLocations';
%[list_names,filter_params,station_position]=load_bowhead_detector_params('.');
try
data=load(list_names);

%%Create annotation file names for output
[Defaults.Description, Defaults.Template, edit_fields]	=	load_default_annotation_template();

%Add custom fields
Nd=length(Defaults.Description);
Defaults.Description{Nd+1}='Airgun Automated Fields';
Defaults.Description{Nd+2}='interval (sec)';
Defaults.Description{Nd+3}='Bearing (deg)';
%Defaults.Description{Nd+4}='Kappa';
%Defaults.Description{Nd+4}='Range (km)';
%Defaults.Description{Nd+5}='Position [easting northing] (m)';
%Defaults.Description{Nd+6}='Index of array position';

Ne=length(edit_fields);
edit_fields{Ne+1}='interval';
edit_fields{Ne+2}='bearing';
%edit_fields{Ne+3}='position';

Defaults.Template.automated=[];
Defaults.Template.interval=[];
Defaults.Template.bearing=0;
%Defaults.Template.kappa=0;
%Defaults.Template.range=0;
%Defaults.Template.position=[0 0];
%Defaults.Template.Istation=1;

Defaults.Events		=	Defaults.Template;
Defaults.edit_fields=edit_fields;

GUI_params	=	[];

%Ns=length(data.goodName);
keyword=data.param.calibration_keyword;

%Create linked annotation names
[PATHSTR,local_names,EXT] = fileparts(list_names);
Idash=min(findstr(local_names,'_'))-1;
goodName=local_names(1:Idash)

%for J=1:Ns
annotation_names{1}	=	sprintf('%s-notes-AirgunAutomated-%s.mat',goodName,keyword);
Data				=	Defaults;


Icount=1;
Npp=length(data.airgun_shots.ctime);
h=waitbar(0,'Please wait')
for J=1:length(data.airgun_shots.ctime)
    if rem(J,500)==0
       waitbar(J/Npp,h); 
        
    end
    newEvent=createEvent(data.airgun_shots,J,Defaults.Template,annotation_names);
    if isempty(newEvent)
        continue
    end
    Data.Events(Icount)=newEvent;
    Icount=Icount+1;
end



    save(annotation_names{1},'Data','GUI_params');

catch
    success_flag=0;
    uiwait(msgbox(sprintf('%s cannot be interpreted as an airgun file',list_names)));
end

end

function newEvent=createEvent(airgun_shots,Ihash,Template,annotation_names)
%%Turn a location object into an annotation event object

%start_time	=	handles.tdate_start...
%    +	datenum(0,0,0,0,0,min(Times));
%min_freq	=	(1000*min(Freq));
%max_freq	=	(1000*max(Freq));
%duration	=	(abs(Times(2) - Times(1)));
%sig_type='FM';
newEvent=[];

start_time=airgun_shots.ctime(Ihash);
if start_time==0
    return
end
min_freq=5;
try
max_freq=airgun_shots.level.max_freq(Ihash);
catch
    max_freq=airgun_shots.level.peakF(Ihash);
end
duration=airgun_shots.level.t_Malme(Ihash);


newEvent=Template;
newEvent.start_time	=	datenum(1970,1,1,0,0,start_time);
newEvent.sig_type		=	'FM';
newEvent.min_freq		=	min_freq;
newEvent.max_freq		=	max_freq;
newEvent.duration		=	duration;
newEvent.hash_tag       =   now+datenum(0,0,Ihash,0,0,0);  %Just need a unique number within a given event
newEvent.link_names    =    cell2mat(annotation_names');
newEvent.interval   =  airgun_shots.ICI(Ihash);
newEvent.bearing    =  airgun_shots.bearing(Ihash);

%Copy over rest of automated data into a separate 'automated' field
names=fieldnames(airgun_shots.level);
for Iname=1:length(names)
    newEvent.automated.(names{Iname})=airgun_shots.level.(names{Iname})(Ihash);
    
end

end
