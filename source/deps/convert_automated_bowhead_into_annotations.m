function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params,station_position)
%function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params)
%  Creates an annotation event from a location structure
%  Note that the Event will have 'automated' and 'localization' fields,
%   which are not present in a standard manual annotation.


try
    
    %debug statements below...
    %list_names='S510G0T20100831T000000_BearingInterval_Huber_FilteredLocations';
    %[list_names,filter_params,station_position]=load_bowhead_detector_params('.');
    
    data=load(list_names);
    
    %%Create annotation file names for output
    [Defaults.Description, Defaults.Template, edit_fields]	=	load_default_annotation_template();
    
    %Add custom fields
    Nd=length(Defaults.Description);
    Defaults.Description{Nd+1}='Bowhead Automated Fields';
    Defaults.Description{Nd+2}='Localization parameters';
    Defaults.Description{Nd+3}='Bearing (deg)';
    Defaults.Description{Nd+4}='Range (km)';
    
    Ne=length(edit_fields);
    edit_fields{Ne+1}='bearing';
    edit_fields{Ne+2}='range';
    
    Defaults.Template.automated=[];
    Defaults.Template.localization=[];
    Defaults.Template.bearing=0;
    Defaults.Template.range=0;
    
    Defaults.Events		=	Defaults.Template;
    Defaults.edit_fields=edit_fields;
    
    GUI_params	=	[];
    
    Ns=length(data.goodName);
    keyword=filter_params.keyword;
    for J=1:Ns
        annotation_names{J}	=	sprintf('%s-notes-BowheadAutomated-%s.mat',data.goodName{J},keyword);
        Data_all{J}				=	Defaults;
        
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
    
    %We cycle through each station at each location to make hashtag processing
    %easier.
    Icount=ones(Ns,1);
    for I=1:length(data.locations)
        hashtags=-1*ones(Ns,1);
        isPresent=false(Ns,1);
        %First create each raw event, including assigning a hashtag
        for J=1:Ns
            
            newEvent=createEvent(data.locations{I},I,J,Defaults.Template,annotation_names,station_position);
            if isempty(newEvent)
                continue
            end
            Data_all{J}.Events(Icount(J))=newEvent;
            Icount(J)=Icount(J)+1;
            hashtags(J)=newEvent.hash_tag;
            isPresent(J)=true;
        end
        
        %Now that hashtag created, link common events together.  Stations not
        %involved have hashtags=-1;  This permits the same index to be used
        %for a station, regardless of the event.
        %  Although assigning a 'annotation names' to each event is wasteful
        %  of memory, may be useful in case individual events are modified.
        
        for J=1:Ns
            if isPresent(J)
                Data_all{J}.Events(Icount(J)-1).link_hashtags=hashtags;
                Data_all{J}.Events(Icount(J)-1).localization.station_position=station_position;
            end
            
        end
     
    end  % I
    
    %Now write events to output file.  Check where folder location is...
    
    
    for J=1:Ns
        Data=Data_all{J};
        %Data.localization.station_position=station_position; %Store DASAR location data in annotation for future plotting.
        save(annotation_names{J},'Data','GUI_params');
        
    end
    success_flag=1;
catch
    success_flag=0;  
end
end

function newEvent=createEvent(location,Ihash,Istation,Template,annotation_names,station_position)

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
newEvent.start_time	=	datenum(1970,1,1,0,0,start_time);
newEvent.sig_type		=	'FM';
newEvent.min_freq		=	min_freq;
newEvent.max_freq		=	max_freq;
newEvent.duration		=	duration;
newEvent.hash_tag       =   now+datenum(0,0,Ihash,0,0,0);  %Just need a unique number within a given event
newEvent.link_names    =    cell2mat(annotation_names');
       
%Copy over rest of automated data into a separate 'automated' field
names=fieldnames(location);
for Iname=1:length(names)
   switch names{Iname}
       case 'feature'
           newEvent.automated.feature=location.feature(Istation);
       case 'equalization'
          
       case 'position'
           newEvent.localization=location.position;
           
           %Since a position is guaranteed for this analysis, can compute
           %range.  But just in case...
           
           try
               newEvent.localization.range=sqrt((station_position.easting(Istation)-location.position.location(1)).^2+ ...
                   (station_position.northing(Istation)-location.position.location(2)).^2);
               newEvent.range=newEvent.localization.range/1000; %km
               newEvent.bearing=location.bearing(Istation);
               newEvent.localization.bearings_all=location.bearing;  %Store all bearings for plotting..
           catch
              fprintf('convert_automated_bowhead_into_annotations: You cannot assign a range to this event: no position associated with this detection.\n'); 
           end
           
           
           
       otherwise
           newEvent.automated.(names{Iname})=location.(names{Iname})(Istation);
   end
    
    
end

end
