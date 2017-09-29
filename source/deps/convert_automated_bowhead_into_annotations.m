function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params,station_position)
%function success_flag=convert_automated_bowhead_into_annotations(list_names,filter_params,station_position)
%  Creates an annotation event from a location structure
%  Note that the Event will have 'automated' and 'localization' fields,
%   which are not present in a standard manual annotation.
%  if filter_params.keyword=''all'' then no position filtering...


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
    %Defaults.Description{Nd+4}='Kappa';
    Defaults.Description{Nd+4}='Range (km)';
    Defaults.Description{Nd+5}='Position [easting northing] (m)';
    Defaults.Description{Nd+6}='Index of array position';
    
    Ne=length(edit_fields);
    edit_fields{Ne+1}='bearing';
    edit_fields{Ne+2}='range';
    edit_fields{Ne+3}='position';
    
    Defaults.Template.automated=[];
    Defaults.Template.localization=[];
    Defaults.Template.bearing=0;
    %Defaults.Template.kappa=0;
    Defaults.Template.range=0;
    Defaults.Template.position=[0 0];
    Defaults.Template.Istation=1;
    Defaults.Template.array_index=[];
    
    
    Defaults.Events		=	Defaults.Template;
    Defaults.edit_fields=edit_fields;
    
    GUI_params	=	[];
    
    Ns=length(data.goodName);
    keyword=filter_params.keyword;
    
    %Create linked annotation names
    for J=1:Ns
        annotation_names{J}	=	sprintf('%s-notes-BowheadAutomated-%s.mat',data.goodName{J},keyword);
        Data_all{J}				=	Defaults;
        
    end
    
    %%Set flags for filtering
    filter_position=0;
    if isfield(filter_params,'northing')&&isfield(filter_params,'easting')&&~(strcmp(lower(keyword),'all'))
        filter_position=1;  %Filter locations by position (rectangle)
    end
    
    if isfield(filter_params,'range')&&isfield(filter_params,'UTM_center')&&~(strcmp(lower(keyword),'all'))
        filter_position=2;  %Filter locations by position (circle)
        rangee=-1*ones(1,length(data.locations));
    end
    
    Igood= true(1,length(data.locations));
    
    %%Restrict annotations to certain UTM positions bounded by filt_param
    %%'northing' and 'easting' restrictions.
    if filter_position==2
        Igood= false(1,length(data.locations));
    
        poss=zeros(length(data.locations),2);
        h=waitbar(0,'Filtering positions.. wait');
        for I=1:length(data.locations)
            
            if rem(I,500)==0
                waitbar(I/length(data.locations),h);
                
            end
            if ~strcmp(data.locations{I}.position.outcome,'successful')
                continue
            end
            
            fmin=min([data.locations{I}.feature.robust_fmin]);
            Igood(I)=(fmin>=filter_params.freq_range(1))&(fmin<=filter_params.freq_range(2));
            %Does succesfull localization fit?
            pos=data.locations{I}.position.location;
            poss(I,:)=pos;
            tmp=[pos(1)-filter_params.UTM_center(1) pos(2)-filter_params.UTM_center(2)];
            rangee(I)=sqrt(sum(tmp.^2));
            Igood(I)=Igood(I) && (filter_params.range>=rangee(I));
            
            
        end
        
        % Check our filtering
        %plot(poss(:,1),poss(:,2),'kx');hold on;grid on; axis('equal');
        %plot(poss(Igood,1),poss(Igood,2),'ro');hold on;grid on; axis('equal');
        close(h)
    elseif filter_position==1
        Igood= false(1,length(data.locations));
    
        poss=zeros(length(data.locations),2);
        h=waitbar(0,'Filtering positions.. wait');
        
        for I=1:length(data.locations)
            
            if rem(I,500)==0
                waitbar(I/length(data.locations),h);
                
            end
            if ~strcmp(data.locations{I}.position.outcome,'successful')
                continue
            end
            
            %Does succesfull localization fit?
            pos=data.locations{I}.position.location;
            poss(I,:)=pos;
            Igood(I)= pos(1)>=filter_params.easting(1) && pos(1)<=filter_params.easting(2);
            Igood(I)=Igood(I)&& (fmin>=filter_params.freq_range(1))&&(fmin<=filter_params.freq_range(2));
           
            Igood(I)= Igood(I) && pos(2)>=filter_params.northing(1) && pos(2)<=filter_params.northing(2);
            
        end
        
        % Check our filtering
        %plot(poss(:,1),poss(:,2),'kx');hold on;grid on; axis('equal');
        %plot(poss(Igood,1),poss(Igood,2),'ro');hold on;grid on; axis('equal');
        close(h)
    end
    
    
    %%Finish filtering locations
    data.locations=data.locations(Igood);
    data.locations_ctime=data.locations_ctime(Igood,:);
    
    %Now convert each location into a series of annotations
    
    %We cycle through each station at each location to make hashtag processing
    %easier.
    Icount=ones(Ns,1);
    h=waitbar(0,'Please wait');
    
    for I=1:length(data.locations)
        hashtags=-1*ones(Ns,1);
        isPresent=false(Ns,1);
        %First create each raw event, including assigning a hashtag
        
        if rem(I,100)==0
            waitbar(I/length(data.locations),h);
            
        end
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
    close(h)
    %Now write events to output file.  Check where folder location is...
    
    %Check folder location, output in same directory as bowhead
    %detections..
    
    
    for J=1:Ns
        Data=Data_all{J};
        %Data.localization.station_position=station_position; %Store DASAR location data in annotation for future plotting.
        save(annotation_names{J},'Data','GUI_params');
        
    end
    success_flag=1;
catch
    success_flag=0;
    MException.last
    %lasterr.identifier
    %lasterr.stack
end
end

function newEvent=createEvent(location,Ihash,Istation,Template,annotation_names,station_position)
%%Turn a location object into an annotation event object

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
                
                %Reject all bearings that do not contribute to a
                %localization...
                newEvent.localization.bearings_all=NaN*ones(size(location.bearing));
                newEvent.localization.kappa=NaN*ones(size(location.bearing));
                Ikeep=location.position.Ikeep;
                
                newEvent.localization.bearings_all(Ikeep)=location.bearing(Ikeep);  %Store all bearings for plotting..
                newEvent.localization.kappa(Ikeep)=location.kappa(Ikeep);  %Store all kappa for plotting..
                
                newEvent.bearing=newEvent.localization.bearings_all(Istation);
                newEvent.position=location.position.location;  %This is an editable field
                newEvent.Istation=Istation; %Useful to identify where in bearings_all we are.
                
                 newEvent.localization.range=sqrt((station_position.easting(Istation)-location.position.location(1)).^2+ ...
                    (station_position.northing(Istation)-location.position.location(2)).^2);
                newEvent.range=newEvent.localization.range/1000; %km
               
            catch
                fprintf('convert_automated_bowhead_into_annotations: You cannot assign a range to this event: no position associated with this detection.\n');
            end
            
            
            
        otherwise
            temp=location.(names{Iname});
            if length(temp)>=Istation
                newEvent.automated.(names{Iname})=location.(names{Iname})(Istation);
            else
                newEvent.automated.(names{Iname})=location.(names{Iname})(end);
                
            end
    end
    
    
end

end
