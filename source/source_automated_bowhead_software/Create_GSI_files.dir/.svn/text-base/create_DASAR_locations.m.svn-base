function create_DASAR_2013_locations

path(path,'../ComputerSpecificScripts.dir');

%[rawdatadir]=load_pathnames('Shell13');
rawdatadir='/Users/thode/Jonah/Shell2013_GSI_Data';
letter='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
strdate='0915';
Zonenumber=[5 6 6 6 7 6];
Ndasar=[7 26 7 15 7 26];
head=0;

for Isite=1:6
    easting=zeros(Ndasar(Isite),1);
    northing=easting;
    for Id=1:Ndasar(Isite)
        clear head
        for JJ=0:1
            
            if Isite==6
                template_str=sprintf('%s/S%i13gsif/S%i13%s%i',rawdatadir,0,0,letter(Id),JJ);
                
            else
                template_str=sprintf('%s/S%i13gsif/S%i13%s%i',rawdatadir,Isite,Isite,letter(Id),JJ);
            end
            try
                fn=dir(sprintf('%s/*%sT*',template_str,strdate));
                if length(fn)==0
                    continue
                end
                head=readgsif_header([template_str '/' fn.name]);
                fprintf('DASAR %i%s  exists in %s!!\n',Isite,letter(Id),fn.name);
                
                easting(Id)=head.UTMX;
                northing(Id)=head.UTMY;
            catch
                fprintf('DASAR %i%s not existing\n',Isite,letter(Id));
            end
        end
    end
    
    
    [easting northing]
    tmp=makesite(easting,northing,Isite);
    
    
    lat=zeros(1,Ndasar(Isite));
    lon=lat;
    for Id=1:Ndasar(Isite)
        
        [lat(Id),lon(Id)]=UTMtoLL(23,northing(Id),easting(Id),Zonenumber(Isite),'W');
    end
    Site{Isite}=tmp;
    Site{Isite}.lat=lat;
    Site{Isite}.lon=lon;
    
end

save DASAR_locations_2013 Site
end

function Site=makesite(easting,northing,I)
strr='abcdefghijklmnopqrstuvwxyz';

Site.easting=easting;
Site.northing=northing;
%Site=rmfield(Site,'stations');
for J=1:length(Site.easting),
    Site.stations{J}=strr(J);
end
end
