%function loc_index=-plot_movie_all_auto_withmap(tabs_range,DASAR_coords,auto,Isite,mask_locations)
%%% Plot locations of all automated detections, without reference to manual
%%% results
%%   Includes an option select localizations from a map, and return
%%   associated indicies...
%   tabs_range=2 element vector of datenumbers
%mask_locations: contains the following fields
%           .on:  If one, apply mask
%           .lats:  latitudes of turning points
%           .longs: longs of turning points
%  Outputs:
%       Iplotted:  Indicies of auto.locations that occur within time range
%       and are localized

function Iplotted=plot_movie_all_auto_withmap(tabs_range,DASAR_coords,auto,Isite,mask_locations)

Iplotted=[];
if isempty(auto)
    return
end
if isempty(auto.locations)||isempty(auto.locations{1})
    return
end
Zonenumber=[5 6 6 6 7 6];
%Ndasar=[11 7 7 7 7 26]; 
plotMarkerSize=6;
color1=[0 0 0];

m_line(DASAR_coords(:,1),DASAR_coords(:,2),'clip','point','marker','^','MarkerSize',8, 'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
if isempty(auto.locations_ctime)
    return
end
auto_tabs=zeros(1,size(auto.locations_ctime{1},1));
Nplotted=auto_tabs;
for Iauto=1:size(auto.locations_ctime{1},1)
    if  isfield(auto.locations{1}{Iauto}.position,'ctime')
        auto_tabs(Iauto)=datenum(1970,1,1,0,0,auto.locations{1}{Iauto}.position.ctime);
    end
end



%%Plot whale data in lat/long
Icall_auto=find(auto_tabs>=tabs_range(1)&auto_tabs<=tabs_range(2));
K=1;
for I=1:length(Icall_auto)
    try

        %Icall_auto=Nlocs_all_match(Icall(I));
        %Iwant2=find(auto.locations_ctime{1}(Icall_auto(I),:)>0);
        if  isfield(auto.locations{1}{Icall_auto(I)}.position,'location')
            Nplotted(Icall_auto(I))=1;
            VM2=[auto.locations{1}{Icall_auto(I)}.position.location];
            % A2=auto.locations{1}{Icall_auto(I)}.position.major;
            %B2=auto.locations{1}{Icall_auto(I)}.position.minor;
            %Baxis2=auto.locations{1}{Icall_auto(I)}.position.ellipse_ang;
            [lat(K),lon(K)]=UTMtoLL(23,VM2(2),VM2(1),Zonenumber(Isite),'W');
            K=K+1;
        end

        %Place on a mask..

    catch
        disp('Localization failed')
    end

end
Iplotted=find(Nplotted>0);
%%Mask locations from distant airguns by removing anything above a cutoff
%%line

if exist('mask_locations','var')
    if mask_locations.on==1
        lat_mask=interp1(mask_locations.longs(1,:),mask_locations.lats(1,:),lon);
        Igoodd=find(lat<lat_mask);
        lat=lat(Igoodd);lon=lon(Igoodd);
        
    end
end

if K>1

    hh=m_line(lon,lat,'clip','point','marker','+','MarkerSize',plotMarkerSize,'MarkerFaceColor',color1,'MarkerEdgeColor',color1,'LineStyle','none');


end
%xlim(xlimm);ylim(ylimm);
hold on;



end



