%function loc_index=-plot_movie_all_manual_withmap(tabs_range,DASAR_coords,manual,Isite,mask_locations)
%%% Plot locations of all manual detections
%%% results
%%   Includes an option select localizations from a map, and return
%%   associated indicies...
%   tabs_range=2 element vector of datenumbers
%mask_locations: contains the following fields
%           .on:  If one, apply mask
%           .lats:  latitudes of turning points
%           .longs: longs of turning points
%  Outputs:
%       Iplotted:  Indicies of manual.locations that occur within time range
%       and are localized

function Iplotted=plot_movie_all_manual_withmap(tabs_range,DASAR_coords,manual,Isite,mask_locations)


Zonenumber=[5 6 6 6 7];
Ndasar=[11 7 7 7 7]; plotMarkerSize=6;
color1=[0 0 0];

if isempty(manual)||isempty(manual.localized)
    Iplotted=[];
    return
end
manual_tabs=datenum(1970,1,1,0,0,manual.localized{1}.ctev);
Nplotted=zeros(size(manual_tabs));




%%Plot whale data in lat/long
Icall_manual=find(manual_tabs>=tabs_range(1)&manual_tabs<=tabs_range(2));
K=1;
for I=1:length(Icall_manual)
    try

        %Icall_manual=Nlocs_all_match(Icall(I));
        %Iwant2=find(manual.locations_ctime{1}(Icall_manual(I),:)>0);
        if  findstr(manual.localized{1}.Outcome{Icall_manual(I)}','uccessful')
            Nplotted(Icall_manual(I))=1;
            VM2=[manual.localized{1}.utmx(Icall_manual(I)) manual.localized{1}.utmy(Icall_manual(I))];
            % A2=manual.locations{1}{Icall_manual(I)}.position.major;
            %B2=manual.locations{1}{Icall_manual(I)}.position.minor;
            %Baxis2=manual.locations{1}{Icall_manual(I)}.position.ellipse_ang;
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

if mask_locations.on==1
    lat_mask=interp1(mask_locations.longs(1,:),mask_locations.lats(1,:),lon);
    Igoodd=find(lat<lat_mask);
    lat=lat(Igoodd);lon=lon(Igoodd);

end

m_line(DASAR_coords(:,1),DASAR_coords(:,2),'clip','point','marker','^','MarkerSize',8, 'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');

if K>1
   
    hh=m_line(lon,lat,'clip','point','marker','+','MarkerSize',plotMarkerSize,'MarkerFaceColor',color1,'MarkerEdgeColor',color1,'LineStyle','none');


end
%xlim(xlimm);ylim(ylimm);
hold on;



end



