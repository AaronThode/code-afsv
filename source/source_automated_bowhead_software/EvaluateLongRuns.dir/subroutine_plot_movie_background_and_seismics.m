%%subroutine_plot_movie_background_and_seismics.m%
% Plots everything except whale locations for a movie
%  Called by master_evaluation_complete_movie.m
%% Iput variable 'Iplot' sets which subplot window to use


color1=[0 0 0];color2=[0.35 0.35 0.35];color3=[0.65 0.65 0.65];color4=[0.85 0.85 0.85];
shotcolor=[0 0 0];lastshotcolor=[.25 .5 1];  % Was [0.5 0.5 1]
aerialSightingcolor=[0 1 0]; lastsightcolor=[0 1 0];
shotcolor=[0 0 1];lastshotcolor=[.25 .5 1];  % Was [0.5 0.5 1]
lastnoshotcolor=[0.2500    0.500    0.5000];


figure(100)
%set(gcf,'pos',[366         536        1149         504]);
%set(gcf,'pos',[366         131        1408         909]);
set(gcf,'pos',[ 161         379        1104         553]);
if Iplot==1
    clf
end
% clf = clear current figure window

if exist('manual')&&run_options.auto_only==0
    subplot(2,1,Iplot)
end
set(gca,'Fontsize',10)      % set = set object properties, where we have property, propertyName, propertyValue.



m_proj('stereographic','lon',mapmidLon,'lat',mapmidLat,'rad',maplimits,'rec','on');
m_grid('color',[.3 .3 .3],'linestyle','-');

m_line(mask_locations.longs,mask_locations.lats)


for Iguns=1:3
    
    %Plot_past_active
    if Icount>1
        
        m_line(past_active(Iguns).long,past_active(Iguns).lat,'clip','point','marker','o','MarkerSize',4.5,'MarkerFaceColor',lastshotcolor,'MarkerEdgeColor','none','LineStyle','none');
        if run_options.shooting_only==0
            m_line(past_notactive(Iguns).long,past_notactive(Iguns).lat,'clip','point','marker','o','MarkerSize',4.5,'MarkerFaceColor',lastnoshotcolor,'MarkerEdgeColor','none','LineStyle','none');
            
            
        end
        
    end
    
    %Plot current active
    if exist('seis')
        Ishoot=find(seis.Seismic(Iguns).tabs>=tabs_range(Ihr)&seis.Seismic(Iguns).tabs<=tabs_range(Ihr+1));
        Igun=find(seis.Seismic(Iguns).guns(Ishoot)>0);
        m_line(seis.Seismic(Iguns).long(Ishoot(Igun)),seis.Seismic(Iguns).lat(Ishoot(Igun)),'clip','point','marker','o','MarkerSize',4.5,'MarkerFaceColor','b','MarkerEdgeColor','none','LineStyle','none');
        
        %Plot current_inactive
        Inone=find(seis.Seismic(Iguns).guns(Ishoot)==0);
        if run_options.shooting_only==0
            m_line(seis.Seismic(Iguns).long(Ishoot(Inone)),seis.Seismic(Iguns).lat(Ishoot(Inone)),'clip','point','marker','o','MarkerSize',4.5,'MarkerFaceColor','g','MarkerEdgeColor','none','LineStyle','none');
        end
        
        
        
        if ~exist('manual')|Iplot==2
            past_active(Iguns).long=seis.Seismic(Iguns).long(Ishoot(Igun));
            past_active(Iguns).lat=seis.Seismic(Iguns).lat(Ishoot(Igun));
            past_notactive(Iguns).long=seis.Seismic(Iguns).long(Ishoot(Inone));
            past_notactive(Iguns).lat=seis.Seismic(Iguns).lat(Ishoot(Inone));
        end
    end
    
end

%m_usercoast('alaska_hi','patch',[.7 .7 .7],'linestyle','none');
m_usercoast('alaska_hi','patch',[0 0 0],'linestyle','none');

%m_coast('patch',[.7 .7 .7],'linestyle','none');

m_line(-142.2,71.1,'marker','o','MarkerSize',4.5,'MarkerFaceColor','b','MarkerEdgeColor','none','LineStyle','none');
m_text(-142.1,71.1,' Airgun activity','color','b','fontweight','normal','fontsize',8);

if run_options.shooting_only==0
    m_line(-142.2,71.3,'marker','o','MarkerSize',4.5,'MarkerFaceColor','g','MarkerEdgeColor','none','LineStyle','none');
    m_text(-142.1,71.3,' No activity','color','g','fontweight','normal','fontsize',8);
end

m_line(-142.2,71.2,'marker','o','MarkerSize',4.5,'MarkerFaceColor',lastshotcolor,'MarkerEdgeColor','none','LineStyle','none');
m_text(-142.1,71.2,' Earlier Activity','color',[0 0 0],'fontweight','normal','fontsize',8);
if Iplot==1
%    title(sprintf('%s; mask on? %i  Active seismic only? %i Locations between %s and %s',Icase, mask_locations.on,run_options.shooting_only, ...
%        datestr(tabs_range(Ihr),0),datestr(tabs_range(Ihr+1),0)));
    % Turn off LaTeX interpreter for correct title syntax.  Make two-line
    % title for legibility. [khkim, 2 Dec 09]
    title1 = sprintf('Locations between %s and %s', ...
        datestr(tabs_range(Ihr),0),datestr(tabs_range(Ihr+1),0));
    title2 = sprintf('%s; Mask on? %i Active seismic only? %i', ...
        Icase, mask_locations.on,run_options.shooting_only);
    titleout={title1;title2};
    title({title1;title2},'Interpreter','None');

else
    title(sprintf('Manual results, locations restricted to %i km of centerlines',run_options.center_dist_limit/1000));
end
