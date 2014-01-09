  function plot_data_boxplot(tabs,time_inc_boxplot,ydata,ylimits,tlabel_inc,label_style)
        
        %tabs vector of datenumbers associated with data
        %ydata, data associated with time series tabs
        %  time_inc_boxplot: datenumber of time interval to compute boxplot over
        % ylimits: two-element vector associated with acceptable values of ydata
        % tlabel_inc:  Number of ticks to skip when labeling time (x) axis
%         t1=datevec(tabs(1));
%         t1(4:6)=0;
%         t1=datenum(t1);  %Floors date to midnight previous day
%         
%         t2=datevec(tabs(end));
%         t2(4:6)=0;
%         t2(3)=t2(3)+1;
%         t2=datenum(t2);  %Ceilings date to midnight next day
%         
        figure;
        t1=tabs(1);
        t2=tabs(end);
        tbin=t1:time_inc_boxplot:t2;  %Shorter time scale
        
        [Ncount,Binn]=histc(tabs,tbin);  %%OK, sort times of observations by category
        data_fin=NaN*ones(max(Ncount),length(Ncount));
        
        for Ibin=1:length(Ncount)
            data=ydata(find(Binn==Ibin));
            data(find(data<ylimits(1)|data>ylimits(2)))=NaN;
            data_fin(1:length(data),Ibin)=data';
            
        end
        make_boxplot(data_fin,tbin,time_inc_boxplot,tlabel_inc,label_style);
      
      function make_boxplot(data_fin,tbin,time_inc_boxplot,tlabel_inc,label_style)
            boxPlot(data_fin)
            set(gca,'fontweight','bold','fontsize',14);
            
            Ibin=get(gca,'xtick');
            Istep=floor(tlabel_inc/time_inc_boxplot);
            Ibin=1:Istep:size(data_fin,2);
            set(gca,'xtick',Ibin);
            set(gca,'xticklabel',datestr(tbin(Ibin),label_style));
            %ylim([80 130]);
            %ylabel('dB re 1uPa^2/Hz');
            %title(sprintf('%s, %s',names{7},dirs{Idir}));
            grid on
            %xlabel('Local time, starting on Aug. 18, 2011')
            
      end
    
  
    end

