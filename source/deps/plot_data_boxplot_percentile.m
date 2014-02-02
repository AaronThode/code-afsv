  function hprint=plot_data_boxplot_percentile(xdata,x_inc_boxplot,ydata,ylimits,xlabel_inc,label_style)
        
        %xdata: vector of data associated with x-axis.  Can be datenumbers
        %ydata: data associated with x-axis data xdata
        %  x_inc_boxplot: interval to compute boxplot over, same units as
        %       xdata
        % ylimits: two-element vector associated with acceptable values of ydata
        % xlabel_inc:  x value to skip between ticks.
        % label_style:  if exist

        if ~exist('label_style')
            label_style=[];
        end
        hprint=figure;
        t1=xdata(1);
        t2=xdata(end);
        tbin=t1:x_inc_boxplot:t2;  %Shorter time scale
        
        [Ncount,Binn]=histc(xdata,tbin);  %%OK, sort times of observations by category
        data_fin=NaN*ones(max(Ncount),length(Ncount));
        
        for Ibin=1:length(Ncount)
            data=ydata(find(Binn==Ibin));
            data(find(data<ylimits(1)|data>ylimits(2)))=NaN;
            data_fin(1:length(data),Ibin)=data';
            
        end
        make_boxplot(data_fin,tbin,x_inc_boxplot,xlabel_inc,label_style);
        ylim(ylimits);
        
      function make_boxplot(data_fin,tbin,x_inc_boxplot,xlabel_inc,label_style)
            boxPlot_percentile(data_fin)
            set(gca,'fontweight','bold','fontsize',14);
            
            Ibin=get(gca,'xtick');
            Istep=floor(xlabel_inc/x_inc_boxplot);
            Ibin=1:Istep:size(data_fin,2);
            set(gca,'xtick',Ibin);
            if ~isempty(label_style)
                set(gca,'xticklabel',datestr(tbin(Ibin),label_style));
            else
                set(gca,'xticklabel',round(tbin(Ibin)));
            end
            %ylim([80 130]);
            %ylabel('dB re 1uPa^2/Hz');
            %title(sprintf('%s, %s',names{7},dirs{Idir}));
            grid on
            %xlabel('Local time, starting on Aug. 18, 2011')
            
      end
    
  
    end

