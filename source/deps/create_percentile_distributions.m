function [output,hprint]=create_percentile_distributions(xdata,ydata,param)

%function hprint=percentile_boxplots(xdata,x_inc_boxplot,ydata,ylimits,xlabel_inc,percentiles,label_style)

% xdata: vector of data associated with x-axis.  Can be datenumbers
% ydata: data associated with x-axis data, rows are for separate plots
%  param:  structure of variables including
%     x_inc: interval to compute boxplot over, same units as
%       xdata
%     y_limits: two-element vector associated with acceptable values of ydata
%     xlabel_inc:  x value to skip between ticks.
%     label_style:  if exists, sets datenumber style on x-axis
%     x_label, y_label: plot labels for x and y axis.
%     percentiles: percentile distribution...
%     title: title to plot, if plot generated
%     plot:  structure with boolean fields
%       'none','boxplot','line','image'
%
% Output:
%   output: fields  x,  percentiles, data [Npercentiles, Nx];
%   hprint:  handle to figure making output

x_inc_boxplot=param.x_inc;
ylimits=param.y_limits;
xlabel_inc=param.xlabel_inc;
if ~isfield(param,'percentiles')
    percentiles=[0.01 0.1 0.25 0.5 0.75 0.9 0.99];
else
    percentiles=param.percentiles;
end

if isfield(param,'label_style')
    label_style=param.label_style;
else
    label_style=[];
end


hprint=[];
output=[];


t1=xdata(1);
t2=xdata(end);

%%Adjust so that time intervals start and end at midnight...
tspan=datevec(t2-t1);
t1_mid=datevec(t1);t1_mid(4:6)=0;t1_mid=datenum(t1_mid);
t2_mid=datevec(t2);t2_mid(4:6)=0;t2_mid(3)=t2_mid(3)+1;t2_mid=datenum(t2_mid);

%If the new span is much greater than old span, shrink the plot bounds.
frac=(t2-t1)./(t2_mid-t1_mid);
    
if frac<0.75 %The events in question only occur for less than a day...round to nearest hour
   t1_mid=datevec(t1);t1_mid(5:6)=0;t1_mid=datenum(t1_mid);
   t2_mid=datevec(t2);t2_mid(5:6)=0;t2_mid(4)=t2_mid(4)+1;t2_mid=datenum(t2_mid);
   frac=(t2-t1)./(t2_mid-t1_mid);    
end


%Placeholder for percentile measures that cover much less than an hour
if frac<0.75 %Data covers less than an hour, round to nearest minute.
   t1_mid=datevec(t1);t1_mid(6)=0;t1_mid=datenum(t1_mid);
   t2_mid=datevec(t2);t2_mid(6)=0;t2_mid(5)=t2_mid(5)+1;t2_mid=datenum(t2_mid);
   frac=(t2-t1)./(t2_mid-t1_mid);    
end



tbin=t1_mid:x_inc_boxplot:t2_mid;  %This is the 'x' output: times for which percentiles computed.

[Ncount,Binn]=histc(xdata,tbin);  %%OK, sort times of observations by category
data_fin=NaN*ones(max(Ncount),length(Ncount));


for Ibin=1:length(Ncount)
    data=ydata(find(Binn==Ibin));
    data(find(data<ylimits(1)|data>ylimits(2)))=NaN;
    data_fin(1:length(data),Ibin)=data';
    
end

if isempty(data_fin)
   uiwait(msgbox('data_fin is empty: you have requested odd times for percentiles')); 
end

if param.plot.boxplot
    hprint=figure;

    data=make_boxplot(data_fin,tbin,x_inc_boxplot,xlabel_inc,percentiles,label_style,ylimits);
    if isfield(param,'x_label')
        xlabel(param.x_label,'fontsize',14,'fontweight','bold');
    end
    if isfield(param,'y_label')
        ylabel(param.y_label,'fontsize',14,'fontweight','bold');
    end

    if isfield(param,'title')
       title(param.title,'fontsize',14,'fontweight','bold'); 
    end
    
    %xlabel(pms.x_label,'fontsize',14,'fontweight','bold');
    %ylabel(pms.y_label,'fontsize',14,'fontweight','bold');
    
    
    
end

if param.plot.line
    hprint=figure;

    data=boxplot_percentile(data_fin,percentiles,1);
    Igood=~any(isnan(data));
    
    plot(tbin(Igood),data(:,Igood)');
    datetick('x',param.label_style);
    ylim(param.y_limits);
    %xlim([min(tbin) max(tbin)]);
    grid on;
    clear legstr
    for II=1:length(param.percentiles)
        legstr{II}=num2str(param.percentiles(II),3);
    end
    legend(legstr)
    title(param.title);
    xlabel('Time','fontweight','bold','fontsize',14);
    ylabel('Integrated power (dB re 1uPa)','fontweight','bold','fontsize',14);
    
    
end

if param.plot.none||param.plot.image
    hprint=-1;
    data=boxplot_percentile(data_fin,percentiles,1);
end

if isempty(data)
    output.data=[];
    output.x=[];
    output.percentiles=[];
    return
end
output.data=data;
output.x=tbin;
output.percentiles=percentiles;
%ylim(ylimits);

    function test=make_boxplot(data_fin,tbin,x_inc_boxplot,xlabel_inc,percentiles,label_style,ylimits)
        test=boxplot_percentile(data_fin,percentiles,0);
        
        if isempty(test)
            return
        end
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
        ylim(ylimits);
        %ylim([80 130]);
        %ylabel('dB re 1uPa^2/Hz');
        %title(sprintf('%s, %s',names{7},dirs{Idir}));
        grid on
        
        
        %xlabel('Local time, starting on Aug. 18, 2011')
        
    end


end

