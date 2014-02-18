function [output,hprint]=create_percentile_distributions(xdata,ydata,param, suppress_output)

%function hprint=percentile_boxplots(xdata,x_inc_boxplot,ydata,ylimits,xlabel_inc,percentiles,label_style)

% xdata: vector of data associated with x-axis.  Can be datenumbers
% ydata: data associated with x-axis data
%  param:  structure of variables including
%     x_inc: interval to compute boxplot over, same units as
%       xdata
%     y_limits: two-element vector associated with acceptable values of ydata
%     xlabel_inc:  x value to skip between ticks.
%     label_style:  if exists, sets datenumber style on x-axis
%     x_label, y_label: plot labels for x and y axis.
%     percentiles: percentile distribution...
%
% suppress_output:  If exists, do not plot data, just return
%      percentile data instead...
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

if ~exist('suppress_output')
    suppress_output=0;
end


t1=xdata(1);
t2=xdata(end);

%%Adjust so that time intervals start and end at midnight...
tspan=datevec(t2-t1);
t1_mid=datevec(t1);t1_mid(4:6)=0;t1_mid=datenum(t1_mid);
t2_mid=datevec(t2);t2_mid(4:6)=0;t2_mid(3)=t2_mid(3)+1;t2_mid=datenum(t2_mid);

%If the new span is much greater than old span, shrink the plot bounds.
frac=(t2-t1)./(t2_mid-t1_mid);
    
if frac<0.75 %This only occurs for less than a day...round to nearest hour
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



tbin=t1_mid:x_inc_boxplot:t2_mid;  %Shorter time scale

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
if suppress_output==0
    hprint=figure;

    make_boxplot(data_fin,tbin,x_inc_boxplot,xlabel_inc,percentiles,label_style,ylimits);
    if isfield(param,'x_label')
        xlabel(param.x_label,'fontsize',14,'fontweight','bold');
    end
    if isfield(param,'y_label')
        ylabel(param.y_label,'fontsize',14,'fontweight','bold');
    end

    if isfield(param,'title')
       title(param.title,'fontsize',14,'fontweight','bold'); 
    end
    
end

data=boxPlot_percentile(data_fin,percentiles,1);
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

    function make_boxplot(data_fin,tbin,x_inc_boxplot,xlabel_inc,percentiles,label_style,ylimits)
        test=boxPlot_percentile(data_fin,percentiles,0);
        
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

