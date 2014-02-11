function draw_data=boxplot_percentile(data0,percentiles,suppress_output, lineWidth, width)
% boxPlot(data0) - plot box-whiskers diagram, accept multiple columns
%  shows max, 25,50,75, max percentiles
% Arguments: data0 -  unsorted data, mxn, m samples, n columns.  If each
%               column has different sample sizes, fill with NaN.
%            lineWidth -  line thickness in the plot default = 1;
%            width -  the width of the box, default = 1;
%            xlab- x-axis to be plotted
%% Fill in extra bins with NaN...
% Returns:
% Notes: each column is considered as a single set


if(nargin < 5)
    width = 1;
end;
if(nargin < 4)
    lineWidth = 1;
end
if nargin<3
   suppress_output=0; 
end
if nargin<2
    percentiles=[0.01 0.1 0.25 0.5 0.75 0.9 0.99];
end

if rem(length(percentiles),2)==0
    uiwait(msgbox('Need an odd number of percentiles to plot!'));
    return
    
end
[m n] = size(data0);

data = sort(data0, 1); % ascend


for i=1:n  %For each column... (need to do because of different numbers of nan's)
    ni=find(isnan(data(:,i)));
    if length(ni)==m
        q2(i)=0;
        q1(i)=0;
        q3(i)=0;
        q_least(i)=0;
        q_most(i)=0;
        continue
    elseif isempty(ni)  %All data good...
        ni(1)=m+1;
        
    end
    mi=ni(1)-1;  %Index of final good values...
    
    Ip=floor(percentiles*mi);
    Ibad=find(Ip==0);
    Ip(Ibad)=1;
    
    Ip_frac=(percentiles*mi-Ip)';
    
    draw_data(:,i)=(1-Ip_frac).*data(Ip,i)+Ip_frac.*data(Ip+1,i); %Average two values in between desired percentile...
    draw_data(Ibad,i)=data(1,i);  %Use minimum value if percentile too small
    
    
    
end
%draw_data = [q_most; q3; q2; q1; q_least];

% adjust the width

if suppress_output==0
    drawBox_percentile(draw_data, lineWidth, width);
    title(sprintf('Percentiles: %s',num2str(percentiles)));grid on
end
return;


function drawBox_percentile(draw_data, lineWidth, width)

n = size(draw_data, 2);
np=size(draw_data,1);
if rem(np,2)==0
    msgbox('Need an odd number of percentiles to plot!');
    return
    
end

unit = (1-1/(1+n))/(1+9/(width+3));

box_width=linspace(0,unit,(np-1)/2);
box_width=[box_width fliplr(box_width)];

hold on;
for i = 1:n  %For each column
    
    v = draw_data(:,i);
    
    
    for J=1:length(box_width)
        
        % draw the horizontal line
        plot([i-box_width(J), i+box_width(J)], [v(J), v(J)],'k', 'LineWidth', lineWidth);
        
        % draw box
        plot([i-box_width(J), i+box_width(J), i+box_width(J), i-box_width(J), i-box_width(J)], [v(J), v(J), v(J+1), v(J+1), v(J)],'k', 'LineWidth', lineWidth);
    end
    plot([i-unit, i+unit], [v(J+1), v(J+1)],'k', 'LineWidth', lineWidth);
    plot([i-unit, i+unit], [v(1), v(1)],'k', 'LineWidth', lineWidth);
    
end;
hold off
return;


function drawBox(draw_data, lineWidth, width)