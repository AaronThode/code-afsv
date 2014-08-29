function boxplot(data0, lineWidth, width)
% boxPlot(data0) - plot box-whiskers diagram, accept multiple columns
%  shows max, 25,50,75, max percentiles
% Arguments: data0 -  unsorted data, mxn, m samples, n columns
%            lineWidth -  line thickness in the plot default = 1;
%            width -  the width of the box, default = 1;
%            xlab- x-axis to be plotted
%% Fill in extra bins with NaN...
% Returns:
% Notes: each column is considered as a single set


if(nargin < 3)
    width = 1;
end;
if(nargin < 2)
    lineWidth = 1;
end;


[m n] = size(data0);

data = sort(data0, 1); % ascend


for i=1:n  %For each column...
    ni=find(isnan(data(:,i)));
    if length(ni)==m
        q2(i)=0;
        q1(i)=0;
        q3(i)=0;
        min_v(i)=0;
        max_v(i)=0;
        continue
    elseif isempty(ni)  %All data good...
        ni(1)=m+1;
        
    end
    mi=ni(1)-1;
    q2(i)=median(data(1:mi,i));
    if(rem(mi,2) == 0)
        upperA = data(1:mi/2,i);
        lowA = data(mi/2+1:mi,i);
    else
        upperA = data(1:round(mi/2), i);
        lowA = data(round(mi/2):mi, i);
    end;
    
    q1(i) = median(upperA, 1);
    q3(i) = median(lowA, 1);
    
    min_v(i) = data(1,i);
    max_v(i) = data(mi,i);
    
end
draw_data = [max_v; q3; q2; q1; min_v];

% adjust the width
drawBox(draw_data, lineWidth, width);


return;


function drawBox(draw_data, lineWidth, width)

n = size(draw_data, 2);

unit = (1-1/(1+n))/(1+9/(width+3));

hold on;
for i = 1:n
    
    v = draw_data(:,i);
    
    % draw the min line
    plot([i-unit, i+unit], [v(5), v(5)],'k', 'LineWidth', lineWidth);
    % draw the max line
    plot([i-unit, i+unit], [v(1), v(1)], 'k','LineWidth', lineWidth);
    % draw middle line
    plot([i-unit, i+unit], [v(3), v(3)], 'k','LineWidth', lineWidth);
    % draw vertical line
    plot([i, i], [v(5), v(4)], 'k','LineWidth', lineWidth);
    plot([i, i], [v(2), v(1)], 'k','LineWidth', lineWidth);
    % draw box
    plot([i-unit, i+unit, i+unit, i-unit, i-unit], [v(2), v(2), v(4), v(4), v(2)],'k', 'LineWidth', lineWidth);
    
end;
hold off
return;