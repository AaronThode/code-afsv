%crude_decimate_uneven.m%%%
%%% Take the mean, median, max, or min of a set of points collected at
%%% uneven time intervals.
%function [xout,tout,Iout,Ncount]=crude_decimate_uneven(x,t,sampling_interval,op_str,trange,averaging_interval, percentile);
% Input:
% x vector of values, assume asociated with
% t: time vector associated with x, same length as x.  Must be
%   monotonically increasing with time..
% sampling_interval: interval between output samples.
% op_str: 'max','min','median','mean','std','max_excursion','percentiles','XXpercent'
%       The 'max_excursion' operation searches for the global maximum and minimum in
%       the time window, then outputs (max(x)-min(x))/(Time separation
%       between the two)
%       'XXpercent' will return the XXth percentile of the data...
%       'percentiles' will return a matrix where the YYth row is the YYth percentile (10% to 90%, in 10%
%       increments)
% trange: [tstart tend] first and end element of tout desired.  
%       Otherwise, t(1) and t(end) will determine bounds of tout.
% averaging_interval: Optional argument: how much time to process per output sample.
%       ...must be same time units as sampling_interval and t, but must be larger than sampling_interval.
%           Default is sampling_interval.
% percentile:  if exists, defines a vector of percentiles to return with
%           op_str='percentiles'.  Values need to be between 0 and 1.
% Output:
% xout,tout: decimated x and t vectors
% Iout: indicies associated with tout in t.  Useful for applying to
% correlated time or other data vectors.
% Ncount:  If exists, export the number of samples used to compute the
%   output value at each time
function [xout,tout,Iout,Ncount]=crude_decimate_uneven(x,t,sampling_interval,op_str,trange,averaging_interval, percentile)

xout=[];
tout=[];
Iout=[];
if nargout==4
    Ncount=[];
end

if isempty(t)||isempty(x)
    disp('crude_decimate_uneven:  t or x is empty, cannot process');
    return;
end

start_flag=0;

if ~exist('averaging_interval','var')||isempty(averaging_interval)
    dT=sampling_interval;
else
    dT=averaging_interval;
end

if exist('percentile','var')&&~strcmp(op_str,'percentiles')
    disp('Warning! crude_decimate_uneven: op_str does not match percentile request');
end

if ~exist('percentile','var')||isempty(percentile)
   percentile=0.1:0.1:0.9; 
end

if exist('trange','var')==1&&~isempty(trange)
    start_flag=1;
    Igood=find(t>=trange(1)&t<=trange(2));
    x=x(Igood);
    t=t(Igood);
    tout_samp=trange(1):sampling_interval:trange(2);
    tspan=trange(2)-trange(1);
else
    tout_samp=min(t):sampling_interval:max(t);
    tspan=max(t)-min(t);
end


%disp(sprintf('Output time vector is %i samples',length(tout_samp)));
if dT>=tspan
    tout_samp=[tout_samp tout_samp(1)+tspan];
end


Nsegs=length(tout_samp)-1;


if Nsegs==0  %Time series is less than dN interval
    [xout,tout,Iout]=getFunction(op_str,x,t,1:length(x));
    if nargout==4
        Ncount=length(x);
    end
else  %process each segment, discarding last fraction...
   if ~strcmp(op_str,'percentiles')  %process each segment, discarding last fraction...
        xout=NaN*ones(Nsegs,1);
    else
        xout=NaN*ones(length(percentile),Nsegs);
   end
   tout=xout(1,:);
   Iout=xout;
   
   if nargout==4
       Ncount=tout;
   end
   for J=1:Nsegs
        %Is=(J-1)*dOvlap+(1:dN);
        Is=find((t>=tout_samp(J))&t<(tout_samp(J)+dT));
        if isempty(Is)
            continue
        end
        if ~strcmp(op_str,'percentiles')
            [xout(J),tout(J),Iout(J)]=getFunction(op_str,x(Is),t(Is),Is,[]);
        else
            [xout(:,J),tout(J),Iout(:,J)]=getFunction(op_str,x(Is),t(Is),Is,percentile);
        end
        if exist('trange','var')
            tout(J)=tout_samp(J);
        end
        if nargout==4
            Ncount(J)=length(Is);
        end
   end
   
   Igood=find(~isnan(tout));
    tout=tout(Igood);
  
    if ~strcmp(op_str,'percentiles')
        xout=xout(Igood);
        Iout=Iout(Igood);
    else
        
        xout=xout(:,Igood);
        Iout=Iout(:,Igood);
    end
end



%Igood = ~isnan(Iout);
%xout=xout(Igood);
%tout=tout(Igood);
%Iout=Iout(Igood);
%if nargout==4
%    Ncount=Ncount(Igood);
%end
end

function [xout,tout,Iout]=getFunction(op_str,x,t,Is,percentile)
switch op_str
    case 'max',
        [xout,Imin]=max(x);
        tout=t(Imin);
        Iout=Is(Imin);
    case 'min'
        [xout,Imin]=min(x);
        tout=t(Imin);
        Iout=Is(Imin);
    case 'median',
        xout=median(x);
        tout=mean(t);
        Iout=floor(mean(Is));
        
    case 'mean',
        xout=mean(x);
        tout=mean(t);
        Iout=floor(mean(Is));
        
    case 'std'
        xout=std(x);
        tout=mean(t);
        Iout=floor(mean(Is));
        
    case 'max_excursion'
        [xmax,Imax]=max(x);
        [xmin,Imin]=min(x);
        duration=abs(Imax-Imin);
        xout=(xmax-xmin)/duration;
        
        tout=mean(t);
        Iout=floor(mean(Is));
    case 'percentiles'
        [xsort,Isort]=sort(x);
        %percentile=0.1:0.1:0.9;
        Iout=round(percentile*length(xsort));
        xout=xsort(Iout)';
        Iout=Isort(Iout)';
        
        tout=mean(t);
        
    otherwise
        
        if ~isempty(findstr('percent',op_str))
            percentage=str2double(op_str(1:2))/100;
            xsort=sort(x);
            Iout=round(percentage*length(xsort));
            xout=xsort(Iout);
            tout=mean(t);
        else
            error('crude_decimate requested operation does not exist');
        end
end
end



