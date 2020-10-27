function percentile = circ_percentile(alpha,dim)
%
% med = circ_median(alpha)
%   Computes the median direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians, from 0 to 2*pi
%     [dim  compute along this dimension, default is 1, must
%           be either 1 or 2 for circ_median]
%
%   Output:
%     mu		mean direction
%
%   circ_median can be slow for large datasets
%
% Update 2012
% PHB 3/19/2009
%
% References:
%   Biostatistical Analysis, J. H. Zar (26.6)
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
    dim = 1;
end

if nargin<2
    per=0.5; %median
end

M = size(alpha);
med = NaN(M(3-dim),1);
for I=1:M(3-dim)
    if dim == 2
        beta = alpha(I,:)';
    elseif dim ==1
        beta = alpha(:,I);
    else
        error('circ_median only works along first two dimensions')
    end
    
    beta = mod(beta,2*pi);  %%Restrict range of angles between 0 and 2*pi
    
    %beta=sort(beta);
    [md,beta2]=get_circ_median(beta);
    
    med(I) = md;
    percentile(4,I)=md; %50th
    
    [percentile(3,I),beta3m]=get_circ_median(beta2{2});  %25th
    [percentile(2,I),beta3m]=get_circ_median(beta3m{2}); %25th-12.5=12.5 percentile
    [percentile(1,I),beta3m]=get_circ_median(beta3m{2}); %12.5-6.25=6.25 percentile
    Icross=find(percentile(1:3,I)>md);
    percentile(Icross,I)=percentile(Icross,I)-2*pi;
    
    
    [percentile(5,I),beta3p]=get_circ_median(beta2{1}); %75th
    [percentile(6,I),beta3p]=get_circ_median(beta3p{1}); %75th+12.5=87.5 percentile
    [percentile(7,I),beta3p]=get_circ_median(beta3p{1}); %87.5+6.2500=93.75 percentile
    Icross=find(percentile(5:7,I)<md);
    II=5:7;
    percentile(II(Icross),I)=percentile(II(Icross),I)+2*pi;
    
  
    
    
end%I


if dim == 2
    med = med';
    percentile=percentile';
end

end

function [md,beta_out]=get_circ_median(beta)

debug=false;
n = size(beta,1);


dd = circ_dist2(beta,beta);  % range between -pi and pi


m1 = sum(dd>=0,1);  %%Angles greater than this angle
m2 = sum(dd<=0,1);

%dm = abs(m1-m2);  %%%Minimum vallue means m1==m2, or median


dm = abs(m1-m2);  %%%Minimum vallue means m1==m2, or median

m = min(dm);
idx = find(dm==m);


if length(idx) > 1  %%Multiple points satisfy median criteria
    %disp('Ties detected.') %#ok<WNTAG>
    
    %%%Keep only points that lie within 90 degrees of each other
    %%%% (get rid of angles that are 180 degrees away; deal with them
    %%%% later)...
    delta_beta=abs(circ_dist2(beta(idx),beta(idx(1))));
    idx=idx(delta_beta<pi/2);
    % if length(delta_beta)>2
    %     fprintf('Many possible medians with delta: %s\n',mat2str(delta_beta',3));
    % end
    
end

dd_med=mean(dd(:,idx),2);
beta_out{1}=beta(dd_med>0);
beta_out{2}=beta(dd_med<0);


md = circ_mean(beta(idx));  %Tentative median
if md<0,md=md+2*pi;end
%%If median is 180 degrees off from center of mass, add 180 degrees
circc_mean=circ_mean(beta);
if abs(circ_dist(circc_mean,md)) > abs(circ_dist(circc_mean,md+pi))
    md = mod(md+pi,2*pi);  %Restrict to 0 to 2*pi
    beta_out=beta_out(2:-1:1);
end

if debug
    figure
    subplot(3,1,1)
    histogram(beta,linspace(0,2*pi,90));
    
    subplot(3,1,2)
    histogram(beta_out{1},linspace(0,2*pi,90));
    subplot(3,1,3);
    histogram(beta_out{2},linspace(0,2*pi,90));
    hold on;
    plot(md,0,'ko');
    pause;close
end

end