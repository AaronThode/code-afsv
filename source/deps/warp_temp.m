function [s_w, Fe_w]=warp_temp(s,Fe,r,c)

% Time warping function using isovelocity waveguide as a warping model.
% Inputs:
% s : signal that will be warped
% Fe : sampling frequency of s
% r : warping parameter (source/receiver range)
% c : warping parameter (water sound speed)
% Outputs :
% s_w : warped signal
% Fe_w : sampling frequency of s_w (i.e. sampling frequency
%        of the warped-time domain)

% NB : warping actually has a single parameter : r/c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%% Important %%%%%%%%%%%%%%%%%%%%%%%%
% Signal s should be defined such that its first sample
% corresponds to the first modal arrival. In an isovelocity or Pekeris
% waveguide, if an impulsive source explodes at time ts, then first sample
% of s should be take at time ts+r/c. Warping is way more sensitive to
% that setting than to the actual parameters r and c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Step 1: preliminary computations
s=real(s);
N=length(s);
dt=1/Fe;

tmax=(N-1)/Fe+r/c;
tmin=r/c;

%% Step 2: new time step
dt_w=iwarp_t(tmax,r,c)-iwarp_t(tmax-dt,r,c);

%% Step 3: new sampling frequency
Fe_w=2/dt_w;

%% Step 4: new number of points
t_w_max=iwarp_t(tmax,r,c);
t_w_min=0;
M= ceil(t_w_max*Fe_w);


%% Step 5: warped signal computation

s_w=zeros(M,1);

for m=1:M
    
    t_w=m/Fe_w;
    % energy conservation
    toto=sqrt(Fe/Fe_w)*sqrt(t_w/sqrt(t_w^2+r^2/c^2));
    
    % linear interpolation of signal from original time domain
    n=(warp_t(t_w,r,c)-r/c)*Fe;
    n1=floor(n);
    n2=n1+1;
    dn=n-n1;
    
    if n2<=N && n1>0
        a=(s(n1)-s(n2))/(n1-n2);
        val=s(n1)+a*dn;
        s_w(m)=toto*val;
    end
    
end


end
