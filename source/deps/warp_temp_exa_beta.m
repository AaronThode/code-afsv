%function [s_w, Fe_w]=warp_temp_exa_beta(s,Fe,r,c,beta)
function [s_w, Fe_w]=warp_temp_exa_beta(s,Fe,r,c,beta)

if iscolumn(s)
    s=s.';
end

%% Step 1: preliminary computations
s=real(s);
N=length(s);
dt=1/Fe;

%tmin=0;
%tmax=(N-1)/Fe;

tmin=r/c;
tmax=(N-1)/Fe+r/c;

%% Step 2: new time step
dt_w=iwarp_t_beta(tmax,r,c,beta)-iwarp_t_beta(tmax-dt,r,c,beta);

%% Step 3: new sampling frequency
Fe_w=2/dt_w;

%% Step 4: new number of points
t_w_max=iwarp_t_beta(tmax,r,c,beta);
M=ceil(t_w_max*Fe_w);


%% Step 5: warped signal computation

% Warped time axis, uniform sampling
t_w=(0:M-1)/Fe_w;

% Warped time axis, non-uniform sampling (starts from r/c)
t_ww=warp_t_beta(t_w,r,c,beta);

% factor for energy conservation
coeff=sqrt(Fe/Fe_w)*sqrt(t_w./t_ww);     

% Start exact interpolation (Shannon)
s_aux=repmat(s,M,1);
aux1=repmat(t_ww',1,N); % size=(M,1) -> repmat(1,N)
aux2=repmat(tmin+(0:N-1)/Fe,M,1);       % size=(1,N) -> repmat(M,1)
aux=sinc(Fe*(aux1-aux2));           % size=(M,N)
% end of exact interpolation --> interpolated signal is sum(s_aux.*aux,2)

% Final warped signal
s_w=coeff'.*sum(s_aux.*aux,2);