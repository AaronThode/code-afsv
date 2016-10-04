%function [s_w, Fe_w,t_w_min]=warp_temp_exa_beta(s,Fe,r,c,beta,dtmax,dt_rc)
% dtmax and dt_rc are times of tmax (maximum evaluation time) and r/c in seconds relative to s(1), r in meters, c in c/sec
%  beta is waveguide invariant.  If zero, use Bonnel's ideal waveguide
%  warping
function [s_w, Fe_w,t_w]=warp_temp_exa_beta(s,Fe,r,c,beta,dtmax,dt_rc)

%%% If beta<0, signal has been time-reversed

if iscolumn(s)
    s=s.';
end

%% Step 1: preliminary computations
s=real(s);
N=length(s);
dt=1/Fe;

%tmin=0;
%tmax=(N-1)/Fe;
%% Step 2: new time step

if beta>0 %dtmax not used
    tmin=r/c;
    tmax=(N-1)/Fe+r/c;
    dt_w=iwarp_t_beta(tmax,r,c,beta)-iwarp_t_beta(tmax-dt,r,c,beta);

else
    %tmax=r/c;
    if isempty(dtmax)
        error('warp_temp_exa_beta:  tmax is empty, although beta is < 0');
        return
    end
    tmax=r/c-(dt_rc-dtmax); %Input tmax is defined relative to start of signal
    tmin=r/c-dt_rc;  %tmin is the time at first sample in Ulysses window
    Npts=floor((tmax-tmin)*Fe);
    dt_w=abs(iwarp_t_beta(tmin,r,c,beta)-iwarp_t_beta(tmin+dt,r,c,beta));

end

%% Step 3: new sampling frequency
Fe_w=2/dt_w;

%% Step 4: new number of points
if beta>0
    t_w_max=iwarp_t_beta(tmax,r,c,beta);
    M=ceil(t_w_max*Fe_w);

else
    t_w_max=iwarp_t_beta(tmax,r,c,beta);  %Should be the same for all j_list
    t_w_min=iwarp_t_beta(tmin,r,c,beta);
    M=ceil((t_w_max-t_w_min)*Fe_w);

end


%% Step 5: warped signal computation

% Warped time axis, uniform sampling
if beta>0
    t_w=(1:M)/Fe_w;
else
    t_w=t_w_min+(1:M)/Fe_w;
    %s=fliplr(s);  %un time-reverse signal. 
end

%  non-uniform sampling of physical time (starts from r/c)
t_ww=warp_t_beta(t_w,r,c,beta);

% factor for energy conservation
coeff=sqrt(Fe/Fe_w)*sqrt(t_w./t_ww);     

% Start exact interpolation (Shannon)
% s_aux=repmat(s,Npts,1);
% aux1=repmat(t_ww',1,Npts); % size=(M,1) -> repmat(1,N)
% aux2=repmat(tmin+(0:Npts-1)/Fe,M,1);       % size=(1,N) -> repmat(M,1)

s_aux=repmat(s,M,1);
aux1=repmat(t_ww',1,N); % size=(M,1) -> repmat(1,N)
aux2=repmat(tmin+(0:N-1)/Fe,M,1); %physical times associated with s
            % size=(1,N) -> repmat(M,1)

aux=sinc(Fe*(aux1-aux2));           % size=(M,N)
% end of exact interpolation --> interpolated signal is sum(s_aux.*aux,2)

% Final warped signal
s_w=coeff'.*sum(s_aux.*aux,2);