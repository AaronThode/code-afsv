function [s_r, Fe]=iwarp_temp_exa(s_w,Fe_w,r,c,Fe,N)

% Inverse (time) warping function using isovelocity waveguide as a warping model.
% Inputs :
% s_w : (warped) signal that will be unwarped
% Fe_w : sampling frequency of s_w
% r : warping parameter (source/receiver distance)
% c : warping parameter (source/receiver distance)
% Fe : sampling frequency of the new unwarped signal s
% N : number of points of the new unwarped signal s
% Sorties :
% s : signal after inverse warping

%% Step 1: preliminary computations
M=length(s_w);

if iscolumn(s_w)
    s_w=s_w';
end

%% Step 5: warped signal computation

% Time axis, uniform sampling (starts from r/c)
t=(1:N)/Fe+r/c;

% Time axis, non-uniform sampling
t_iw=iwarp_t(t,r,c);

% factor for energy conservation
coeff=sqrt(Fe_w/Fe)*sqrt(t./t_iw); % Energy conservation     

% Start exact interpolation (Shannon)
s_aux=repmat(s_w,N,1);  % initial signal replicated N times (N rows)
aux1=repmat(t_iw',1,M);
aux2=repmat((0:M-1)/Fe_w,N,1);
aux=sinc(Fe_w*(aux1-aux2));
% end of exact interpolation --> interpolated signal is sum(s_aux.*aux,2)

% Final warped signal
s_r=real(coeff'.*sum(s_aux.*aux,2));
      