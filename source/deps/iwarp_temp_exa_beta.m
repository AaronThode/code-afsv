function [s_r]=iwarp_temp_exa_beta(s_w,t_w,Fe_w,r,c,beta,Fe)
%%%%t_w is evenly-sampled warped time
%% Step 1: preliminary computations
M=length(s_w);

if iscolumn(s_w)
    s_w=s_w';
end

%% Step 5: warped signal computation


t_real_uneven=warp_t_beta(t_w,r,c,beta);
N=floor((max(t_real_uneven)-min(t_real_uneven))*Fe);
t_real_even=linspace(min(t_real_uneven),max(t_real_uneven),N);
% Time axis, uniform sampling (starts from r/c)
%t=tmin+(0:N-1)/Fe;

% Warped axis, non-uniform sampling
t_warp_uneven=iwarp_t_beta(t_real_even,r,c,beta);

% factor for energy conservation
coeff=sqrt(Fe_w/Fe)*sqrt(t_real_even./t_warp_uneven); % Energy conservation     .

% Start exact interpolation (Shannon)
s_aux=repmat(s_w,N,1);  % initial signal replicated N times (N rows)
aux1=repmat(t_warp_uneven',1,M);
aux2=repmat(t_w,N,1);
aux=sinc(Fe_w*(aux1-aux2));
% end of exact interpolation --> interpolated signal is sum(s_aux.*aux,2)

% Final warped signal
s_r=real(coeff'.*sum(s_aux.*aux,2));

 if iscolumn(s_r)
     s_r=s_r';
 end
      