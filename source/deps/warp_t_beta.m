function [ t_w ] = warp_t_beta(t,r,c,beta)

%t_w=t.^(1+beta);


% t0=abs(beta).*(t_w-r/c);
% t=(t0.^(1/(1+beta)));

t0=(t).^(1+(beta));
t_w=r/c+t0/(beta);