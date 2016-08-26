function [ t ] = iwarp_t_beta(t_w,r,c,beta)

t0=(beta).*(t_w-r/c);
t=(t0.^(1/(1+beta)));