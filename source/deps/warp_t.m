function [ t_w ] = warp_t(t,r,c)

t_w=sqrt(t.^2+r^2/c^2);

%t_w=r/c + (t.^(1+beta))/beta;

end