function [ t_w ] = warp_t(t,r,c)

t_w=sqrt(t.^2+r^2/c^2);

end