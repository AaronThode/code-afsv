function [ t ] = iwarp_t(t_w,r,c)

t=sqrt(t_w.^2-r^2/c^2);

%t=(beta*(t-r/c)).^(1/(1+beta));
end
