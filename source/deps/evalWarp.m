function [val] = evalWarp(RTF,x_rect,y_rect)

ym=min(y_rect);
yM=max(y_rect);
f=RTF(ym:yM,min(x_rect):max(x_rect));
val=[norm(f,2) yM-ym];

end

