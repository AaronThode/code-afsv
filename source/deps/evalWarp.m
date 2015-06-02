function [val] = evalWarp(RTF,x_rect,y_rect,Fe_w,M)

ym=min(y_rect);
yM=max(y_rect);
f=RTF(ym:yM,min(x_rect):max(x_rect));

Nt=length(f(1,:));
curve=zeros(1,Nt);

for ii=1:Nt
    [~,curve(ii)]=max(f(:,ii));
end

interp=polyval(polyfit(1:Nt,curve,1),1:Nt);
plot((min(x_rect)+(0:Nt-1))/Fe_w,(interp+ym)*Fe_w/M,'k')
val=[norm(f,2) yM-ym interp(2)-interp(1)];

end

