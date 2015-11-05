%%%pekeris group and phase speed.m%%%
clear
close all
f=1:0.5:150;
c1=[1465 1465 1465];
c2=[1550 1600 1650];
rho1=1;
rho2=[1.8 1.8 1.8];
D=51;

ptchc='.ox+*sd';

for Icase=1:length(c1)
    
    vp=zeros(10,length(f));
    vg=zeros(10,length(f));
    for If=1:length(f)
        k1=2*pi*f(If)/c1(Icase);
        k2=2*pi*f(If)/c2(Icase);
        kr=linspace(k2,k1,1000);
        kz1=sqrt(k1.^2-kr.^2);
        kz2=sqrt(kr.^2-k2.^2);
        
        ff=sin(kz1.*D).*rho1.*kz2+rho2(Icase)*kz1.*cos(kz1.*D);
        plot(kr,ff);grid on
        ff=sign(ff);
        Igood=find(ff(2:end).*ff(1:(end-1))<1);
        
        if isempty(Igood)
            continue
        end
        kr=kr(Igood);
        vp(1:length(Igood),If)=2*pi*f(If)./fliplr(kr);
        
        
        term1=D.*((kr.^2-k2^2)).^1.5./(cos(D*sqrt(k1.^2-kr.^2)).^2);
        term2=(rho2(Icase)./rho1)*(k1.^2-k2.^2);
        Vn=kr.*(term1+term2);
        term3=(rho2(Icase)./rho1)*(kr.^2-k2.^2);
        term4=(rho2(Icase)./rho1).*k2*(k1.^2-kr.^2)./c2(Icase);
        Vn=Vn./((k1./c1(Icase)).*(term1+term3)+term4);
        vg(1:length(Igood),If)=fliplr(Vn);
        
    end
    
    plot(f,vp',ptchc(Icase),f,vg',ptchc(Icase));ylim([1200 1800]);grid on;hold on
    set(gca,'xtick',0:10:max(f))
    
end