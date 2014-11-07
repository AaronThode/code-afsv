%%%%%%%%%%%%%%%%%%%%%linearcalc.m%%%%%%%%%%%%%%%%%%%
%Given a column vector x containing [r, z, rho, eta], return
%  the quantity dx/ds for a linear profile.

function xprime=linearcalc(s,x)

r=1;   %x(1) is r
z=2;   %x(2) is z
rho=3; %x(3) is rho
eta=4; %x(4) is eta

dcdr=0;

zmin=1300;
cmin=1500;c0=1560;
g=(cmin-c0)/zmin;

c=getclin(x(z));
if (x(z)>zmin)
	dcdz=-g;
else	
	dcdz=g;
end

xprime(r,1)=c*x(rho);
xprime(z,1)=c*x(eta);
xprime(rho,1)=-dcdr/c.^2;
xprime(eta,1)=-dcdz/c.^2;

%disp(xprime)
