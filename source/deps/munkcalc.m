%%%%%%%%%%%%%%%%%%%%%munkcalc.m%%%%%%%%%%%%%%%%%%%
%Given a column vector x containing [r, z, rho, eta], return
%  the quantity dx/ds for the munk profile.

function xprime=munkcalc(~,x)

r=1;   %x(1) is r
z=2;   %x(2) is z
rho=3; %x(3) is rho
eta=4; %x(4) is eta
tau=5;

c=getcmunk(x(z));
dcdr=0;
zbar=2*(x(z)-1300)/1300;etap=0.00737;
dcdz=3000*etap*(1-exp(-zbar))/1300;

xprime(r,1)=c*x(rho);
xprime(z,1)=c*x(eta);
xprime(rho,1)=-dcdr/c.^2;
xprime(eta,1)=-dcdz/c.^2;
xprime(tau,1)=1./c;
%disp(xprime)
