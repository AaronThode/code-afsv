%%%%%%%%%%%%%%%%%%%%%generalcalc.m%%%%%%%%%%%%%%%%%%%
%Given a column vector x containing [r, z, rho, eta], return
%  the quantity dx/ds for a linear profile.

function xprime=generalcalc(s,x)



dcdr=0;

%[c,dcdz]=getcgeneral(x(2));        % Uses a polynomial fit for rapid
%       calculations
[c,dcdz]=getcgeneral_crude(x(2));  %Interpolates from input profile

xprime(1,1)=c*x(3);
xprime(2,1)=c*x(4);
xprime(3,1)=-dcdr/c.^2;
xprime(4,1)=-dcdz/c.^2;
xprime(5,1)=1./c;
%disp([xprime x])
