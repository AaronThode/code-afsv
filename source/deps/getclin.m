
%%%%%%%%%%%%%%%%%%%%getc.m%%%%%%%%%%%%%%%%
% Create a linear sound speed profile
function c=getclin(z);
cmin=1500;zmin=1300;c0=1560;
g=(cmin-c0)/zmin;
c=c0+g*z;
deep=find(z>zmin);
c(deep)=cmin-g*(z(deep)-zmin);


