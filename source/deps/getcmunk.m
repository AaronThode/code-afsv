
%%%%%%%%%%%%%%%%%%%%getcmunk.m%%%%%%%%%%%%%%%%
% Create a munk sound speed profile
function c=getcmunk(z);

eta=0.00737;
zbar=2*(z-1300)/1300;
c=1500*(1+eta*(zbar-1+exp(-zbar)));


