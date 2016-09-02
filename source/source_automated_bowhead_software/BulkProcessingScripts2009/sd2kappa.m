function k = sd2kappa(s)
%SD2KAPPA Converts bearing standard deviation (degrees) to kappa.
%  K = SD2KAPPA(S) converts a vector of bearing standard deviations to
%  a corresponding vector of kappa values (kappa is the concentration
%  parameter of the Von Mises distribution).  The formula for A (first
%  line of code) is from White and Garrott (1990), equation 4.3,
%  page 59.  The formula for K is originally from Lenth (1981),
%  equation 2.10, page 150.

A = exp(-0.5*(s*pi/180).^2);
k = 1./(2*(1-A)+(((1-A).^2).*(0.48794-0.82905*A-1.3915*A.^2))./A);