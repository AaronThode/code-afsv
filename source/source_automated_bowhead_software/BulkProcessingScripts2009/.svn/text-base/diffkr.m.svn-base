function d = diffkr(k,r,n)
%DIFFKR Difference of a function of Von Mises kappa and R, mean vector length.
%  D = DIFFKR(K,R) is used by Matlab function FMINBND to find the maximum
%  likelihood estimate of kappa, the concentration parameter of the Von Mises
%  distribution. The MLE of kappa is the value K such that |R - A(K)| is
%  minimum, where R is the length of the mean vector, A(K) =  I1(K)/I0(K), and
%  I1 and I0 are modified bessel functions of order 1 and 0, respectively.
%
%  A bias correction for small sample size (see Batschelet, 1981, p. 47) is
%  included so that K is sought for min{|R - A(K)/A(KRN)|} , where N is the
%  sample size.
%
%  Rather than call the Matlab function, BESSELI, faster polynomial
%  approximations from Abramowitz and Stegun (1965, p. 378) are used.  Also,
%  the approximations allow arguments > 700 (which cause overflow in BESSELI).
%  In the code below, I0 and I1 represent functions of I0(X) and I1(X),
%  respectively, depending on the value of the argument X.  For X>3.75, the
%  leading factor,  cancels in the ratio allowing calculation of A(X) without
%  numeric overflow.

%  Abramowitz, M. and I.A. Stegun. 1965.  Handbook of Mathematical Functions.
%  Dover Publications, New York.

%  Batschelet, E. 1981. Circular Statistics in Biology. Academic Press, London.

krn = k*r*n;
t = krn/3.75;
if krn<=3.75,
  I0 = 1 + 3.5156229*t^2 + 3.0899424*t^4 + 1.2067492*t^6 + 0.2659732*t^8 + ...
       0.0360768*t^10 + 0.0045813*t^12;
  I1 = krn * (0.5 + 0.87890594*t^2 + 0.51498869*t^4 + 0.15084934*t^6 + ...
       0.02658733*t^8 + 0.00301532*t^10 + 0.00032411*t^12);
  Akrn = I1/I0;
else
  I0 = 0.39894228 + 0.01328592/t + 0.00225391/t^2 - 0.00157565/t^3 + ...
       0.00916281/t^4 - 0.02057706/t^5 + 0.02635537/t^6 - 0.01647633/t^7 + ...
       0.00392377/t^8;
  I1 = 0.39894228 - 0.03988024/t - 0.00362018/t^2 + 0.00163801/t^3 - ...
       0.01031555/t^4 + 0.02282967/t^5 - 0.02895312/t^6 + 0.01787654/t^7 - ...
       0.00420059/t^8;
  Akrn = I1/I0;
end
t = k/3.75;
if k<=3.75,
  I0 = 1 + 3.5156229*t^2 + 3.0899424*t^4 + 1.2067492*t^6 + 0.2659732*t^8 + ...
       0.0360768*t^10 + 0.0045813*t^12;
  I1 = k * (0.5 + 0.87890594*t^2 + 0.51498869*t^4 + 0.15084934*t^6 + ...
       0.02658733*t^8 + 0.00301532*t^10 + 0.00032411*t^12);
  A = I1/I0;
else
  I0 = 0.39894228 + 0.01328592/t + 0.00225391/t^2 - 0.00157565/t^3 + ...
       0.00916281/t^4 - 0.02057706/t^5 + 0.02635537/t^6 - 0.01647633/t^7 + ...
       0.00392377/t^8;
  I1 = 0.39894228 - 0.03988024/t - 0.00362018/t^2 + 0.00163801/t^3 - ...
       0.01031555/t^4 + 0.02282967/t^5 - 0.02895312/t^6 + 0.01787654/t^7 - ...
       0.00420059/t^8;
  A = I1/I0;
end
d = abs(r - A./Akrn);
