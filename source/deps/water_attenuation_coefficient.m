%%%function [alphadB, alpha,alphadB_simple]=water_attenuation_coefficient(f,S,T,pH)
%%%f is in kilohertz
%%% alphadB is in units of dB/km
%%%  alpha is in nepers (i.e. the imaginary part of medium wavenumber k)

function [alphadB, alpha,alphadB_simple]=water_attenuation_coefficient(f,S,T,pH)

%%%%Attenuation coefficient

if nargin==1 %frequency only
    S=35;
    T=23;
    pH=8;
    D=0.1;  %km depth %%%??%%%
end
%R(1)=2.65-0.99;
%R(2)=3.52-0.99;

%AdB=10*log10(X099_av);
%AdB(AdB<Noise_floor)=NaN;

f1=0.78*sqrt(S/35).*exp(T/26);
f2=42*exp(T/17);
alpha_boric=0.106*exp((pH-8)/0.56).*(f1*f.^2)./(f1.^2+f.^2);
alpha_magnesium=0.52*(1+T/43)*(S/35).*exp(-D/6).*(f2.*f.^2)./(f2.^2+f.^2);
alphadB=alpha_boric+alpha_magnesium+0.00049*exp(-T/27-D/17).*f.^2;
alpha=alphadB/8686;

alphadB_simple=3.3e-3+0.11*(f.^2)./(1+f.^2)+(44*f.^2)./(4100+f.^2)+3e-4*f.^2;

if 1==0
    f=logspace(log10(0.01),log10(10000),500);
    [alphadB, alpha,alphadB_simple]=water_attenuation_coefficient(f);
    loglog(f*1000,alphadB,'k',f*1000,alphadB_simple,'r');grid on
end
%Ar265=AdB'-R(1).*alpha;
%Ar352=AdB'-R(2).*alpha;