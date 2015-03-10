
function [s_r]=iwarp_temp(s_w,Fe_w,r,c,Fe,N)

% Inverse (time) warping function using isovelocity waveguide as a warping model.
% Inputs :
% s_w : (warped) signal that will be unwarped
% Fe_w : sampling frequency of s_w
% r : warping parameter (source/receiver distance)
% c : warping parameter (source/receiver distance)
% Fe : sampling frequency of the new unwarped signal s
% N : number of points of the new unwarped signal s
% Sorties :
% s : signal after inverse warping


%%

M=length(s_w);
s_r=zeros(N,1);

for n=1:N
    
    t=n/Fe+r/c;
    % Energy conservation
    toto=sqrt(Fe_w/Fe)*sqrt(t/sqrt(t^2-r^2/c^2));
    
    % Interpolation
    m=(iwarp_t(t,r,c))*Fe_w;
    m1=floor(m);
    m2=m1+1;
    dm=m-m1;
    
    if m2<=M && m1>0
        a=(s_w(m1)-s_w(m2))/(m1-m2);
        val=s_w(m1)+a*dm;
        s_r(n)=toto*val;
    end
    
end


end

