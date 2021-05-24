function [area,a,b,ang,angax] = ellipsparms(S,critval,MD,VM)
%ELLIPSPARMS  Calculates parameters of bivariate normal ellipse from covariance.
%   [AREA,A,B,ANG,angax] = ELLIPSPARMS(S,CRITVAL,MD,VM) where
% S is the covariance matrix,
% CRITVAL is the critical value from the chi-squared distribution,
% MD is the mean position of the DASAR array,
% VM is the estimated location.
% It returns
%   AREA of the ellipse,
%   A and B are lengths of major and minor full axes respectively, and
% ANG is angle in radians ccw of  major semi-axis.
%   The angle is referenced to an imaginary vector emanating from the center of the
%   DASAR array and directed toward the estimated location.
%   Note: index matrices using vec representation.
% ANGAX is angle in radians of major axis measured ccw from +X axis
%
% The adjustments to ANG (last set of if/then statements) reflect the
%   angle of the major axis 180 degrees about the ray directed from MD to
%   VM (if necessary) so that the major axis is always directed away from
%   MD.  The absolute orientation of the ellipse is unchanged; only its
%   "direction" relative to the center of the DASAR array may be adjusted.

dp = VM-MD;
theta = atan2(dp(2),dp(1));     % 4-quadrant inverse tangent: angle from MD to VM
detS = det(S);
if isnan(S(1)) || (detS<=0) || min(diag(S))<0,
    area = nan;  a = nan;  b = nan;  ang = nan;  angax = nan;
else
    area = pi*sqrt(detS)*critval;
    [V,L] = eig(S);
    d = 2*sqrt(critval*diag(L));  % Lengths of both axes.
    [d,i] = sort(d);
    a = d(2);  b = d(1);          % a & b: major & minor axis lengths, respectively.
    v = V(:,i(2));                % dominant eigenvector
    angax = atan(v(2)/v(1));      % Absolute orientation of ellipse major axis
    if angax<0,                   % Returns an angle between 0 and pi following
        angax = angax+pi;           %    math convention (0 radians = 90 degrees, East;
    end                           %    pi radians = 270 degrees, West).
    ang = angax;
    if ((theta>-pi/2) && (theta<pi/2)  && (ang>(pi/2+theta))), % Adjust angle for
        ang = ang-pi;                                          %   estimated position
    elseif ((theta>=pi/2) && (ang<(theta-pi/2))),             %   relative to center
        ang = ang+pi;                                          %   of DASAR array.
    elseif ((theta<=-pi/2) && (ang<(3*pi/2+theta))),
        ang = ang+pi;
    end
end
end
