function [VM,Qhat,w,outcome] = vmmle_r(angle,dasar,r,k)
%[VM,QHAT,W,OUTCOME] = VMMLE(ANGLE,DASAR,R,K)
%   Von Mises MLE estimate of location of animal target, based on:
%
%   Lenth,R.V. 1981.  On finding the source of a signal.
%     Technometrics 23:149-154.
%
%   [VM,QHAT,W,OUTCOME] = VMMLE(ANGLE,DASAR,R,K) processes a single observation
%   comprised of compass ANGLE (an n-by-1 column vector of bearings), DASAR
%   (an n-by-2 matrix of coordinates for the DASAR source of each
%   bearing), and K (an n-by-1 vector of bearing standard error
%   estimates, expressed as the von Mises concentration parameter kappa,
%   associated with each DASAR).  R is a character code for method -
%   either 'm' for the standard (non-robust) MLE, 'a' for the Andrews
%   robust procedure, or 'h' for the Huber robust procedure.
%   VM is the von Mises estimate of location, QHAT is the associated
%   covariance matrix estimate and W is a 1-by-n row vector of weights:
%   for the MLE, either all 1's if the procedure converged or 0's if not;
%   for the Andrews, generally a mix of 1's and 0's; and for the Huber,
%   generally a mix of values between 0 and 1.  OUTCOME is a character
%   array identifying whether a good solution was found, and if not,
%   the reason for failure.

if nargin<3,
    r = 'm';
end
if nargin<4||all(k==0)          % Estimate kappa from the data
    kest = true;
else
    kest = false;
end
r = lower(r);
robust = strfind('mah',r);
if isempty(robust)
    error('Input r must be either "m", "a", or "h"');
end
failed = {'less than 2 bearings','negative variance estimates',...
    'solution behind DASARs','failed to converge'};
n = length(angle);

if n<=1  %If only one set of bearings present...
    outcome = failed{1};
else
    cond_num = 1e15;       % For test of singularity.
    tc = 1.5;              % Tuning constant for robust versions.
    dist1 = 1;             % Was 0.1
    maxiter = 50;
    x = dasar(:,1);
    y = dasar(:,2);
    iter = 0;
    theta = (90-angle)*pi/180;  %Mathematical angle defined
    theta = (theta<-pi)*2*pi + theta;
    s = sin(theta);
    c = cos(theta);
    z = s.*x - c.*y;
    sstar = s';
    cstar = c';
    w = ones(1,n);
    M1 = [sstar; -cstar]*[s -c];
    converge = 0;
    if cond2(M1)<cond_num,
        M2 = [sstar; -cstar]* z;
        xyhat = M1\M2;
        while (~converge)&&(iter<maxiter),
            iter = iter+1;
            xyold = xyhat;
            % d = dist(xyhat', [x y]');
            
            for JJ=1:length(x)
                d(JJ)=sqrt((xyhat(1)-x(JJ)).^2+(xyhat(2)-y(JJ)).^2);
            end
            if (robust>1) && (n>2), % Need 3 or more bearings to calculate weights for
                dxy = repmat(xyhat',n,1)-[x y];      % robust methods
                muhat = cart2pol(dxy(:,1),dxy(:,2));
                Cd = cos([theta-muhat]');
                if kest,
                    Cbar = abs(w*Cd'/sum(w))^(n/(n-2));  % Abs is ad hoc but may avoid temporary numeric problems in iteration
                    k = inv(2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar);
                end
                t = sqrt(2*k'.*(1-Cd));
                if robust==2,
                    phit = tc*sin(t/tc).*(abs(t)<tc*pi);
                else
                    phit = sign(t).*min(abs(t),tc);
                end
                same = (Cd==1);     % Take care when estimated & observed bearings
                t(same) = 1;        %   are identical; avoid division by zero.
                phit(same) = 1;     %   Usually occurs with just 2 intersecting bearings.
                w = phit./t;
            end    %  if robust>1
            sstar = w.*(xyhat(2)-y')./(d.^3);
            cstar = w.*(xyhat(1)-x')./(d.^3);
            M1 = [sstar; (-cstar)]*[s -c];
            if ((n-sum(~w))>1) && (cond2(M1)<cond_num),
                M2 = [sstar; -cstar]*z;
                xyhat = M1\M2;
                converge = sum(abs(xyhat-xyold)<dist1)==2;
            else
                break    % If either condition above occurs, convergence will very
            end        %    likely fail, so break out of while loop.
        end  % while
    end    % if cond2
    if converge,
        dxy = repmat(xyhat',n,1)-[x y];
        muhat = cart2pol(dxy(:,1),dxy(:,2));
        Cd = cos([theta-muhat]');
        if kest && (n>2),
            Cbar = (w*Cd'/sum(w))^(n/(n-2));   % Exponent is small sample size correction
            k = inv(2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar);
        end
        if Cd*w'>0,             % Weighted average of cosine differences (check
            VM = xyhat';          %   on bad solutions behind DASARs)
            if kest && (n==2),     % Cannot estimate Qhat with only 2 bearings
                Qhat = nan*ones(2);
                outcome = 'successful; 2 bearings; no Kappa';
            else
                k = k';
                cv = -(k.*sstar*c + k.*cstar*s)/2;
                M3 = [k.*sstar*s cv; cv k.*cstar*c];
                Qhat = inv(M3);
            end
            if ~kest || (n>2),
                if all(diag(Qhat)>0),
                    outcome = 'successful';       % Successful solution
                else
                    outcome = failed{2};        % Implausible variance estimate(s)
                end
            end
        else
            outcome = failed{3};          % Bad solution behind DASARs
        end % if all(Cd>0)
    else
        outcome = failed{4};            % No convergence
    end   % if converge
end     % if n<=1


if isempty(strfind(outcome,'successful'));
    %if ~isempty(strfind(failed,outcome))
    VM = [nan nan];
    Qhat = nan*ones(2);
    w = zeros(1,n);
end
end %vmmle_r
