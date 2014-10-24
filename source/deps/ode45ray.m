%%%%%%%%%%%%%%%%%%%%%%%ode45ray.m%%%%%%%%%%%%%%%%%%%%%%%
function [tout, yout] = ode45(ypfun, t0, tfinal, y0, tol, trace, D)
%Modified version of ode45.m that allows reflection of a ray off a surface.
% I use the symbol # to mark off changes I have made.
%
%function [tout, yout] = ode45(ypfun, t0, tfinal, y0, tol, trace)
%
%# We have to include ocean depth in our new ode solver
%
%ODE45	Solve differential equations, higher order method.
%	ODE45 integrates a system of ordinary differential equations using
%	4th and 5th order Runge-Kutta formulas.
%	[T,Y] = ODE45('yprime', T0, Tfinal, Y0) integrates the system of
%	ordinary differential equations described by the M-file YPRIME.M,
%	over the interval T0 to Tfinal, with initial conditions Y0.
%	[T, Y] = ODE45(F, T0, Tfinal, Y0, TOL, 1) uses tolerance TOL
%	and displays status while the integration proceeds.
%
%	INPUT:
%	F     - String containing name of user-supplied problem description.
%	        Call: yprime = fun(t,y) where F = 'fun'.
%	        t      - Time (scalar).
%	        y      - Solution column-vector.
%	        yprime - Returned derivative column-vector; yprime(i) = dy(i)/dt.
%	t0    - Initial value of t.
%	tfinal- Final value of t.
%	y0    - Initial value column-vector.
%	tol   - The desired accuracy. (Default: tol = 1.e-6).
%	trace - If nonzero, each step is printed. (Default: trace = 0).
%
%	OUTPUT:
%	T  - Returned integration time points (column-vector).
%	Y  - Returned solution, one solution column-vector per tout-value.
%
%	The result can be displayed by: plot(tout, yout).
%
%	See also ODE23, ODEDEMO.

%	C.B. Moler, 3-25-87, 8-26-91, 9-08-92.
%	Copyright (c) 1984-93 by The MathWorks, Inc.

% The Fehlberg coefficients:
alpha = [1/4  3/8  12/13  1  1/2]';
beta  = [ [    1      0      0     0      0    0]/4
    [    3      9      0     0      0    0]/32
    [ 1932  -7200   7296     0      0    0]/2197
    [ 8341 -32832  29440  -845      0    0]/4104
    [-6080  41040 -28352  9295  -5643    0]/20520 ]';
gamma = [ [902880  0  3953664  3855735  -1371249  277020]/7618050
    [ -2090  0    22528    21970    -15048  -27360]/752400 ]';
pow = 1/5;
if nargin < 5, tol = 1.e-6; end
if nargin < 6, trace = 0; end

% Initialization
t = t0;
hmax = (tfinal - t)/16;
h = hmax/8;
y = y0(:);
f = zeros(length(y),6);
chunk = 128;
tout = zeros(chunk,1);
yout = zeros(chunk,length(y));
k = 1;
tout(k) = t;
yout(k,:) = y.';

if trace
    clc, t, h, y
end

% The main loop
%# Added an '>=' expression to t+h>=t
zold=0;
hold=0;
while (t < tfinal) & (t + h >= t)
    if t + h > tfinal, h = tfinal - t; end
    
    % Aaron Thode correction...
    %#If ray is at boundry, reverse sign of
    %#dc/dz, or eta term, in state vector.
    %#It is important to put this section first.
    z=y(2); 		%#
    if (z>=D)		%#
        disp('z greater than D')
        %disp([y(4) z])		%#
        if y(4)>0
            y(4)=-y(4);
        end
       % h=hmax/8;
    elseif (z<=0)	%#
        disp('z is less than 0')
        %disp([y(4) z]);
        y(4)=abs(y(4));	%#
    end			%#
    
    % Compute the slopes
    temp = feval(ypfun,t,y);
    
    %#Correction for a step size that would go past a boundry
    ceta=temp(2); %# the second element of our derivative vector
    %# is c*eta or c*dz/ds
    if (z+h*ceta>D)	%#
        disp('Next step will exceed D')
        %disp([h z ceta])
        h=(D-z)/ceta;	%#
        disp(['h is now ' num2str(h)]);
        %pause
    elseif (z+h*ceta<0)%#
        %fprintf('Next step will surface, z is %6.2f,',z);
        
        %fprintf('h was %6.2f,',h);
        h=max([0.01 -z/ceta]);
        %fprintf('h is now %6.2f \n',h);
        if abs(hold-h)<0.01 && abs(zold-z)<0.1
            %keyboard
            %h=h/2;
        end
        hold=h;
        zold=z;
        
        %pause	%#
    elseif (z==D||z==0)
        %disp('Resetting h')
        %h=hmax/8
    end			%#
    
    f(:,1) = temp(:);
    for j = 1:5
        temp = feval(ypfun, t+alpha(j)*h, y+h*f*beta(:,j));
        f(:,j+1) = temp(:);
    end
    
    % Estimate the error and the acceptable error
    delta = norm(h*f*gamma(:,2),'inf');
    tau = tol*max(norm(y,'inf'),1.0);
%     while delta>=tau
%         h=h/2;
%         delta = norm(h*f*gamma(:,2),'inf');
%         tau = tol*max(norm(y,'inf'),1.0);
%     end
    
    % Update the solution only if the error is acceptable
    % Problem, may cause hang up
    if delta <= tau
        t = t + h;
        y = y + h*f*gamma(:,1);
        k = k+1;
        if k > length(tout)
            tout = [tout; zeros(chunk,1)];
            yout = [yout; zeros(chunk,length(y))];
        end
        tout(k) = t;
        yout(k,:) = y.';
    end
    if trace
        home, t, h, y
    end
    
    % Update the step size
    if delta ~= 0.0 && ~isnan(delta)
        h = min(hmax, 0.8*h*(tau/delta)^pow);
    end
    
end

if (t < tfinal)
    disp('Singularity likely.')
    t
end

tout = tout(1:k);
yout = yout(1:k,:);

