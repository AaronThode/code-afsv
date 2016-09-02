%	Function That returns the envelope of a time series (multiple
%	simultaneously
%
%	Jit Sarkar
%	MPL-SIO
%	07/17/2005
%
%
%	Input Variables:
%		x		:	Data series to be processed	[t,m]
%						time must be first dimension
%
%	Output Variables:
%		y		:	Envelope of input			[t,m]
%


function	[ya, yp]	=	envelope(x)

if isreal(x)
	y	=	hilbert(x);
else
	y	=	x;
end

ya	=	abs(y);
yp	=	angle(y);