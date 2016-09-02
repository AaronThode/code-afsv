% Function to find all the peaks within a given signal, relative to some
% specified peak width, and SNR - noise is taken to be the median of the
% entire signal
% Modified from Matt Dzieucheu's pkpfit code
%
% Jit Sarkar
% MPL/SIO
% 06/29/2006
%
%
% Function form:
%	Peaks	=	peak_finder(Signal, Time, b_width, SNR)
%
%
% Input Variables:
%	Signal	=	[dbl/cplx](t)	Signal to be searched for peaks
%	Time	=	[dbl](t)		Time axis for Signal				{s}
%	b_width	=	[dbl](1)		Full width of signal in frequency	{Hz}
%	SNR		=	[dbl](1)		Ratio of peak amplitude to noise floor
%
% Output variables:
%
%	Peaks	=	[struct](np)	Structure containing peaks and all relevant
%								information

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Changed to include a Time base instead of sample rate
% Must modify data processing functions accordingly
function	[Peaks, NL]	=	peak_finder(Signal, Time, b_width, SNR, NL)



%%	Global variables for debugging
%global	DJ_debug;


iwidth	=	3;
%%	derive required quantities from Time axis (assume regularly spaced)
dt		=	mean(diff(Time));
tmin	=	min(Time);

%	find pulse width in time and samples
p_width	=	1./(b_width);		%	I specify bandwidth as the full width
p_width	=	floor(p_width / dt);
if	~(p_width > 0)
	error('Sample rate is too small to resolve this bandwidth');
end


%%	We need to use the envelope for this operation, so...
%	Check if it's complex
%if	~isreal(Signal)
%	Sig_env	=	abs(Signal);
%else	%	Otherwise use hilbert transform to get envelope
	Sig_env	=	envelope(Signal);
%end

%	Define noise level
if	~exist('NL','var') || isempty(NL)
NL	=	median(Sig_env);
if	~(NL	>	0)
	NL	=	mean(Sig_env);
	disp('JIT:noise','Noise floor median is 0, using mean instead');
end
end

%	Define search space reduced by 1 pulse width, so that you don't grab
%	spurius results at start and end
Sig_0	=	Sig_env(1+p_width:end-p_width);

%	Define logical indexing vector of acceptable points
Pks_lgcl	=	Sig_0 > SNR*NL;

%	a local maximum has to be larger than all the points around it within
%	the peak width
for	ii	=	1:p_width
	Sig_prev	=	Sig_env(1+p_width-ii : end-p_width-ii);
	Sig_next	=	Sig_env(1+p_width+ii : end-p_width+ii);
	
	Pks_lgcl	=	Pks_lgcl & (Sig_0 > Sig_prev) & (Sig_0 > Sig_next);
end
%	Convert logical indexing into real indices
Pks_i	=	find(Pks_lgcl) + p_width;		%	for indexing into Signal
%N_p		=	length(Pks_i);					%	# of peaks found
% should double check that it's length 1?
%	Sort by magnitude of peaks
% [B, IX]	=	 sort(Sig_env(Pks_i), 'descend');
% Pks_i	=	Pks_i(IX);
% clear	B IX;

%	Find peak position in time relative to start of signal
pk_Times	=	dt*(Pks_i - 1);

iold = 0;
if (iold == 1)
  % old way
%	Correct time and and value using parabolic fitting
%	p	=	a + b*t + c*t^2
t0	=	Pks_i;
p0	=	Sig_env(t0);		%	Peaks.Hill(:,2)	=	Signal(t0);
p1	=	Sig_env(t0-iwidth);	%	Peaks.Hill(:,1)	=	Signal(t0-1);
p2	=	Sig_env(t0+iwidth);	%	Peaks.Hill(:,3)	=	Signal(t0+1);

%	Make times relative to defined peak
%t0	=	0;
t1	=	-iwidth;
t2	=	+iwidth;

%	Work out equation constants
a	=	p0;
%	c	=	( (p1./t1 - p2./t2) - a.*(1./t1 - 1./t2) ) ./(t1 - t2);
b	=	( (t2./t1).*(p1-p0) - (t1./t2).*(p2-p0) ) ./ (t2-t1);
%	b	=	(p1 - a - c.*(t1.^2))./t1;
c	=	(p2 - a - b.*t2)./ (t2.^2);
else
% new: 
% workspace for more than 3 time points:
%	Correct time and and value using parabolic fitting
%	p	=	a + b*t + c*t^2
%       p_1     =   [ 1 t_1 t_1^2 ]
%       p_2     =   [ 1 t_2 t_2^2 ]
%             ....
%       p_n     =   [ 1 t_n t_n^2 ]
% d = Gm; m = (G'G)\G'd
% 5 unique elements of gram matrix (G'G)
%      g1 = n;
%      g2 = sum(t);
%      g3 = sum(t.*t);
%      g4 = sum(t.*t.*t);
%      g5 = sum(t.*t.*t.*t);

% choose a range of times around the max sample (same for all peaks)
it0 = [-iwidth:1:iwidth].';

% data are that range of pressure values for each peak: makes a matrix
% length(it0) x length(Pks_i);
nt0 = length(it0);
npks = length(Pks_i);
% make the gram matrix elements
g1 = nt0;
g2 = sum(it0);
g3 = sum(it0.*it0);
g4 = sum(it0.*it0.*it0);
g5 = sum(it0.*it0.*it0.*it0);
gtg = [g1 g2 g3; g2 g3 g4; g3 g4 g5];
%gtgi = inv(gtg);
% initialize parameter vectors: same names as old way
a = zeros(npks,1);
b = zeros(npks,1);
c = zeros(npks,1);
% make the projected data; have to do a loop
% !!try to do without loop; all by matrices
%it = repmat(it0',1,npks) + repmat(Pks_i,nt0,1);
%pt = Sig_env(it);
% compute G'*d for all peaks
for ipk = 1:npks
  it = it0 + Pks_i(ipk);
  pt = Sig_env(it);
  % compute G'*d
  gtd1 = sum(pt);
  gtd2 = sum(it0.*pt);
  gtd3 = sum(it0.*it0.*pt);
  % estimate a,b,c
  phat = gtg \ [gtd1; gtd2; gtd3];
  a(ipk) = phat(1);
  b(ipk) = phat(2);
  c(ipk) = phat(3);
end
% end of if-then
end
%	dp/dt	=	b + 2*c*t  ( = 0 at max)
tm			=	-b ./ (2*c);
% now compute peak times and values
pk_Vals		=	a + b.*tm + c.*(tm.^2);
pk_Times	=	pk_Times + dt.*tm;
% for angle, just take the value at the peak.
% could also average???
pk_Theta	=	angle(Signal(Pks_i));
  
%	matt's code, does the same thing......
% pk_Times	=	pk_Times + dt*(p1 - p2)./(p1 - 2*p0 + p2)/2;
% pk_Vals		=	p0 - (p2 - p1).*(p2 - p1)./(p2 - 2*p0 + p1)/8;


%	Copy outputs to return structure
Peaks.Index		=	Pks_i;						%	Indices of peaks
Peaks.Time_z	=	pk_Times + tmin;			%	Interpolated time
Peaks.Z			=	pk_Vals.*exp(1i*pk_Theta);	%	Interpolated value
Peaks.P			=	Signal(Peaks.Index);		%	Peak
Peaks.Time		=	dt*(Peaks.Index - 1) + tmin;%	indexed time
Peaks.dt		=	dt;							%	sample time step size
Peaks.Hill(:,1)	=	Signal(Peaks.Index-1);		%	Pre-peak point
Peaks.Hill(:,2)	=	Signal(Peaks.Index);		%	Peak
Peaks.Hill(:,3)	=	Signal(Peaks.Index+1);		%	Post-peak point
Peaks.SNR		=	SNR;						%	SNR used
Peaks.NL		=	NL;							%	Noise Level used

%	Holder variables for subsequent peak tracking through a series of
%	records
Peaks.tracks.prev_i		=	[];
Peaks.tracks.prev_cf	=	[];
Peaks.tracks.next_i		=	[];
Peaks.tracks.next_cf	=	[];




return;
