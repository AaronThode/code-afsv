function	SD_metrics	=	gsi_shipdet_metrics(file_name, t_win, t_int, bwidth, f_range, Nfft, overlap)


%%	Functional form of ship detection metrics set,
%	to process a whole GSI file in one go
%
%	Jit Sarkar
%	MPL/SIO
%	2013/03/15
%
%
% Function form:
%	SD_metrics	=	gsi_shipdet_metrics(file_name, t_win, t_int, bwidth, f_range)
%
% Inputs:
%	file_name	-	{str}(...}	Full/relative path of GSI file to process
%	t_win		-	{dbl}(1)	Spectrogram window size, in seconds
%	t_int		-	{dbl}(1)	Time interval between each sample window
%	bwidth		-	{dbl}(1)	Bandwidth for peak detection in # of FFT bins
%	f_range		-	{dbl}(2)	Frequency range to limit statistical
%								operations too, in Hz
%								Default cuts off upper/lower 10%
%	Nfft		-	{int}(1)	FFT length
%	overlap		-	{dbl)(1)	FFT window overlap size, as fraction 0 to 1
% Outputs:
%


%%	Parameter checking


%%	Load GSI file
%	Load header info first
[~,~,head]	=	readgsi(file_name);

Tmin	=	head.ctbc;
Tmax	=	head.ctec;
Tlen	=	Tmax - Tmin;
Fs		=	head.Fs;
%DTstart	=	datenum(1970,1,1,0,0,Tmin);
keyword	=	'DASARC';


%	Read the whole file
[X, Tx, ~] =	readgsi(file_name,0,Tlen);
X	=	X.';
Tx	=	Tx.';

%	Calibrate data
X	=	calibrate_GSI_signal(X, keyword);

%	Keep only sound data
X	=	X(:,1);

%	Size of record
Ntx	=	length(Tx);


%%	Calculate required sample lengths
it_win	=	floor(Fs*t_win);
it_inc	=	floor(Fs*t_int);

%	# of records to be produced
Nit		=	floor(Ntx/it_inc)+1;

%	Preallocate storage
SD_metrics.T			=	Tx(1:it_inc:Ntx);
SD_metrics.P_tot		=	zeros(size(SD_metrics.T));
SD_metrics.Time_entropy =	SD_metrics.P_tot;
SD_metrics.Time_kurtosis=	SD_metrics.P_tot;
SD_metrics.Freq_entropy =	SD_metrics.P_tot;
SD_metrics.Freq_kurtosis=	SD_metrics.P_tot;
SD_metrics.avg_CorrCoef =	SD_metrics.P_tot;
SD_metrics.avg_Conf		=	SD_metrics.P_tot;
SD_metrics.Max_Peaks	=	cell(size(SD_metrics.T));
SD_metrics.Min_Peaks	=	SD_metrics.Max_Peaks;

if	Nit ~= length(SD_metrics.T)
	disp('???');
	keyboard;
end



%%	Loop over record and calculate statistics
it_start	=	1;
it			=	1;

while it_start < Ntx
	%	Fetch current segment
	it_end		=	it_start + it_win - 1;
	if it_end > Ntx
		it_end	=	Ntx;
	end
	sig		=	X(it_start:it_end);
	
	%	Compute spectrogram for current window
	[~,F,T,PSD] = spectrogram(sig,hanning(Nfft),round(overlap*Nfft),Nfft,Fs);

	B		=	abs(10*log10(abs(PSD)));
	
	%	Limited frequency range for statistical operations
	Nt	=	length(T);
	Nf	=	length(F);
	dF	=	mean(diff(F));
	f_min	=	f_range(1);
	f_max	=	f_range(2);
	[~,if_min]	=	min(abs(F - f_min));
	[~,if_max]	=	min(abs(f_max - F));
	if_min	=	if_min(1);
	if_max	=	if_max(1);
	
	%	Calculate total power within window
	PSD_avg		=	sum(PSD,2)/Nt;
	P_tot		=	sum(PSD_avg)*dF;
	SD_metrics.P_tot(it)	=	P_tot;
		
	%	Average power vs time
	P	=	sum(B,1)/Nf;
	[E, Emax]	=	sentropy(P);
	K			=	kurtosis(P);
	SD_metrics.Time_entropy(it)	=	E/Emax * 100;
	SD_metrics.Time_kurtosis(it)=	K;
	
	%	Average power vs freq
	P	=	sum(B,2)/Nt;
	[E, Emax]	=	sentropy(P(if_min:if_max));
	K			=	kurtosis(P(if_min:if_max));
	SD_metrics.Freq_entropy(it)	=	E/Emax * 100;
	SD_metrics.Freq_kurtosis(it)=	K;
	
	%	Find peaks within this sum spectrum and overlay
	NL	=	harmmean(P(if_min:if_max));
	SNR	=	1.05;
	[Max_Peaks, ~]	=	peak_finder(P(if_min:if_max), F(if_min:if_max), 1./(bwidth*dF), SNR, NL);
	Pm		=	1./P;
	NLm		=	harmmean(Pm(if_min:if_max));
	[Min_Peaks, ~]	=	peak_finder(Pm(if_min:if_max), F(if_min:if_max), 1./(bwidth*dF), SNR, NLm);
	
	%	Rename fields to be spectrogram appropriate
	Max_Peaks.Freq_z	=	Max_Peaks.Time_z;			%	Interpolated freq point
	Max_Peaks.Freq		=	Max_Peaks.Time;				%	Indexed freq point
	Max_Peaks.df		=	Max_Peaks.dt;				%	Sample freq step size
	Max_Peaks	=	rmfield(Max_Peaks, 'Time_z');
	Max_Peaks	=	rmfield(Max_Peaks, 'Time');
	Max_Peaks	=	rmfield(Max_Peaks, 'dt');
	SD_metrics.Max_Peaks{it}	=	Max_Peaks;
	
	Min_Peaks.Freq_z	=	Min_Peaks.Time_z;			%	Interpolated freq point
	Min_Peaks.Freq		=	Min_Peaks.Time;				%	Indexed freq point
	Min_Peaks.df		=	Min_Peaks.dt;				%	Sample freq step size
	Min_Peaks	=	rmfield(Min_Peaks, 'Time_z');
	Min_Peaks	=	rmfield(Min_Peaks, 'Time');
	Min_Peaks	=	rmfield(Min_Peaks, 'dt');
	SD_metrics.Min_Peaks{it}	=	Min_Peaks;
	
	%	Correlation coefficients
	[R, P]	=	corrcoef(abs(B));
	SD_metrics.avg_CorrCoef(it)		=	mean(R(:));
	SD_metrics.avg_Conf(it)			=	mean(P(:));
		
	
	%	Increment window position
	it_start	=	it_start + it_inc;
	it			=	it + 1;
	
end
end

%%	Subfunctions
function	[E, Emax]	=	sentropy(S)
	S	=	abs(S);
	S	=	S./sum(S);
	E	=	-sum(S .* log2(S));
	
	Emax	=	log2(length(S));
	
	if	~isreal(E)
		warning('Entropy value not real');
		keyboard;
	end
end

function x=calibrate_GSI_signal(xin, keyword,RawFileName)

%keyboard
if strcmp(keyword,'short')||~isempty(findstr(keyword,'DASAR2007'))
    [numd,dend] = DASAR_Shell_2007_equalization(1000,0);
    filt.a=dend;
    filt.b=numd;
    amp_Scale = (2.5/65535)*(10^(149/20));
    Nchan=size(xin,2);
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
        %x(:,I) = amp_Scale*filter(filt.a,filt.b,xt);
    end
elseif strcmp(keyword,'short')||~isempty(findstr(keyword,'DASARC'))
    %     x=xin.*(2.5/(2^16-1));  %x in Volts
    Nchan=size(xin,2);
    %      %%This is a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    amp_Scale = (2.5/65535)*(10^(149/20));
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
    end
elseif strcmp(keyword,'filter')
    error('calibrate_GSI_signal: filter no longer a valid keyword...')
elseif ~isempty(findstr(keyword,'NorthStar08')),
    
    %%This has a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    if ~isempty(findstr(RawFileName,'NS08A0'))|~isempty(findstr(RawFileName,'NA08Cx'))
        amp_Scale = (2.5/65535)*(10^(134/20));
    else
        amp_Scale = (2.5/65535)*(10^(148.8/20));
    end
    %amp_Scale=1;
    x = amp_Scale*filter(filt.b,filt.a,xin);
    %
    %     The Northstar deployments consist of 14 units, 12 of which are DASAR-Cs (built in 2008), and the other two DASAR-As (built in 2003).
    % These DASAR-As are identical to the Liberty DASAR-As with which you are already familiar.
    % By location, the breakdown is as follows (all units are DASAR-C unless marked otherwise):
    %
    % NA08A0 SN45
    % NA08B0 SN51
    % NA08C0 SN36
    % NA08D0 SN37
    % NA08E0 SN48
    % NA08F0 SN47
    % NA08G0 SN65
    % NA08H0 SN52
    % NA08I0 SN49
    % NA08J0 SN50
    % NS08A0 SN2 DASAR-A
    % NS08B0 SN58
    % NS08C0 SN59
    % NA08Cx SN1 DASAR-A
    % Response equalization for both types of DASAR (in the band 10 to 450 Hz) is very similar, with the only difference being a scalar gain value.
    % The hydrophone to ADC gain is
    % -149 dB V/uPa @ 100 Hz for the DASAR-C
    % -134 dB V/uPa @ 100 Hz for the DASAR-A
    %
    % You have the code to equalization both of these. For your convenience I have copied the relevant email message below...
    %
    % ------------------------
    % Aaron -
    %
    % For equalization of the omni channel on the BP DASAR-Cs, at sample rate 1 kHz and good for 10 Hz to 450 Hz analysis, you can use the
    %same IIR filter and coefficient values as were used on the Liberty DASAR-As. The only major difference is a change in gross sensitivity:
    % the DASAR-Cs are about 1/5 as sensitive as the DASAR-As. The same equalization can be used for all DASAR-Cs. Measurements in the Greeneridge
    % in-air loudspeaker test box showed all 65 units to be fairly close in sensitivity, with deviations looking no worse than the uncertainty in the calibration itself.
    %
    % The DASAR-Cs start rolling off below about 4 Hz (with a high-pass characteristic), so I don't recommend doing any analysis below that frequency without different equalization. You should be OK with these coef values since you are staying between 10 and 450 Hz.
    %
    % This equalization as shown includes a second order high-pass filter with break freq 1 Hz. You have the code to change it to 10 Hz if desired.
    %
    % Again, only one line changes going from DASAR-A to DASAR-C...
    %
    % om = x * 2.5/65536; % convert from counts to V
    % oml = 10^(148.8/20)*om; % uPa    convert from V to uPa   (WAS: 134 for DASAR-As, IS: 148.8 for DASAR-Cs)
    %
    % % this cascades the integrator with a 1 Hz high pass to cut low freqs
    % % -0.2 db at 100 Hz (ideally 0 db)
    % % -8.3 db at 250 Hz (ideally 20*log10(100/250)= -7.96
    % b = [0.53503845807280  -0.98982114743468   0.37452692065096   0.08025576871092]; % 11/13/2005
    % a = [1.00000000000000  -2.99111429220165   2.98226788807059  -0.99115359586894]; % 11/13/2005
    % % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % oml =  filter(b,a,oml); % uPa     equalize
elseif ~isempty(findstr(keyword,'Liberty08')>0)
    %om = xin * 2.5/65535; % convert from counts to V
    %om = (10^(134/20))*om; % uPa    convert from V to uPa
    %
    %     % this cascades the integrator with a 1 Hz high pass to cut low freqs
    %     % -0.2 db at 100 Hz (ideally 0 db)
    %     % -8.3 db at 250 Hz (ideally 20*log10(100/250)= -7.96
    % filt.b = [0.53503845807280  -0.98982114743468   0.37452692065096   0.08025576871092]; % 11/13/2005
    %filt.a = [1.00000000000000  -2.99111429220165   2.98226788807059  -0.99115359586894]; % 11/13/2005
    
    
    %This has a 10 Hz high pass filter
    
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    amp_Scale = (2.5/65535)*(10^(134/20));
    %amp_Scale=1;
    x = amp_Scale*filter(filt.b,filt.a,xin);
    
    %     keyboard;
end

if 1==0,
    Fs=1000;
    Nfft=256;
    A=fft(filt.a,Nfft);
    B=fft(filt.b,Nfft);
    F=linspace(0.1,Fs,Nfft);
    F=F(1:(Nfft/2));
    W=20*log10(abs(B./A));
    W=W(1:(Nfft/2));
    freqz(filt.b,filt.a,Nfft,Fs)
    hold on;plot(F,W,'r')
    
    figure;semilogx(F,W);xlim([1 1000]);grid on
    
end
end

function [numd,dend] = DASAR_Shell_2007_equalization(Fs,plot_on)
% DASAR_Shell_2007_equalization
%       numd=b; dend=a;
%   Returns filter coefficients of a digital equalization filter to flatten
%   the (very) low-frequency response of data recorded by the
%   omnidirectional channel in the Shell 2007 DASAR.
%
%   The 2007 Shell DASAR was a 4-channel unit, with a PZT flexural disk
%   hydrophone and three orthogonal geophones. The lower band edge of the
%   hydrophone channel is controlled by two cascaded (and independent)
%   high-pass filters, one formed by the shunt resistor across the PZT
%   ceramic, the other by the preamp (a single non-inverting opamp stage)
%
%   This equalization is only needed if it is desired to reconstruct data
%   below 5 (or 8 Hz): for the "usual" (10 Hz and above) analysis, this
%   equalization to the recorded data is not necessary. The user is advised
%   that attempting to equalize more than a decade below these frequencies
%   (i.e., 0.5 Hz) should be attempted with caution, since ambient noise
%   can be large (especially pressure fluctuations from surface waves in
%   shallow water), and self-noise of the instrument can become visible.
%
%   Unlike the DASAR-A and DASAR-Cs, which employed modified sonobuoy heads
%   having a differentiator-like response (rising ~6 dB/octave), the
%   in-band response of all four channels of the Shell 2007 DASAR were
%   flat, so no equalization is needed in the passband of 10 to 450 Hz

% R Norman Mar 22, 2011

%%
if(Fs ~= 1000)
    error('designed & tested for Fs = 1000 Hz only')
end
%plot_on=1;
%% equalization filter for the high-pass filter formed by PZT ceramic and shunt resistor
R1 = 2e6; % shunt resistance, Ohms
C1 = 15e-9; % ceramic capacitance, F
f1a = 1/(2*pi*R1*C1); % break frequency, Hz

p1a = 2*pi*f1a; % zero location, rad/s

% pole location, rad/s (here, arbitrarily placed one-decade below the zero)
% Don't place this pole at a frequency any lower than necessary)
c = 2.5;
p1b = p1a/c; % pole location, rad/s

nums1 = c * [1/p1a 1]; % s-plane numerator coefs
dens1 = [1/p1b 1]; % s-plane denominator coefs

[numd1,dend1] = bilinear(nums1,dens1,Fs);

if(plot_on)
    hFVT = fvtool(numd1,dend1);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'First Equalization Filter')
    pause
end

%% equalization filter for the high-pass filter formed by preamplifier (non-inverting gain stage)
R2 = 200; % resistance, Ohms
C2 = 100e-6; % capacitance, F
f2a = 1/(2*pi*R2*C2);

p2a = 2*pi*f2a; % zero location, rad/s

% pole location, rad/s (here, arbitrarily placed one-decade below the zero)
% Don't place this pole at a frequency any lower than necessary)
p2b = p2a/c; % pole location, rad/s

nums2 = c * [1/p2a 1];
dens2 = [1/p2b 1];

[numd2,dend2] = bilinear(nums2,dens2,Fs);

if(plot_on)
    hFVT = fvtool(numd2,dend2);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'Second Equalization Filter')
    pause
end

%% cascade the two filters
numd = conv(numd1,numd2);
dend = conv(dend1,dend2);

% future fancier
% H1 = dfilt.df2(numd1,dend1); % df2t for transposed
% H2 = dfilt.df2(numd2,dend2); % df2t for transposed
% Hcascade = dfilt.cascade(H1,H2)
% [n1,d1] = tf(Hcascade)
% fvtool(Hcascade)

%%
if(plot_on)
    % freqress(nums,dens,logspace(-3,2,1000))
    hFVT = fvtool(numd,dend);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'Final (composite) Equalization Filter')
    orient landscape
    print -djpeg DASAR2007_equalization.jpg
end

end
