%calibrate_GSI_signal.m
% Convert raw A/D value of presure sensor in DASAR to uPa
% using sensitivity of 150 dB re 1uPa/V
% and peak A/D voltage of 2.5 V
%function x=calibrate_GSI_signal(xin, keyword,RawFileName)

function x=calibrate_GSI_signal(xin, keyword,RawFileName)

%%Check orientation
if size(xin,2)>100 %data is horizontal time series
    xin=xin';
end
%keyboard
if strcmp(keyword,'short')|~isempty(findstr(keyword,'DASAR2007'))
    [numd,dend] = DASAR_Shell_2007_equalization(1000,0);
    filt.a=numd;
    filt.b=dend;
    amp_Scale = (2.5/65535)*(10^(149/20));
    Nchan=size(xin,2);
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
    end
elseif strcmp(keyword,'short')|~isempty(findstr(keyword,'DASARC'))
    Nchan=size(xin,2);
    
    %      %%This is a 10 Hz high pass filter, computed using get_DASARA_filter below.
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    amp_Scale = (2.5/65535)*(10^(149/20));
    
    for I=1:Nchan
        xt=xin(:,I)-(2^15);%disp('subratc medi')
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

function  [a_eqC2, b_eqC2]=get_DASARA_filter(f_hp,plot_data)
% DASAR_A_equalization.m and DASAR_C_equalization
%  From Bob Norman, Oct. 16, 2008 11:33 PM...
% Aaron -
% 
% For equalization of the omni channel on the BP DASAR-Cs, at sample rate 1 kHz and good for 10 Hz to 450 Hz analysis, 
% you can use the same IIR filter and coefficient values as were used on the Liberty DASAR-As. 
% The only major difference is a change in gross sensitivity: the DASAR-Cs are about 1/5 as sensitive as the DASAR-As. 
% The same equalization can be used for all DASAR-Cs. Measurements in the Greeneridge in-air loudspeaker test box 
% showed all 65 units to be fairly close in sensitivity, with deviations looking no worse than the uncertainty in the calibration itself.
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
% oml = filter(b,a,oml); % uPa     equalize
% 
%  
% Regards
% 
% - Bob

%%
f_hp=10; %Hz.  Use 1 Hz to generate the b,a values in Bob's email...
Fs = 1000; % sample rate, Hz

% equalizer rev C
% used to equalize Sparton DIFAR sonobuoy head omni phone response
% The raw sonobuoy response is like a differentiator, with a +19 to +20dB/dec slope.
% This equalizer is mostly an integrator, with some shaping near Nyquist to get rid
% of aliasing effects and a high-pass filter cascaded to reduce low frequency gain

f1 = 100; % Hz
a1 = 1; % gain @ f1 (ignoring high-pass), V/V

% integrator
a = [1 -1]; % gain of integrator is cos(w/2)/sin(w), w = [0,2*pi] around unit circle

% hf zero
k = (f1/Fs)*2*pi;
b = a1 * sin(k)/cos(k/2) * [1 0.15]/1.15; % 0.15 eyeball adhoc
%b=a1;
% b = a1 * sin(k)/cos(k/2);

% high-pass
%f_hp = 10; % high-pass -3 dB break freq, Hz
[bb,aa] = butter(2,f_hp/(Fs/2),'high');

a_eqC2 = conv(a,aa);
b_eqC2 = conv(b,bb);

if plot_data==1,
    ff = [logspace(-2,log10(500),1000)]';
    [h_eqC2,w] = freqz(b_eqC2,a_eqC2,ff*pi*(2/Fs));
    figure
    semilogx(w/pi *(Fs/2),20*log10(abs(h_eqC2)),'k')
    grid on
    xlabel('Frequency, Hz')
    ylabel('gain, dB V/V')
end

function [numd,dend] = DASAR_Shell_2007_equalization(Fs,plot_on)
% DASAR_Shell_2007_equalization
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

%% equalization filter for the high-pass filter formed by PZT ceramic and shunt resistor
R1 = 2e6; % shunt resistance, Ohms
C1 = 15e-9; % ceramic capacitance, F
f1a = 1/(2*pi*R1*C1) % break frequency, Hz

p1a = 2*pi*f1a; % zero location, rad/s

% pole location, rad/s (here, arbitrarily placed one-decade below the zero)
% Don't place this pole at a frequency any lower than necessary)
c = 10;
p1b = p1a/c; % pole location, rad/s

nums1 = c * [1/p1a 1]; % s-plane numerator coefs
dens1 = [1/p1b 1]; % s-plane denominator coefs

[numd1,dend1] = bilinear(nums1,dens1,Fs);

if(plot_on)
    hFVT = fvtool(numd1,dend1);
    set(hFVT,'NumberofPoints',8192,'FrequencyScale','Log')
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'First Equalization Filter')
end

%% equalization filter for the high-pass filter formed by preamplifier (non-inverting gain stage)
R2 = 200; % resistance, Ohms
C2 = 100e-6; % capacitance, F
f2a = 1/(2*pi*R2*C2)

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
end

