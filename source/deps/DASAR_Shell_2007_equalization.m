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
