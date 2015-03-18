function [cable_factor,sens]=get_ADAT24_cable_factor
%%Can we calibrate the data?
%%  load_wav often normalizes the data so the peak value is 1.
%sens0=157; %dB re 1 unit of wav entry
%sens=input('Enter sensitivity of entire system (flat-spectrum calibration) [180 dB re 1 unit]:');
%if isempty(sens)
%sens=sens0;
%end
%sens=10^(sens/20);

%%Calibration for ADAT 24 attached to Sonotech hydrophone array
%//ADAT HD24 has 6.9 V Rms max input for 24 bit data
%//Factor of 2 from differential inputs...
%//Hydrophone gain set to Sonotech array -157 dB
%// Don't have cable attenuation  here...
sens=(6.9*sqrt(2)/16777215)*0.5*10^(157.0 / 20.0);
%freq_cal[i]=(1.0+2*Math.PI*f*(110e-9)*140.0);

% Nov 15, 2011
%         Hi Aaron,
%
% I thumbed through my notes for a few minutes to refresh my memory on this project.  The attenuation problem can be hugely simplified by eliminating many of the parameters right away.  The preamp output impedance is ordinarily very small, and the cable insulation resistance and receiver input impedance will be very large.  All of the above parameters can generally be ignored.
%
% This really only leaves the cable resistance and capacitance and so becomes a fairly simple voltage divider problem.  From my notes, I measured the cable (the entire 800+ length) resistance to be ~70 ohms and the capacitance to be ~110 nF.  The model would look like a series resistance with a shunt capacitance, so:
%
% R = 140 ohms (round trip)
% X = 1/(2*pi*f*C)
%
% So, the 6 dB down frequency will be when R= X, so:
% f = 1/(2*pi*C*R) = ~10 kHz.
%
% And the attenuation at any frequency:
% X/(R + X)  or 1/(1+ 2*pi*f*C*R)
%
% The cable attenuation seen by the vector sensor will essentially look like a single pole low pass filter with Fc = 10 kHz:
% 5   kHz: -3 dB
% 10 kHz: -6 dB
% 20 kHz: -9 dB
%
% This is consistent with my notes where I measured the cable attenuation to roll off by 3 dB/octave starting at 4-5 kHz.  The above is for the case of the vector sensor driving the entire length of the array.  For the other hydrophones driving shorter lengths of cable, the calculations will be the same but you will need to use different R and C values in the equation 20 log (1/(1+ 2*pi*f*C*R)) to compute the new attenuation values vs frequency.
%
% Jeff

cable_factor=2*pi*(110e-9)*140.0;  %Unit resistance 140 ohm, capacitance 110 nF
end
