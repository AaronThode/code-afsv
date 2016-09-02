% DASAR_A_equalization.m

%%
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
f_hp = 10; % high-pass -3 dB break freq, Hz
[bb,aa] = butter(2,f_hp/(Fs/2),'high');

a_eqC2 = conv(a,aa);
b_eqC2 = conv(b,bb);

ff = [logspace(-2,log10(500),1000)]';
[h_eqC2,w] = freqz(b_eqC2,a_eqC2,ff*pi*(2/Fs));
figure
semilogx(w/pi *(Fs/2),20*log10(abs(h_eqC2)),'k')
grid on
xlabel('Frequency, Hz')
ylabel('gain, dB V/V')