function [Hd] = brickwall_lpf(bandwidth,Fs,plot_on)
%  brickwall_lpf   designs a "brickwall" lowpass FIR filter 
%   Hd = brickwall_lpf(Fs,plot_on) creates a sharp bandpass filter
%   bandwidth: two-element vector specifying transition band in Hz.
%   Fs: sampling frequency Fs (Hz)
%   if intermediate plot is desired, set plot_on = 1

% R Norman / Greeneridge Sciences 3/24/2011, Thode modified 12/29/2013

if length(bandwidth)~=2
    disp('brickwall_lpf:  filter_bandwidth parameter length not set correctly');
    return
end

Hd = lpf(Fs);
%Hd_hpf = hpf(Fs);
%Hd = dfilt.cascade(Hd_lpf,Hd_hpf);
if(plot_on)
    hFVT = fvtool(Hd);
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'brickwall lpf')
end

%Bfilt = conv(Hd_lpf.Numerator,Hd_hpf.Numerator);


    function Hd = lpf(Fs)
        % LOWPASS
        % Equiripple Lowpass filter designed using Parks-McClellan optimal
        % equiripple FIR filter design
        
        %Fpass = 450;                % Passband Frequency, Hz
        %Fstop = 451;                % Stopband Frequency, Hz
        
        %Fpass=bandwidth(2);
        %Fstop=Fpass+transition_zone;
        
        Fpass = bandwidth(1);                % Passband Frequency, Hz
        Fstop = bandwidth(2);                % Stopband Frequency, Hz
        
        rs = 40;                    % Stopband attenuation, dB
        %rs=20;
        Dstop = 10^(-rs/20);
        
        rp = 1;                     % Passband ripple, dB (peak to trough)
        Dpass = (10^(rp/20)-1) / (10^(rp/20)+1);
        
        dens  = 20;                 % Density Factor
        
        % Calculate the order from the parameters using FIRPMORD.
        [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
        
        % Calculate the coefficients using the FIRPM function.
        b  = firpm(N, Fo, Ao, W, {dens});
        Hd = dfilt.dffir(b);
    end
end
