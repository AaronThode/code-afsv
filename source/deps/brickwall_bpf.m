function [Hd,Bfilt] = brickwall_bpf(bandwidth,Fs,plot_on)
%  brickwall_bpf   designs a "brickwall" bandpass FIR filter (passband 10 to 450 Hz)
%   Hd = brickwall_bpf(Fs,plot_on) creates a sharp bandpass filter
%   at sampling frequency Fs (Hz)
%   if intermediate plot is desired, set plot_on = 1

% R Norman / Greeneridge Sciences 3/24/2011

% form a bandpass by cascading a lpf with a hpf
if length(bandwidth)==2
    fprintf('brickwall_bpf:  transition zone not specified, using default of 1 Hz');
    bandwidth=[bandwidth(1)-1 bandwidth bandwidth(2)+1];
elseif length(bandwidth)~=4
    disp('brickwall_bpf:  filter_bandwidth parameter length not set correctly in param.airgun');
end
% I changed this because with the AURALs and a Fs of 32kHz, the orders of
% firpm were ~46000 and the filters couldn't be made.  This keeps the
% transition_zone at 1Hz for DASARs, but increases to 33 Hz for AURALs.
%                           -Alex Conrad July 20, 2012
% Aaron changed this so that the default value is 1 Hz unless another value is substituted
Hd_lpf = lpf(Fs);
Hd_hpf = hpf(Fs);
Hd = dfilt.cascade(Hd_lpf,Hd_hpf);
if(plot_on)
    hFVT = fvtool(Hd)
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'brickwall bpf')
end

Bfilt = conv(Hd_lpf.Numerator,Hd_hpf.Numerator);


    function Hd = hpf(Fs)
        % HIGH-PASS
        % Equiripple Highpass filter designed using Parks-McClellan optimal
        % equiripple FIR filter design
        
        %Fpass = 9;                % Passband Frequency, Hz
        %Fstop = 10;                % Stopband Frequency, Hz
        
        Fstop = bandwidth(1);                  % Stopband Frequency, Hz
        Fpass = bandwidth(2);                 % Passband Frequency, Hz
        
        rs = 40;                    % Stopband attenuation, dB
        %rs=20;
        Dstop = 10^(-rs/20);
        
        rp = 1;                     % Passband ripple, dB (peak to trough)
        Dpass = (10^(rp/20)-1) / (10^(rp/20)+1);
        
        dens  = 20;                 % Density Factor
        
        % Calculate the order from the parameters using FIRPMORD.
        [N, Fo, Ao, W] = firpmord([Fstop, Fpass]/(Fs/2), [0 1], [Dstop, Dpass]);
        
        % Calculate the coefficients using the FIRPM function.
        b  = firpm(N, Fo, Ao, W, {dens});
        Hd = dfilt.dffir(b);
    end

    function Hd = lpf(Fs)
        % LOWPASS
        % Equiripple Lowpass filter designed using Parks-McClellan optimal
        % equiripple FIR filter design
        
        %Fpass = 450;                % Passband Frequency, Hz
        %Fstop = 451;                % Stopband Frequency, Hz
        
        %Fpass=bandwidth(2);
        %Fstop=Fpass+transition_zone;
        
        Fpass = bandwidth(3);                % Passband Frequency, Hz
        Fstop = bandwidth(4);                % Stopband Frequency, Hz
        
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
