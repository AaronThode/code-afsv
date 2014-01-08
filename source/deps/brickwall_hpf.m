function [Hd] = brickwall_hpf(bandwidth,Fs,plot_on)
%  brickwall_hpf   designs a "brickwall" lowpass FIR filter 
%   Hd = brickwall_hpf(Fs,plot_on) creates a sharp bandpass filter
%   bandwidth: two-element vector specifying transition band in Hz.
%   Fs: sampling frequency Fs (Hz)
%   if intermediate plot is desired, set plot_on = 1

% R Norman / Greeneridge Sciences 3/24/2011, Thode modified 12/29/2013

if length(bandwidth)~=2
    disp('brickwall_hpf:  filter_bandwidth parameter length not set correctly');
    return
end

Hd = hpf(Fs);
%Hd_hpf = hpf(Fs);

if(plot_on)
    hFVT = fvtool(Hd);
    set(hFVT,'NormalizedFrequency','off','Fs',Fs)
    legend(hFVT,'brickwall hpf')
end



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

    
end
