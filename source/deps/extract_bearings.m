function [thet, kappa, sd, x] = extract_bearings(app, y, bufferTime, Nfft, Fs, fmin, fmax, Nsamples, f_transition)
%function [thet,kappa,sd,x]=extract_bearings(y,bufferTime,Nfft,Fs,fmin,fmax,Nsamples)
% Input:
%    y: time series, with channels arranged as columns
%    bufferTime: how much time exists before and after signal proper.
%           Needed because time-domain filtering requires a signal buffer.
%    Nfft: FFT size to use if CSDM is to be estimated.
%    Fs: Sampling frequency, Hz
%    fmin,fmax: minimum and maximum frequency of signal, Hz.
%    algchc:  String containing one of three possibilities:
%           'spectrogram: make spectrogram, estimate bearing from each
%                   time-frequency cell of contour, provide standard
%                   deviation
%           'sel_ratio':  time-domain method that takes rms value of
%           bandpassed signals.
%           'sel_ratio_FFT': takes the FFT of all signal components and
%               computes active intensity in frequency domain.
%           'CSDM':  With a small Nfft value, construct multiple shapshots
%               of signal structure, yielding a CSDM that in turn gives
%               and active intensity.  May also be modified to permit
%               robust beamforming for weak signals.
%    Nsamples:  number of bootstrap samples for estimating kappa and sd
%   f_transition:  frequency over which the pressure/velocity are 90 out of
%       phase, so need to use a hilbert transform to get the bearing..
%
%  Output: thet: angle in degrees, increasing clockwise from true north...,
%  sd standard deviation in degrees
%   x: filtered time series, columns are channels
kappa=[];sd=[];thet=[];

filter_chc='butter';

if strcmp(filter_chc,'FIR')
    transband=0.1*(fmax-fmin); %Hz
    filter_min=max([0.5*transband 0.8*fmin]);
    filter_max=min([500-0.5*transband 1.2*fmax]);
    [n,fo,mo,w] = firpmord( [filter_min+0.5*transband*[-1 1] filter_max+0.5*transband*[-1 1]], [0 1 0], [0.01 0.1 0.01], Fs );
    b = firpm(n,fo,mo,w);

    % design and apply filter to pass only between e1 and e2

    for I=1:3,
        y(:,I)=y(:,I)-mean(y(:,I));
        x(:,I)=filtfilt(b,1,y(:,I));
    end

else
    w1=max([0.01 fmin*2/Fs]);
    w2=min( [fmax*2/Fs 0.99]);

    [B,A]=butter(2,[w1 w2]);
    for I=1:3
        y(:,I)=y(:,I)-mean(y(:,I));
        x(:,I)=filter(B,A,y(:,I));
    end
end
if bufferTime>0
    x=x(ceil(1+bufferTime*Fs):(end-floor(bufferTime*Fs)),:);
end

vx=(x(:,1).*x(:,2));
vy=(x(:,1).*x(:,3));
if fmin>=f_transition
    vx=(imag(hilbert(x(:,1))).*x(:,2));
    vy=(imag(hilbert(x(:,1))).*x(:,3));

end
%thet=atan2(sum(vx),sum(vy))*180/pi;


[thet,kappa,sd]=get_vmests([vx vy],Nsamples);

if isnan(thet)
    keyboard
end
end
