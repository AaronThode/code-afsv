function modes= warp_transform(plot_chc,s_ok,r_guess,c1,Fe, Nfft, N_window,N_window_w,t_trim, ...
    t_sweep,fstart,fend,spectro_mask_in,rd) % You need a long window to see the warped modes)
% The time-frequency representations that are used in this script are computed using a free
% toolbox that you can find here : http://tftb.nongnu.org/
% Be sure to add the tftb-01 folder to your matlab file if you want this to
% work !
% Inputs:
%  plot_chc:  If 1, make plots
%  s_ok: input signal, single channel
%  r_guess: initial guess of signal range, in meters
%  c1: initial estimate of water sound speed, m/s
%  Fe: sampling rate, Hz
%  Nfft: FFT size
%  N_window:  number of FFT window points for original time series
%  N_window_w: number of FFT window points for warped time series
%  t_trim: number of sec to remove from time series after source deconvolution estimate
%  fstart: vector of start frequencies of piecemeal linear FM sweeps, Hz
%  fend: vector of end frequencies of pieceal FM sweeps, Hz
%  t_sweep: duration of these sweeps, sec
%  spectro_mask_in:  masking window for filtering, externally provided from processing a different channel.
% pt: propagated signal in a Pekeris waveguide with parameters
%       c1, c2, rho1, rho2: sound speed / density
%       D: depth
%       r: range
%       zs, zr: source/receiver depth
% Fe: sampling frequency
%  Output:
%       modes:
%           s_w:  warped deconvolved input signal
%           s_filt: IFFT of mode-selections
%           spectro_w: spectrogram of s_w, evaluated for every time sample
%           spectro_mask: blocked out region of interest for inverse filtering.  3D arrays
%           time_w, freq_w: time and frequency scale for spectro_w
%           Fe_w:  Effective sampling rate in Hz for s_w


%N_window=31; % you need a short window to see the modes
%N_window_w=301; % You need a long window to see the warped modes
%t_trim=400;% Don't forget to change that if you consider other data

if size(s_ok,2)>1
    s_ok=s_ok.';
end

if ~exist('rd')
    rd=[];
end
N_ok=length(s_ok);
Nfft_ok=2^nextpow2(N_ok);
% Looking at the received signal, you can guess that the source is a piecewise linear FM
% downsweep. This is wrong but it shows you that you don't need to know
% the source exactly
%source_est=fmlin(160,0.5,0.01);

iflaw=[];
for I=1:length(fstart)
    iflaw_seg=linspace(fstart(I)/Fe,fend(I)/Fe,round(Fe*t_sweep(I)));
    iflaw=[iflaw; iflaw_seg(2:end)'];
end

%source_est=fmlin(round(t_sweep*Fe),fstart/Fe,fend/Fe);
source_est=fmodany(iflaw);

source_est_f=fft(source_est,Nfft_ok);
phi=angle(source_est_f);

% Phase correction for estimated source signal

%Trim front of signal to ensure first mode arrival is time zero, if necessary...
%s_ok=s_ok(round(Fe*t_trim):end);

sig_prop_f=fft(s_ok,Nfft_ok);
sig_rec_f=sig_prop_f.*exp(-1i*phi); % note that only source phase is deconvoluted
sig_rec_t=ifft(sig_rec_f,Nfft_ok); % Signal in time domain after source deconvolution


% STFT computation of original data
tfr=tfrstft(s_ok,1:N_ok,Nfft,hamming(N_window));
spectro=abs(tfr).^2;
% Time and frequency axis of the original signal
time=0:1/Fe:(N_ok-1)/Fe;
freq=0:Fe/Nfft:Fe-Fe/Nfft;

%Trim front of signal to ensure first mode arrival is time zero, if necessary...
sig_rec_t=sig_rec_t(round(Fe*t_trim):end);

tfr_sig_rec=tfrstft(sig_rec_t,1:N_ok,Nfft,hamming(N_window));

% Figure

if plot_chc==1
    figure(10)
    subplot(2,2,1);
    imagesc(time, freq, 10*log10(spectro));caxis([80 130])
    set(gca,'fontweight','bold','fontsize',14);
    axis xy
    ylim([0 Fe/2])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title(sprintf('Original signal at %2.0f m depth',rd))
    colorbar
    
    subplot(2,2,2)
    imagesc(time, freq, 20*log10(abs(tfr_sig_rec)))
    set(gca,'fontweight','bold','fontsize',14);
    
    axis xy
    ylim([0 Fe/2])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    title('Signal after source deconvolution')
    colorbar
    orient landscape
    caxis([80 130])
    %print('-djpeg',sprintf('DeconvSpectrogram_depth%2.0fm',rd));
    
end

for Iguess=1:length(r_guess)
    time=r_guess(Iguess)/c1:1/Fe:r_guess(Iguess)/c1+(N_ok-1)/Fe;
    
    % The warped signal will be s_w
    [s_w, Fe_w]=warp_temp(sig_rec_t,Fe,r_guess(Iguess),c1);
    M=length(s_w);
    
    tfr_w=tfrstft(s_w,1:M,Nfft,hamming(N_window_w));
    spectro_w=abs(tfr_w).^2;
    % Time and frequency axis of the warped signal
    time_w=0:1/Fe_w:(M-1)/Fe_w;
    freq_w=0:Fe_w/Nfft:Fe_w-Fe_w/Nfft; % this is actually not exact, but precise enough for a rough example
    % Figure
    
    if plot_chc==1
        subplot(2,2,3)
        imagesc(time_w, freq_w, spectro_w)
        set(gca,'fontweight','bold','fontsize',14);
        
        axis xy
        ylim([0 50])
        xlabel('Warped time (sec)')
        ylabel('Corresponding warped frequency (Hz)')
        title(sprintf('Warped signal square amp. assuming range %6.2f',r_guess(Iguess)))
        caxis([1e10 5e11]);colorbar
        
        subplot(2,2,4)
        imagesc(time_w, freq_w, 10*log10(abs(spectro_w)))
        set(gca,'fontweight','bold','fontsize',14);
        
        axis xy
        ylim([0 50])
        xlabel('Warped time (sec)')
        ylabel('Corresponding warped frequency (Hz)')
        title('dB value');
        caxis([80 120]);colorbar;
        
        orient landscape
        print('-djpeg','-r300',sprintf('WarpedSpectrogram_range%2.0fkm_depth%2.0fm',r_guess(Iguess)/1000,rd));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%

% To make it easier, filtering will be done by hand using the roipoly tool.
% See Matlab help for more information

% We work on a smaller TFR to make it easier

if ~exist('spectro_mask_in')||isempty(spectro_mask_in)
    %figure(5)
    spectro_mask_in=[];
    %RTF=abs(tfr_w).^2;%RTF=RTF/max(max(RTF));
    figure(10);
    %imagesc(time_w, freq_w, 10*log10(RTF));axis xy;caxis([90 120]); colorbar; ylim([0 50])
    BW_flag=0;
    Nroi=input('Enter number of modes to select, or return to skip:');
    modes.s_filt=[];
    modes.spectro_mask=[];
    if isempty(Nroi)
        return
    end
else
    BW_flag=1;
    spectro_mask=spectro_mask_in;
    Nroi=size(spectro_mask,3);
    fprintf('Input spectro_mask detected with %i modes.\n',Nroi);
end

s_filt=zeros(length(s_ok),Nroi);
for J=1:Nroi
    
    if BW_flag==0
        figure(10);
        BW0 = roipoly ;
        spectro_mask(:,:,J)=BW0;
    else
        BW0=squeeze(spectro_mask(:,:,J));
    end
    % To create a time-frequency mask, click several times on the spectrogram
    % to define the region you want to filter. The last click must be a
    % double-click to validate the mask creation
    %masque=zeros(size(tfr_w));
    %masque(1:n_limit-1,:)=spectro_mask;
    
    % Note that the mask is applied on the STFT (not on the spectrogram)
    mode_rtf_warp=BW0.*tfr_w;
    
    %%The trick below works if the FFT spectrogram is computed for every sample...
    [mode_temp_warp]=real(sum(mode_rtf_warp,1));
    s_filt(:,J)=iwarp_temp(mode_temp_warp,Fe_w,r_guess,c1,Fe,N_ok)';
    
    if plot_chc==1
        figure(11);
        
        subplot(Nroi,2,2*J-1)
        %RTF=abs(mode_rtf_warp/max(max(abs(tfr_w))));RTF=RTF./max(max(RTF));
        RTF=20*log10(abs(mode_rtf_warp));
        imagesc(time_w, freq_w, RTF);axis xy;ylim([0 50]);
        caxis([90 120]);colorbar
        
        subplot(Nroi,2,2*J);
        [~,FF,TT,BB] = spectrogram(s_filt(:,J),hanning(32),round(0.5*32),32,Fe);
        imagesc(TT,FF,10*log10(BB));%
        axis('xy');
        caxis([90 120]);colorbar
        title(sprintf('Recovered mode %i at depth: %6.2f m',J,rd));
        if J==Nroi
            orient tall
            print('-djpeg','-r300',sprintf('warp_tranform_%2.0fm',round(rd)))
        end
    end
end

modes.s_w=s_w;
modes.s_filt=s_filt;
modes.spectro_w=spectro_w;
modes.time_w=time_w;
modes.freq_w=freq_w;
modes.spectro_mask=spectro_mask;
modes.Fe_w=Fe_w;

end  %warp_transform
