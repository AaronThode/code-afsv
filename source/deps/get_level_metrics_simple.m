%%%get_level_metrics_simple.m%%%%
%function features=get_level_metrics_simple(x,Fs,bufferTime,debug,x2mean)
%   bandwidth, freq_third_octave removed
%
% Given a transient signal in a time series, output a series of measures
%  of signal 'level'.

%  Aaron Thode Dec 29, 2013
% Input:
%   x: time series, assumed centered around zero (mean subtracted).  Units typically in uPa.
%           Assumed already bandpass filtered.
%   Fs: sampling frequency, Hz
%   bufferTime: time in seconds used to estimate background noise level at start of x
%   debug: if exist, debug plots, assuming debug is a ctime
%   x2mean: if exist, use this value instead of the computed value of noise
%           using bufferTime
% Output:
%   If no values computed for a field, -1 is returned.
%   features: structure varible that containts following fields:
%
%      msg: 'success'
%      peak: 2.1899e+06
%     p10dB: [1×1 struct]
%     p20dB: [1×1 struct]
%     cumSE: [1×1 struct]
%     noise: [1×1 struct]
%     tpeak: 0.1810
%       rms: [1×1 struct]
%       
%       msg:  indicates whether feature calculation was succesful, and if
%           not, why not.
%       peak: peak amplitude of pulse envelope (abs of hilbert xform)
%           (20*log10(peak)=peak power in dB re 1 uPa
%       tpeak: time (sec) measured from start of x (includes buffer) at
%           which peak occurs.
%       p10dB,p10dB: duration metrics computed by measuring 10 and 20 db down from
%           envelope peak
%       cumSE:  duration metrics computed by using cumulative sound exposure.
%       noise:  metrics of the noise sample in the buffer
%       rms:  root-mean-square time duration.
% Fields associated with p10dB, p20dB, and cumSE:
%          duration: 0.0100
%           SE: 8.8674e+09
%          rms: 9.4167e+05
%           SNR: 10
%           tstart: 0.1780
%           Istart: 178
%         Iend: 188
%         tend: 0.1880
%          FFT: [1×1 struct]
%
%       duration: estimated pulse duration in sec.
%       SE: cumulative sound exposure (int(x.^2-xnoise.^2));
%           (use  10*log10(SE)  to get dB re 1 uPa^2-s)
%       rms: RMS amplitude of pulse (is (SE/duration));
%               (actually square of rms, so use
%               10log10 to get dB re 1 uPa (rms))
%       SNR:  ratio of pulse rms level to noise rms level
%             (take 10log10(SNR) to get SNR in dB)
%       tstart, tend: start and end time of pulse in seconds, measured
%              from start of x (includes buffer).
%       Istart, Iend, indcies of tstart and tend in x
%       FFT:  peakF:  frequency of peak powerspectral density
%       
%NOISE properties:
%       noise.duration: length of time of final noise sample.  Varies as
%               program tries to locate stationary intervel
%       noise.rms: rms level of noise (actually square of rms, so use
%               10log10 to get dB re 1 uPa (rms)
%       noise.SE: SE of noise inegrated over noise.duration (use
%           10*log10(SE)  to get dB re 1 uPa^2-s
%       NOTE THAT RMS IS AN AMPLITUDE, THUS 20*log10(RMS) REQUIRED, WHILE SE
%           IS AN INTENSITY MEASURE, THUS 10*LOG10(SE)

function features=get_level_metrics_simple(x,Fs,bufferTime,debug,x2mean)

%persistent B
if ~exist('debug','var')
    debug=0;
end
strr='xs';
upper_energy_limit=0.95;

if size(x,2)>1 %make x a column vector
    x=x.';
end

features.msg='success';
Ibuffer=round(bufferTime*Fs);
dt=1/Fs;
tt=dt*(1:length(x));

features.peak=-1;
levels=[10 20];

for I=1:length(levels)
    fname=sprintf('p%idB',levels(I));
    
    features.(fname).duration=-1;
    features.(fname).SE=-1;
    features.(fname).rms=-1;
    features.(fname).tstart=-1;
    features.(fname).Istart=-1;
    features.(fname).Iend=-1;
    features.peak=-1;
end

features.cumSE.duration=-1;
features.cumSE.SE=-1;
features.cumSE.rms=-1;

features.noise.rms=-1;
features.noise.SE=-1;
features.noise.duration=-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Estimate duration using dB relative to peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove DC bias (we now assume bias has been removed)
if ~exist('x2mean','var')
    x=x-mean(x);
end

x2=x.^2;
xhilb=abs(hilbert(x)).^2;  %Envolope of pulse
[max_xhilb,Imax]=max(xhilb(Ibuffer:end));
Imax=Imax+(Ibuffer-1);
features.peak=sqrt(max_xhilb);
features.tpeak=tt(Imax);  %%%Peak time relative to start of x (includes buffer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate mean rms background noise level, including stationarity check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbuff=min([length(x) Ibuffer]);
[x2mean0,Istart,Iadjust,x2eq_vec]=get_noise_est;
if ~exist('x2mean','var')
    x2mean=x2mean0;
end
features.noise.rms=(x2mean);
features.noise.duration=length(x2eq_vec)*dt;
features.noise.SE=features.cumSE.duration*(features.noise.rms);

%%%%%%%%%%%%%%%%%%%%%%%%%
%Find level dB points
for I=1:length(levels)
    fname=sprintf('p%idB',levels(I));
    scale=(10.^(-levels(I)/10));
    I20p=Imax-1+find(xhilb(Imax:end)<=scale.*max_xhilb, 1 );
    I20m=find(xhilb(1:Imax)<=scale*max_xhilb, 1, 'last' );
    if ~isempty(I20p)&~isempty(I20m)
        features.(fname).duration=dt*(I20p-I20m);
        features.(fname).SE=dt*trapz(x2(I20m:I20p));
        features.(fname).rms=(features.(fname).SE/features.(fname).duration);
        features.(fname).SNR=features.(fname).rms/features.noise.rms;
         %features.(fname).SNRse=features.(fname).SE-features.noise.SE;
    end
    features.(fname).Istart=I20m;
    features.(fname).Iend=I20p;
    features.(fname).tstart=tt(I20m);
    features.(fname).tend=tt(I20p);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Estimate duration using rms measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

features.rms.duration=sqrt(trapz((tt.^2)*x2)./trapz(x2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Estimating sound exposure pulse length%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Noise values...
cumSE=zeros(size(x));
cumSE(Ibuffer:end)=dt*cumsum(x2(Ibuffer:end)-x2mean); %cumulative SE, subtracting out background noise estimate
cumSE=cumSE/max(cumSE);
%close all;

Ip=[];Im=[];
Ip=find(cumSE>=upper_energy_limit, 1 );
%%Robust estimation of Im works backwards from Ip
Im=find(cumSE(1:Ip)<=1-upper_energy_limit, 1, 'last' );

features.cumSE.Istart=Im;
features.cumSE.Iend=Ip;
features.cumSE.tstart=tt(Im);
features.cumSE.tend=tt(Ip);

if Ip<Im
    features.msg='cumSE decreasing';
    if exist('debug','var')
        
        figure(1)
        subplot(3,1,1);specgram(x(Istart:end),128,Fs,[],96);
        title(ctime2str(debug));caxis([80 120]);
        
        tt=dt*(1:length(x2(Istart:end)));
        subplot(3,1,2);plot(tt,x2(Istart:end)-x2mean);title('x2-x2mean');
        subplot(3,1,3);plot(tt,cumSE);title('CUMULATIVE SE DECREASING');
        hold on;plot(tt(Ip-Istart+1),cumSE(Ip-Istart+1),'go');plot(tt(Im-Istart+1),cumSE(Im-Istart+1),'ro');
        pause;
        hold off
        
    end
    return
end


if isempty(Ip)||Ip<Fs*bufferTime
    features.msg='peak not reached in proper time window';
    return
end
% if isempty(Im)||Im<bufferTime*Fs
%     features.msg='pulse start begins in noise';
%     return
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Pulse duration and broadband SE and RMS calculation%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

features.cumSE.duration=dt*(Ip-Im);
%  Assumes noise is uncorrelated with pulse.
%April 23, 2009: Subtract background noise from pulse.

features.cumSE.SE=dt*trapz(x2(Im:Ip))-features.cumSE.duration*x2mean;
features.cumSE.rms=(features.cumSE.SE/features.cumSE.duration);
features.cumSE.SNR=features.cumSE.rms/features.noise.rms;

%%%%%%%%%%%%%%%%%%%%
%%%autocorrelation
%%%%%%%%%%%%%%%%%%%%

[xc,laggs]=xcorr(x(Ibuffer:end),0.1*Fs);


if debug>0
    if debug==1
        plocumSE_calculation;
    elseif debug==2
        plot_cum_only;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Compute FFT metrics%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


metrics={'p10dB','p20dB','cumSE'};
Nfft=0;
for KK=1:length(metrics)
    Im=features.(metrics{KK}).Istart;
    Ip=features.(metrics{KK}).Iend;
    Nfft=max([Nfft Ip-Im]);
end
hplot=[];
signn=-1;
if debug>0
    signn=1;
end
for KK=1:length(metrics)
    Im=features.(metrics{KK}).Istart;
    Ip=features.(metrics{KK}).Iend;
    if isempty(Ip)||isempty(Im)||Ip-Im==1
        continue
    end
    try
        [features.(metrics{KK}).FFT.peakF,features.(metrics{KK}).FFT.rms_band,hplot(KK)]=get_FFT_metrics(x(Im:Ip),Fs,Nfft,signn*KK);
    catch
        disp('get_level_metrics:  failure to compute signal from get_FFT_metrics');
    end
end
hold off
if debug>0
    legend(hplot,{'cumSum','p10dB','p20dB'})
    disp('pause');
    pause;
end


    function [peakF,rms_band,hplot]=get_FFT_metrics(x,Fs,Nfft0,debug)
        
        peakF=-1;
        Nx=length(x);
        if exist('Nfft0','var')
            Nfft=2^ceil(log10(Nfft0)/log10(2));
        else
            Nfft=2^ceil(log10(Nx)/log10(2));
        end
        tpulse=Nx/Fs;
        
        
        Nzero=floor((Nfft-Nx-1)/2);
        xpad=[zeros(Nzero,1); hanning(Nx).*x; zeros(Nzero,1)];
        if length(xpad)~=Nfft
            xpad=[xpad; zeros(Nfft-length(xpad),1)];
        end
        X=fft(xpad,Nfft); %% March 25, 2011: Now padding zeros on each end
        
        wpwr=sum(hanning(Nfft).^2)/Nfft;
        normm=2./(wpwr*Nfft*Fs);  %3/24/11 made it two-sided power spectral density
        if 2*floor(Nfft/2)-Nfft~=0
            keyboard
        end
        PSDD=(abs(X(1:(Nfft/2))).^2)*normm;  %Units of power spectral density (e.g. dB 1 uPa^2/Hz)
        freq=((0:(Nfft/2-1))*Fs/Nfft)';
        
        [~,Ipeak]=max(PSDD);
        peakF=(Ipeak-1)*Fs/Nfft;
        rms_band=sqrt(trapz((freq.^2).*PSDD)./trapz(PSDD));
        
        if debug>0
            strr_color='krb';
            figure(2);
            hplot=plot(freq,10*log10(PSDD),strr_color(debug));grid on
            xlabel('Frequency(Hz)');ylabel('dB re 1 uPa2/Hz');
            hold on
            plot(freq(Ipeak),10*log10(PSDD(Ipeak)),['o' strr_color(debug)]);
            set(gca,'fontweight','bold','fontsize',14);
        else
            hplot=-1;
            
        end
        
    end


    function [x2mean,Istart,Iadjust,x2eq_vec]=get_noise_est
        Ntrials=4;
        x2eq_vec=[];
        x2mean=[];
        
        x2eq=(x2(1:nbuff));
        %Safety check
        x2eqhalf1=median(x2eq(1:round((nbuff/2))));
        x2eqhalf2=median(x2eq(round((nbuff/2)):end));
        Iadjust=0;
        Istart=1;
        %fprintf('nbuff is  %i, xeq1 and xeq2 are %6.2g, %6.2g\n',nbuff,x2eqhalf1,x2eqhalf2);
        
        while Iadjust<Ntrials&&abs(x2eqhalf1-x2eqhalf2)/min([x2eqhalf2 x2eqhalf1])>1  %Noise has to be within 3 dB of each other
            Iadjust=Iadjust+1;
            nbuff=round(nbuff/2);
            bufferTime=nbuff/Fs;
            
            if x2eqhalf1>x2eqhalf2 %If unwanted pulse in first half of data
                Istart=Istart+nbuff;
            end
            x2eq=x2(Istart+(1:nbuff));
            
            x2eqhalf1=median(x2eq(1:round(nbuff/2)));
            x2eqhalf2=median(x2eq(round(nbuff/2):end));
           % fprintf('nbuff is now %i, xeq1 and xeq2 now %6.2g, %6.2g\n',nbuff,x2eqhalf1,x2eqhalf2);
            
        end
        
        x2eq_vec=x2eq;
        x2mean=mean(x2eq);
        
        if Iadjust==Ntrials
            features.msg='Equalization not stationary';
           % disp(features.msg);
            return
        end
        
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%plocumSE_calculation inner function%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_cum_only
        
        plot(tt,cumSE);
        hold on;
        set(gca,'fontweight','bold','fontsize',14);
        plot(tt(Ip),cumSE(Ip),'ro','markerfacecolor','r','markersize',8);
        plot(tt(Im),cumSE(Im),'go','markerfacecolor','g','markersize',8);
        xlimm=xlim;xlimm(1)=0;xlim(xlimm);grid on
        
        xlabel('Time (s)');ylabel('I(t)/I_{max}(t)','interp','tex');
        grid on
        title(sprintf('I(t)=\\int{(p^{2}-p_{noise,rms}^{2})dt}, Duration is %6.2f sec',dt*(Ip-Im)),'interp','tex');
        for II=1:length(levels)
            fname=sprintf('p%idB',levels(II));
            plot(features.(fname).tstart,cumSE(features.(fname).Istart),[strr(II) 'g'],'markerfacecolor','w','markersize',8);
            plot(features.(fname).tend,cumSE(features.(fname).Iend),[strr(II) 'r'],'markerfacecolor','w','markersize',8);
        end
        
        %set(gcf,'pos',[ 119   664   560   420]);
        hold off
        
    end

    function plocumSE_calculation
        
        figure(1)
        subplot(4,1,1);
        plot(tt,x);hold on;
        line(tt(Ibuffer)*[1 1],[min(x(Ibuffer:end)) max(x(Ibuffer:end))],'color','k','linewidth',5)
        %specgram(x,128,Fs,[],96);%title(sprintf('tstart: %i',(Istart-1)/Fs));
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (s)');ylabel('Amplitude (uPa)');
        xlimm=xlim;xlimm(1)=0;xlim(xlimm);grid on
        hold off
        
        %%%
        subplot(4,1,2);plot(tt,x2-x2mean,tt,xhilb);
        title(sprintf('x2-x2mean(black) and abs(hilbert(x)) (red), %i adjustments to x2mean',Iadjust));
        xlim(xlimm);%ylimm=ylim;ylimm(1)=0;ylim(ylimm);
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (s)');ylabel('p^{2}-p_{noise,rms}^{2} (Pa^2)','interp','tex');
        grid on
        %title(ctime2str(debug));
        
        hold on
        
        for II=1:length(levels)
            fname=sprintf('p%idB',levels(II));
            plot(features.(fname).tstart,x2(features.(fname).Istart)-x2mean,[strr(II) 'g'],'markerfacecolor','g','markersize',6);
            plot(features.(fname).tend,x2(features.(fname).Iend)-x2mean,[strr(II) 'r'],'markerfacecolor','r','markersize',6);
        end
        
        %line([tt_all(I20m) tt_all(I20p)],[1 1]*features.peak.^2,'color','r')
        hold off
        
        %%%
        subplot(4,1,3);
        plot(tt,cumSE);
        hold on;
        set(gca,'fontweight','bold','fontsize',14);
        plot(tt(Ip),cumSE(Ip),'ro','markerfacecolor','r','markersize',8);
        h(1)=plot(tt(Im),cumSE(Im),'go','markerfacecolor','g','markersize',8);
        xlim(xlimm);%ylimm=ylim;ylimm(1)=0;ylim(ylimm);
        xlabel('Time (s)');ylabel('I(t)/I_{max}(t)','interp','tex');
        grid on
        title(sprintf('I(t)=\\int{(p^{2}-p_{noise,rms}^{2})dt}, Duration is %6.2f sec',dt*(Ip-Im)),'interp','tex');
        for II=1:length(levels)
            fname=sprintf('p%idB',levels(II));
            h(II+1)=plot(features.(fname).tstart,cumSE(features.(fname).Istart),[strr(II) 'g'],'markerfacecolor','w','markersize',8);
            plot(features.(fname).tend,cumSE(features.(fname).Iend),[strr(II) 'r'],'markerfacecolor','w','markersize',8);
        end
        hold off
        legend(h,{'cumSum','p10dB','p20dB'})
        
        subplot(4,1,4)
        plot(laggs/Fs,xc/max(abs(xc)),'k');grid on
        xlabel('Lag (sec)');
        set(gca,'fontweight','bold','fontsize',14);
        hold off
        
        %disp('pause in get_level_metrics_simple');
        %pause
    end


end



