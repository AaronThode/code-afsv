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
%       peak: peak amplitude of square of signal values (peak power)
%       t_20dB: peak width, computed by taking 20 dB points on either side of the peak
%       SE_20dB: Sound Exposure  (SE) computed using t20dB duration
%       t_Malme: peak width, computed by taking 5% and 95% levels of cumulative equalized SEL values.
%       SE_Malme:  SE computed using t_Malme duration.  Estimated
%           background noise subtracted.
%
%       rms_Malme: rms values =sqrt(SE_Malme/t_Malme):  Thus background
%           noise estimate has been subtracted.
%       %NOISE properties:
%       noise.duration: length of time of final noise sample.  Varies as
%               program tries to locate stationary intervel
%       noise.rms: rms level of noise
%       noise.SE: SE of noise inegrated over noise.duration
%       NOTE THAT RMS IS AN AMPLITUDE, THUS 20*log10(RMS) REQUIRED, WHILE SE
%           IS AN INTENSITY MEASURE, THUS 10*LOG10(SE)

function features=get_level_metrics_simple(x,Fs,bufferTime,debug,x2mean)

%persistent B

if size(x,2)>1
    x=x.';
end

features.msg='success';

features.peak=-1;
levels=[10 20];

for I=1:length(levels)
    fname=sprintf('p%idB',levels(I));
    
    features.(fname).duration=-1;
    features.(fname).SE=-1;
    features.(fname).tstart=-1;
    features.(fname).Istart=-1;
    features.(fname).Iend=-1;
    features.peak=-1;
end

features.t_Malme.duration=-1;
features.t_Malme.SE=-1;
features.t_Malme.rms=-1;

features.noise.rms=-1;
features.noise.SE=-1;
features.noise.duration=-1;


nbuff=min([length(x) round(bufferTime*Fs)]);
dt=1/Fs;
tt_all=dt*(1:length(x));

%Remove DC bias (we now assume bias has been removed)
if ~exist('x2mean','var')
    x=x-mean(x);
end
x2=x.^2;
xhilb=abs(hilbert(x)).^2;

[features.peak,Imax]=max(x2);
features.peak=sqrt(features.peak);

%Find level dB points
max_xhilb=max(xhilb);
for I=1:length(levels)
    fname=sprintf('p%idB',levels(I));
    I20p=Imax-1+find(xhilb(Imax:end)<=0.01*max_xhilb, 1 );
    I20m=find(xhilb(1:Imax)<=0.01*max_xhilb, 1, 'last' );
    if ~isempty(I20p)&~isempty(I20m)
        features.(fname).duration=dt*(I20p-I20m);
        features.(fname).SE=dt*trapz(x2(I20m:I20p));
    end
    features.(fname).Istart=I20p;
    features.(fname).Iend=I20m;
    features.(fname).tstart=tt_all(I20m);
    
end

%Estimate mean rms background noise level, including stationarity check
if ~exist('x2mean','var')
    [x2mean,Istart,Iadjust,x2eq_vec]=get_noise_est;
    
end

%%%%%%%Estimating rms pulse length%%%%%%%%%%
%Noise values...
cumSE=dt*cumsum(x2(Istart:end)-x2mean); %cumulative SE, subtracting out background noise estimate
cumSE=cumSE/max(cumSE);
%close all;

%[junk,Ip]=min(abs(cumSE-0.95));
%[junk,Im]=min(abs(cumSE-0.05));
%toll=0.005;
Ip=[];Im=[];
%toll=2*toll;
%Ip=Istart-1+min(find(abs(cumSE-0.95)<toll));
%Im=Istart-1+min(find(abs(cumSE-0.05)<toll));
Ip=Istart-1+find(cumSE>=0.95, 1 );

%%Robust estimation of Im works backwards from Ip
Im=Istart-1+find(cumSE(1:(Ip-Istart+1))<=0.05, 1, 'last' );
%Im=Istart-1+min(find(cumSE>=0.05));



if Ip<Im
    features.msg='cumSE decreasing';
    if exist('debug')==1
        
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



if exist('debug','var')
    if debug==1
        plot_Malme_calculation;
    elseif debug==2
        plot_cum_only;
    end
    
end

if isempty(Ip)||Ip<Fs*bufferTime
    features.msg='peak not reached in proper time window';
    return
end
if isempty(Im)||Im<bufferTime*Fs
    features.msg='pulse start begins in noise';
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Pulse duration and broadband SE and RMS calculation%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

features.t_Malme.duration=dt*(Ip-Im);
%April 23, 2009: Substract background noise from pulse.
%  Assumes noise is uncorrelated with pulse.
features.t_Malme.SE=dt*trapz(x2(Im:Ip))-features.t_Malme.duration*x2mean;
%features.t_Malme.SE=dt*trapz(x2(Im:Ip));
features.t_Malme.rms=sqrt(features.t_Malme.SE/features.t_Malme.duration);

features.noise.rms=sqrt(x2mean);
features.noise.duration=length(x2eq_vec)*dt;
features.noise.SE=features.t_Malme.duration*(features.noise.rms.^2);


% try
% [peakF]=get_FFT_metrics(x(1:Im),Fs);
% %features.noise.SE_FFT=SE_FFT;  %NO!  need to mutiply rms_FFT^2*features.t_Malme.duration;
% features.noise.peakF=peakF;
% catch
%    disp('get_level_metrics:  failure to compute noise from get_FFT_metrics');
% end
%
% try
% [peakF]=get_FFT_metrics(x(Im:Ip),Fs);
%
% features.peakF=peakF;
%
% catch
%    disp('get_level_metrics:  failure to compute signal from get_FFT_metrics');
% end


    function [peakF]=get_FFT_metrics(x,Fs)
        %function [SE_FFT,rms_FFT,SE_band,rms_band,freq_bandwidth,peakF]=get_FFT_metrics(x,bandwidth,freq_third_octave,Fs)
        
        peakF=-1;
        Nx=length(x);
        tpulse=Nx/Fs;
        
        Nfft=2^ceil(log10(Nx)/log10(2));
        
        Nzero=floor((Nfft-Nx-1)/2);
        xpad=[zeros(Nzero,1); hanning(Nx).*x; zeros(Nzero,1)];
        if length(xpad)~=Nfft
            xpad=[xpad; zeros(Nfft-length(xpad),1)];
        end
        X=fft(xpad,Nfft); %% March 25, 2011: Now padding zeros on each end
        
        wpwr=sum(hanning(Nfft).^2)/Nfft;
        normm=2./(wpwr*Nfft*Fs);  %3/24/11 made it two-side power spectral density
        PSDD=(abs(X).^2)*normm;
        
        
        
        [~,Ipeak]=max(PSDD);
        peakF=(Ipeak-1)*Fs/Nfft;
        
        
    end


    function [x2mean,Istart,Iadjust,x2eq_vec]=get_noise_est
        x2eq_vec=[];
        
        x2eq=(x2(1:nbuff));
        %Safety check
        x2eqhalf1=median(x2(1:round((nbuff/2))));
        x2eqhalf2=median(x2(round((nbuff/2)):nbuff));
        Iadjust=1;
        Istart=1;
        
        while nbuff>16&&abs(x2eqhalf1-x2eqhalf2)/x2eqhalf1>1  %Noise has to be within 3 dB of each other
            Iadjust=Iadjust+1;
            fprintf('nbuff is now %i\n',nbuff);
            if x2eqhalf1>x2eqhalf2 %bump in signal in last half of buffer
                Istart=Istart+round(nbuff/2);
            end
            x2eq=x2(Istart+(1:(nbuff/2)));
            nbuff=round(nbuff/2);
            bufferTime=nbuff/Fs;
            
            x2eqhalf1=median(x2eq(1:round(nbuff/2)));
            x2eqhalf2=median(x2eq(round(nbuff/2):end));
        end
        
        if nbuff<16
            features.msg='Equalization not stationary';
            return
        end
        
        x2eq_vec=x2eq;
        x2mean=mean(x2eq);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%plot_Malme_calculation inner function%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_cum_only
        
        tt=Istart/Fs+dt*(1:length(x2(Istart:end)));
        tt_all=dt*(1:length(x2));
        
        plot(tt,cumSE);%title('CUMULATIVE SE');
        hold on;
        set(gca,'fontweight','bold','fontsize',14);
        plot(tt(Ip-Istart),cumSE(Ip-Istart+1),'ro','markerfacecolor','r','markersize',8);
        plot(tt(Im-Istart),cumSE(Im-Istart+1),'go','markerfacecolor','g','markersize',8);
        xlimm=xlim;xlimm(1)=0;xlim(xlimm);grid on
        
        ylimm=ylim;ylimm(1)=-0.1;ylim(ylimm);
        
        xlabel('Time (s)');ylabel('I(t)/I_{max}(t)','interp','tex');
        grid on
        title(sprintf('I(t)=\\int{(p^{2}-p_{noise,rms}^{2})dt}, Duration is %6.2f sec',dt*(Ip-Im)),'interp','tex');
        plot(tt_all(I20m),cumSE(I10m-Istart+1),'gx','markerfacecolor','r','markersize',8);
        plot(tt_all(I20p),cumSE(I10p-Istart+1),'rx','markerfacecolor','g','markersize',8);
        hold off
        
    end

    function plot_Malme_calculation
        tt=Istart/Fs+dt*(1:length(x2(Istart:end)));
        
        figure(1)
        subplot(3,1,1);
        plot(tt,x(Istart:end));
        %specgram(x,128,Fs,[],96);%title(sprintf('tstart: %i',(Istart-1)/Fs));
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (s)');ylabel('Frequency (Hz)');
        xlimm=xlim;xlimm(1)=0;xlim(xlimm);grid on
        
        subplot(3,1,2);plot(tt,x2(Istart:end)-x2mean,tt_all,xhilb);title(sprintf('x2-x2mean, %i adjustments to x2mean',Iadjust));
        xlim(xlimm);%ylimm=ylim;ylimm(1)=0;ylim(ylimm);
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (s)');ylabel('p^{2}-p_{noise,rms}^{2} (Pa^2)','interp','tex');
        grid on
        %title(ctime2str(debug));
        
        hold on
        plot(tt_all(I20m),x2(I20m)-x2mean,'rx','markerfacecolor','r','markersize',8);
        plot(tt_all(I20p),x2(I20m)-x2mean,'gx','markerfacecolor','g','markersize',8);
        
        
        line([tt_all(I20m) tt_all(I20p)],[1 1]*features.peak.^2,'color','r')
        hold off
        
        subplot(3,1,3);plot(tt,cumSE);title('CUMULATIVE SE');
        hold on;
        set(gca,'fontweight','bold','fontsize',14);
        plot(tt(Ip-Istart),cumSE(Ip-Istart+1),'ro','markerfacecolor','r','markersize',8);
        plot(tt(Im-Istart),cumSE(Im-Istart+1),'go','markerfacecolor','g','markersize',8);
        xlim(xlimm);%ylimm=ylim;ylimm(1)=0;ylim(ylimm);
        xlabel('Time (s)');ylabel('I(t)/I_{max}(t)','interp','tex');
        grid on
        title(sprintf('I(t)=\\int{(p^{2}-p_{noise,rms}^{2})dt}, Duration is %6.2f sec',dt*(Ip-Im)),'interp','tex');
        plot(tt_all(I20m),cumSE(I20m-Istart+1),'rx','markerfacecolor','r','markersize',8);
        plot(tt_all(I20p),cumSE(I20p-Istart+1),'gx','markerfacecolor','g','markersize',8);
        
        
        %set(gcf,'pos',[ 119   664   560   420]);
        hold off
        disp('pause in get_level_metrics_simple');
        pause
    end


end



