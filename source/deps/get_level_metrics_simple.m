%%%get_level_metrics.m%%%%
%function features=get_level_metrics_simple(x,Fs,bufferTime,debug)
%   bandwidth, freq_third_octave removed
%
% Given a transient signal in a time series, output a series of measures
%  of signal 'level'
%  Aaron Thode Dec 29, 2013
% Input:
%   x: time series, assumed centered around zero (mean subtracted).  Units typically in uPa.  Assumed bandpass filtered.
%   Fs: sampling frequency, Hz
%   bufferTime: time in seconds used to estimate background noise level at start of x
%   debug: if exist, debug plots, assuming debug is a ctime
% Output:
%   If no values computed for a field, -1 is returned.
%   features: structure varible that containts following fields:
%       peak: peak amplitude of square of signal values (peak power)
%       t20dB: peak width, computed by taking 20 dB points on either side of the peak
%       SEL20dB: Sound Exposure Level (SEL) computed using t20dB duration
%       t_Malme: peak width, computed by taking 5% and 95% levels of cumulative equalized SEL values.
%       SEL_Malme:  SEL computed using t_Malme duration.  Estimated
%           background noise subtracted.
%
%       rms_Malme: rms values =sqrt(SEL_Malme/t_Malme):  Thus background
%           noise estimate has been subtracted.
%       %NOISE properties:
%       noise.duration: length of time of final noise sample.  Varies as
%               program tries to locate stationary intervel
%       noise.rms: rms level of noise
%       noise.SEL: SEL of noise inegrated over noise.duration
%       NOTE THAT RMS IS AN AMPLITUDE, THUS 20*log10(RMS) REQUIRED, WHILE SEL
%           IS AN INTENSITY MEASURE, THUS 10*LOG10(SEL)

function features=get_level_metrics_simple(x,Fs,bufferTime,debug)

%persistent B


features.peak=-1;
features.t20dB=-1;
features.SEL20dB=-1;
features.t_Malme=-1;
features.SEL_Malme=-1;
features.rms_Malme=-1;
features.msg='success';
features.noise.rms=-1;
features.noise.SEL=-1;
features.noise.duration=-1;


if size(x,2)>1
    x=x.';
end


nbuff=min([length(x) round(bufferTime*Fs)]);
dt=1/Fs;

%Remove DC bias
x=x-mean(x);
x2=x.^2;

[features.peak,Imax]=max(x2);
features.peak=sqrt(features.peak);

%Find 20 dB point
I20p=Imax-1+min(find(x2(Imax:end)<0.01*features.peak));
I20m=max(find(x2(1:Imax)<0.01*features.peak));
if ~isempty(I20p)&~isempty(I20m)
    features.t20dB=dt*(I20p-I20m);
    features.SEL20dB=dt*trapz(x2(I20m:I20p));
end

%Subtract background, assuming noise uncorrelated with signal
x2eq=(x2(1:nbuff));
%Safety check
x2eqhalf1=median(x2(1:round((nbuff/2))));
x2eqhalf2=median(x2(round((nbuff/2)):nbuff));
Iadjust=0;
Istart=1;


while nbuff>16&&abs(x2eqhalf1-x2eqhalf2)/x2eqhalf1>1  %Noise has to be within 3 dB of each other
    Iadjust=Iadjust+1;
    %disp(sprintf('nbuff is now %i',nbuff));
    if x2eqhalf1>x2eqhalf2, %bump in signal in last half of buffer
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

%Noise values...

cumSEL=dt*cumsum(x2(Istart:end)-x2mean); %cumulative SEL, subtracting out background noise estimate
cumSEL=cumSEL/max(cumSEL);
%close all;

%[junk,Ip]=min(abs(cumSEL-0.95));
%[junk,Im]=min(abs(cumSEL-0.05));
%toll=0.005;
Ip=[];Im=[];
%toll=2*toll;
%Ip=Istart-1+min(find(abs(cumSEL-0.95)<toll));
%Im=Istart-1+min(find(abs(cumSEL-0.05)<toll));
Ip=Istart-1+find(cumSEL>=0.95, 1 );

%%Robust estimation of Im works backwards from Ip
Im=Istart-1+find(cumSEL(1:(Ip-Istart+1))<=0.05, 1, 'last' );
%Im=Istart-1+min(find(cumSEL>=0.05));



if Ip<Im
    features.msg='cumSEL decreasing';
    if exist('debug')==1
        
        figure(1)
        subplot(3,1,1);specgram(x(Istart:end),128,Fs,[],96);
        title(ctime2str(debug));caxis([80 120]);
        
        tt=dt*(1:length(x2(Istart:end)));
        subplot(3,1,2);plot(tt,x2(Istart:end)-x2mean);title('x2-x2mean');
        subplot(3,1,3);plot(tt,cumSEL);title('CUMULATIVE SEL DECREASING');
        hold on;plot(tt(Ip-Istart+1),cumSEL(Ip-Istart+1),'go');plot(tt(Im-Istart+1),cumSEL(Im-Istart+1),'ro');
        pause;
        hold off
        
    end
    return
end



if exist('debug')==1
    plot_Malme_calculation;
    
end

if isempty(Ip)||Ip<Fs*bufferTime
    features.msg='peak not reached in proper time window';
    return
end
if isempty(Im)||Im<length(x2eq)/2
    features.msg='pulse start begins in noise';
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Pulse duration and broadband SEL and RMS calculation%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

features.t_Malme=dt*(Ip-Im);
features.SEL_Malme=dt*trapz(x2(Im:Ip))-features.t_Malme*x2mean;  %April 23, 2009: Substract background noise from pulse.
%features.SEL_Malme=dt*trapz(x2(Im:Ip));
features.rms_Malme=sqrt(features.SEL_Malme/features.t_Malme);

features.noise.rms=sqrt(x2mean);
features.noise.duration=length(x2eq_vec)*dt;
features.noise.SEL=features.t_Malme*(features.noise.rms.^2);


% try
% [peakF]=get_FFT_metrics(x(1:Im),Fs);
% %features.noise.SEL_FFT=SEL_FFT;  %NO!  need to mutiply rms_FFT^2*features.t_Malme;
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
        %function [SEL_FFT,rms_FFT,SEL_band,rms_band,freq_bandwidth,peakF]=get_FFT_metrics(x,bandwidth,freq_third_octave,Fs)
        
        peakF=-1;
        Nx=length(x);
        tpulse=Nx/Fs;
        
        Nfft=2^ceil(log10(Nx)/log10(2));
        
        Nzero=floor((Nfft-Nx-1)/2);
        xpad=[zeros(Nzero,1); hanning(Nx).*x; zeros(Nzero,1)];
        if length(xpad)~=Nfft;
            xpad=[xpad; zeros(Nfft-length(xpad),1)];
        end
        X=fft(xpad,Nfft); %% March 25, 2011: Now padding zeros on each end
        
        wpwr=sum(hanning(Nfft).^2)/Nfft;
        normm=2./(wpwr*Nfft*Fs);  %3/24/11 made it two-side power spectral density
        PSDD=(abs(X).^2)*normm;
        
    
        
        [~,Ipeak]=max(PSDD);
        peakF=(Ipeak-1)*Fs/Nfft;
        
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%plot_Malme_calculation inner function%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_Malme_calculation
        tt=Istart/Fs+dt*(1:length(x2(Istart:end)));
        
        figure(1)
        subplot(3,1,1);specgram(x,128,Fs,[],96);%title(sprintf('tstart: %i',(Istart-1)/Fs));
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (s)');ylabel('Frequency (Hz)');
        xlimm=xlim;xlimm(1)=0;xlim(xlimm);
        
        subplot(3,1,2);plot(tt,x2(Istart:end)-x2mean);title(sprintf('x2-x2mean, %i adjustments to x2mean',Iadjust));
        xlim(xlimm);ylimm=ylim;ylimm(1)=0;ylim(ylimm);
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Time (s)');ylabel('p^{2}-p_{noise,rms}^{2} (Pa^2)','interp','tex');
        grid on
        %title(ctime2str(debug));
        
        subplot(3,1,3);plot(tt,cumSEL);title('CUMULATIVE SEL');
        hold on;
        set(gca,'fontweight','bold','fontsize',14);
        plot(tt(Ip-Istart),cumSEL(Ip-Istart+1),'ro','markerfacecolor','r','markersize',8);
        plot(tt(Im-Istart),cumSEL(Im-Istart+1),'go','markerfacecolor','g','markersize',8);
        xlim(xlimm);ylimm=ylim;ylimm(1)=0;ylim(ylimm);
        xlabel('Time (s)');ylabel('I(t)/I_{max}(t)','interp','tex');
        grid on
        title(sprintf('I(t)=\\int{(p^{2}-p_{noise,rms}^{2})dt}, Duration is %6.2f sec',dt*(Ip-Im)),'interp','tex');
        
        set(gcf,'pos',[ 119   664   560   420]);
        hold off
      %pause
    end


end



