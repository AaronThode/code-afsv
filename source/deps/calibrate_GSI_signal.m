function x=calibrate_GSI_signal(xin, keyword,RawFileName)

%calibrate_GSI_signal.m
% Convert raw A/D value of presure sensor in DASAR to uPa
% using sensitivity of 150 dB re 1uPa/V
% and peak A/D voltage of 2.5 V
%keyboard

if isempty(xin)
    x=[];
    return
end
if strcmp(keyword,'short')||~isempty(strfind(keyword,'DASAR2007'))
    [numd,dend] = DASAR_Shell_2007_equalization(1000,0);
    filt.a=dend;
    filt.b=numd;
    amp_Scale = (2.5/65535)*(10^(149/20));
    Nchan=size(xin,2);
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
        %x(:,I) = amp_Scale*filter(filt.a,filt.b,xt);
    end
elseif strcmp(keyword,'short')||~isempty(strfind(keyword,'DASARC'))
    %     x=xin.*(2.5/(2^16-1));  %x in Volts
    Nchan=size(xin,2);
    %      %%This is a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    amp_Scale = (2.5/65535)*(10^(149/20));
    
    for I=1:Nchan
        %xt=xin(:,I)-median(xin(:,I));%disp('subratc medi')
        xt=xin(:,I)-(2^15);
        x(:,I) = amp_Scale*filter(filt.b,filt.a,xt);
    end
elseif strcmp(keyword,'filter')
    error('calibrate_GSI_signal: filter no longer a valid keyword...')
elseif ~isempty(strfind(keyword,'NorthStar08')),
    
    %%This has a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    if ~isempty(strfind(RawFileName,'NS08A0'))||~isempty(strfind(RawFileName,'NA08Cx'))
        amp_Scale = (2.5/65535)*(10^(134/20));
    else
        amp_Scale = (2.5/65535)*(10^(148.8/20));
    end
    %amp_Scale=1;
    x = amp_Scale*filter(filt.b,filt.a,xin);
    %
    %     The Northstar deployments consist of 14 units, 12 of which are DASAR-Cs (built in 2008), and the other two DASAR-As (built in 2003).
    % These DASAR-As are identical to the Liberty DASAR-As with which you are already familiar.
    % By location, the breakdown is as follows (all units are DASAR-C unless marked otherwise):
    %
    % NA08A0 SN45
    % NA08B0 SN51
    % NA08C0 SN36
    % NA08D0 SN37
    % NA08E0 SN48
    % NA08F0 SN47
    % NA08G0 SN65
    % NA08H0 SN52
    % NA08I0 SN49
    % NA08J0 SN50
    % NS08A0 SN2 DASAR-A
    % NS08B0 SN58
    % NS08C0 SN59
    % NA08Cx SN1 DASAR-A
    % Response equalization for both types of DASAR (in the band 10 to 450 Hz) is very similar, with the only difference being a scalar gain value.
    % The hydrophone to ADC gain is
    % -149 dB V/uPa @ 100 Hz for the DASAR-C
    % -134 dB V/uPa @ 100 Hz for the DASAR-A
    %
    % You have the code to equalization both of these. For your convenience I have copied the relevant email message below...
    %
    % ------------------------
    % Aaron -
    %
    % For equalization of the omni channel on the BP DASAR-Cs, at sample rate 1 kHz and good for 10 Hz to 450 Hz analysis, you can use the
    %same IIR filter and coefficient values as were used on the Liberty DASAR-As. The only major difference is a change in gross sensitivity:
    % the DASAR-Cs are about 1/5 as sensitive as the DASAR-As. The same equalization can be used for all DASAR-Cs. Measurements in the Greeneridge
    % in-air loudspeaker test box showed all 65 units to be fairly close in sensitivity, with deviations looking no worse than the uncertainty in the calibration itself.
    %
    % The DASAR-Cs start rolling off below about 4 Hz (with a high-pass characteristic), so I don't recommend doing any analysis below that frequency without different equalization. You should be OK with these coef values since you are staying between 10 and 450 Hz.
    %
    % This equalization as shown includes a second order high-pass filter with break freq 1 Hz. You have the code to change it to 10 Hz if desired.
    %
    % Again, only one line changes going from DASAR-A to DASAR-C...
    %
    % om = x * 2.5/65536; % convert from counts to V
    % oml = 10^(148.8/20)*om; % uPa    convert from V to uPa   (WAS: 134 for DASAR-As, IS: 148.8 for DASAR-Cs)
    %
    % % this cascades the integrator with a 1 Hz high pass to cut low freqs
    % % -0.2 db at 100 Hz (ideally 0 db)
    % % -8.3 db at 250 Hz (ideally 20*log10(100/250)= -7.96
    % b = [0.53503845807280  -0.98982114743468   0.37452692065096   0.08025576871092]; % 11/13/2005
    % a = [1.00000000000000  -2.99111429220165   2.98226788807059  -0.99115359586894]; % 11/13/2005
    % % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % oml =  filter(b,a,oml); % uPa     equalize
elseif ~isempty(strfind(keyword,'Liberty08')>0)
    %om = xin * 2.5/65535; % convert from counts to V
    %om = (10^(134/20))*om; % uPa    convert from V to uPa
    %
    %     % this cascades the integrator with a 1 Hz high pass to cut low freqs
    %     % -0.2 db at 100 Hz (ideally 0 db)
    %     % -8.3 db at 250 Hz (ideally 20*log10(100/250)= -7.96
    % filt.b = [0.53503845807280  -0.98982114743468   0.37452692065096   0.08025576871092]; % 11/13/2005
    %filt.a = [1.00000000000000  -2.99111429220165   2.98226788807059  -0.99115359586894]; % 11/13/2005
    
    
    %This has a 10 Hz high pass filter
    
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    amp_Scale = (2.5/65535)*(10^(134/20));
    %amp_Scale=1;
    x = amp_Scale*filter(filt.b,filt.a,xin);
    
    %     keyboard;
end

if 1==0,
    Fs=1000;
    Nfft=256;
    A=fft(filt.a,Nfft);
    B=fft(filt.b,Nfft);
    F=linspace(0.1,Fs,Nfft);
    F=F(1:(Nfft/2));
    W=20*log10(abs(B./A));
    W=W(1:(Nfft/2));
    freqz(filt.b,filt.a,Nfft,Fs)
    hold on;plot(F,W,'r')
    
    figure;semilogx(F,W);xlim([1 1000]);grid on
    
end
end
