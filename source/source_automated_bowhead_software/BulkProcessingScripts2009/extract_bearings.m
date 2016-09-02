%function [thet,kappa,sd]=extract_bearings(y,bufferTime,Nfft,Fs,fmin,fmax,algchc,Nsamples)
%
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
%  Output: thet: angle in degrees, increasing clockwise from true north...,
%  sd standard deviation in degrees
function [thet,kappa,sd]=extract_bearings(y,bufferTime,Nfft,Fs,fmin,fmax,algchc,Nsamples)
kappa=[];sd=[];thet=[];
switch algchc
    case 'spectrogram'
        %In this case, y, is a template cut spectrogram.

        %        for I=1:3,
        %            subplot(3,1,I);
        %            imagesc(10*log10(abs(y{I})));
        %        end

        Cpvx=2*real(y{1}.*conj(y{2}) );
        Cpvy=2*real(y{1}.*conj(y{3}) );

        Igood=find(abs(Cpvx)>0);
       
        %   keyboard;
        [thet,kappa,sd] = get_vmests([Cpvx(Igood) Cpvy(Igood)],500);
%         subplot(3,1,1)
%         imagesc(10*log10(abs(y{1})));
%         subplot(3,1,2)
%         imagesc(10*log10(abs(Cpvx)));
%         subplot(3,1,3)
%         imagesc(10*log10(abs(Cpvy)));
%         pause;


    case 'sel_ratio'
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
        if bufferTime>0,
            x=x(ceil(1+bufferTime*Fs):(end-floor(bufferTime*Fs)),:);
        end

        vx=(x(:,1).*x(:,2));
        vy=(x(:,1).*x(:,3));
        %thet=atan2(sum(vx),sum(vy))*180/pi;
        
     
        [thet,kappa,sd]=get_vmests([vx vy],Nsamples);

    case 'sel_ratio_FFT'
        if bufferTime>0,
            x=y(ceil(1+bufferTime*Fs):(end-floor(bufferTime*Fs)),:);
        else
            x=y;
        end
        for I=1:3,
            x(:,I)=x(:,I)-mean(x(:,I));
        end

        Nfft_signal=2^nextpow2(size(x,1));


        X=fft(x,Nfft_signal);

        Ifrange=round((Nfft_signal/Fs)*[fmin fmax]);
        %Ifrange=max([1 Ifrange(1)]):min([Nfft_signal/2-1 Ifrange(2)]);

        %Try adding a buffer around frequencies...
        Ifrange=max([1 Ifrange(1)-2]):min([Nfft_signal/2-1 Ifrange(2)+2]);

        Cpvx=2*real(X(Ifrange,1)'*X(Ifrange,2));
        Cpvy=2*real(X(Ifrange,1)'*X(Ifrange,3));


        %Note, if vector is in form (x,y),
        % then atan2(x,y) will express bearing realtive to north( +y),
        %   increasing clockwise.  Negative angles increase
        %   counterclockwise.  Thus add 2*pi if resultant less than 0 to
        %   get traditional compass measure.
        thet=atan2(Cpvx,Cpvy)*180/pi;
    case 'CSDM'
        ovlap=7/8;
        if bufferTime>0,
            x=y(ceil(1+bufferTime*Fs):(end-floor(bufferTime*Fs)),:);
        else
            x=y;
        end
        for I=1:3,
            x(:,I)=x(:,I)-mean(x(:,I));
        end

        [Ksall,Ks_eig,Nsnap,freq_all,power,SNRest]=extractKsexact(x, ...
            ovlap,Nfft,1:3,[fmin fmax],Fs,-1,Nfft);

        Cpvx=2*real(sum(squeeze(Ksall(2,1,:))));
        Cpvy=2*real(sum(squeeze(Ksall(3,1,:))));

        %Cpvx=2*sum(real(squeeze(Ks_eig(2,1,:))));
        %Cpvy=2*sum(real(squeeze(Ks_eig(3,1,:))));

        thet=atan2(Cpvx,Cpvy)*180/pi;




end


end
