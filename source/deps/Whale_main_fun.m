function [filt] = Whale_main_fun(x,Fs,r_guess,c1,NFFT,Nwindow1,Ncontour,Nhydro)

filt=struct('mask',[]);

set(0,'DefaultAxesFontName', 'Californian FB')
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultAxesFontWeight', 'bold')

% Overall parameters that may have to be changed depending on the whale call.

%NFFT_sel=512;           % Lower NFFT to make the time-frequency selection faster.

next=zeros(1,4)+1;      % A boolean vector used to go back and forward in the program

while next(1)
    sig_r=x(:,Nhydro);
    Nx=size(x,1);

%% Signal selection

    [x_dec,f_min,f_max]=sourceSelect_dB(sig_r,Fs,NFFT,Nwindow1);
    %clear sig_r
    next(1:end)=1;

%% Decimation
    % new sampling frequency is 2.1*fmax (with fmax=y_max*Fs as selected earlier)
%     decim_fact=ceil(1/2.1/y_max);
%     x_dec=decimate(x_sel,decim_fact);
%     clear x_sel

    % New parameters (because of the decimation)
    %Fe=Fs/decim_fact;
    
    Fe=Fs;
    
    %[y_min,y_max]=deal(decim_fact*y_min,decim_fact*y_max);
    %N_dec=length(x_dec);
    %Nwindow=1+2*ceil((ceil(Nwindow1/decim_fact))/2);
    
    freq=(0:NFFT-1)/NFFT*Fe;

    %% Source deconvolution

    while next(2)
        [x_deconv,source] = sourceDeconv_dB(x_dec,Fs,NFFT,Nwindow1,f_min,f_max,Ncontour);
        next(2:end)=1;
        filt.source=source;

    %% Time selection (for warping)

        while next(3)

            disp('Time selection')
            sig_whale=abs(tfrstft(x_deconv,1:length(x_deconv),NFFT,hamming(Nwindow1)));
            fig=imagescFun((1:Nx)/Fs,Fs*(1:NFFT)/NFFT,10*log10(sig_whale),[f_min f_max],...
                'Signal after source deconvolution: choose the first time instant (for warping) and last time instant');
            clear sig_whale
            [j,~]=ginput(2);
            close(fig)
            
            x_min=max(1,round(Fs*min(j)));
            x_max=min(length(x_dec),round(Fs*max(j)));
            x_ok=x_deconv(x_min:end);
            x_ok(x_max:end)=0;
            next(3:end)=1;
            disp('End of time selection')
            
            filt.lims=[x_min,x_max,f_min/Fs,f_max/Fs];

        %% Warping

            while next(4)
                N=length(x_ok);
                time=(0:N-1)/Fe;

                [modes,Nmode,choice] = modeSelect(x_ok,Fe,N,r_guess,c1,'samp',filt);


                if choice ~=0

                    if choice<4
                        next(choice+1:end)=0;
                        break
                    end

                    next(1:end)=0;

                %% Source reconvolution (add the initial source phase that was removed)

                    if Nmode~=0

                        modes_f=fft(modes,N,2);
                        modes_f_ok=zeros(Nmode,N);
                        source_phase=angle(fft(source,N));

                        for m=1:Nmode
                            modes_f_ok(m,:)=modes_f(m,:).*exp(1i*source_phase).';
                        end

                        modes_ok=ifft(modes_f_ok,N,2,'symmetric');
                        
                        filt.spectro_w=modes_f_ok;
                        clear modes_f_ok

                    %% Extraction FI

                        tm=zeros(Nmode,NFFT);

                        for m=1:Nmode  
                            [~, mode_rtfr,~]=tfrrsp(hilbert(modes_ok(m,:).'),1:N,NFFT,hamming(Nwindow1));
                            [tm(m,:),~]=momftfr(mode_rtfr,1,N,time);
                        end
                          
                    filt.s_filt=modes_ok;

                    %% Last figures
                        % Legend
                        leg=repmat(' ',Nmode,7);
                        figure,
                        for m=1:Nmode
                            leg_i=['Mode ' num2str(m)];
                            leg(m,1:length(leg_i))=leg_i;
                            subplot(ceil(Nmode/2),ceil(Nmode/2),m)
                            tfr_modes=tfrstft(modes_ok(m,:)',1:N,NFFT,hamming(Nwindow1));
                            imagesc(time,freq,abs(tfr_modes));axis('xy')
                            ylim([f_min f_max])
                            title(leg_i)
                        end

                        % Plot received TFR+dispersion curves
                        tfr=tfrstft(x_dec(x_min:end),1:N,NFFT,hamming(Nwindow));
                        TFR=abs(tfr);
                        FI=imagescFun(time,freq,10*log10(TFR),[y_min y_max]*Fe,'Received signal and extracted dispersion curves');
                        hold on
                        plot(tm,freq,'linewidth',2)
                        xlabel('Time [s]')
                        ylabel('Frequency [Hz]')
                        legend(leg)

                        temp=str2double(input(['\nAre you satisfied with the mode selection ? \n' ...
                        '1: Back to signal selection\n'...
                        '2: Back to source deconvolution\n'...
                        '3: Back to time selection (for warping)\n'...
                        '4: Save and exit\n'...
                        '5: Exit\n'],'s'));

                        if temp==4
                            save([soundsampPath soundsamp(1:end-4) '_filtered_modes.mat'])
                            try
                                saveas(FI,[soundsampPath soundsamp(1:end-4) '_FI'])
                            catch ME
                                disp('The figure was closed before it was saved')
                            end
                        end

                        if temp<4
                            next(1:temp)=1;
                            next(temp+1:end)=0;
                            break
                        else
                            next(1:end)=0;
                        end

                    end
                end
            end
        end
    end
end

end

