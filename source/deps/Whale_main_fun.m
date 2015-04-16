function [filt] = Whale_main_fun(x,Fs,r_guess,c1,NFFT,Nwindow,Ncontour,Nhydro,flims)

filt=struct('mask',[]);

set(0,'DefaultAxesFontName', 'Californian FB')
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultAxesFontWeight', 'bold')

% Overall parameters that may have to be changed depending on the whale call.

next=zeros(1,4)+1;      % A boolean vector used to go back and forward in the program

while next(1)
    sig_r=x(:,Nhydro);

%% Signal selection

    [x_dec,f_min,f_max,x_max]=sourceSelect_dB(sig_r,Fs,NFFT,Nwindow,1000*flims);
 
    clear sig_r
    next(1:end)=1;
    
    freq=(0:NFFT-1)/NFFT*Fs;
    
    %% Source deconvolution

    while next(2)
        [x_deconv,source] = sourceDeconv_dB(x_dec,Fs,NFFT,Nwindow,f_min,f_max,Ncontour);
        next(2:end)=1;
        filt.source=source;

    %% Time selection (for warping)

        while next(3)
            Nx=length(x_deconv);

            disp('Time selection')
            sig_whale=abs(tfrstft(x_deconv,1:Nx,NFFT,hamming(Nwindow)));
            fig=imagescFun((1:Nx)/Fs,Fs*(1:NFFT)/NFFT,10*log10(sig_whale),[f_min f_max],...
                'Signal after source deconvolution: choose the first time instant (for warping)');
            clear sig_whale
            
            j=ceil(Fs*ginput(1));
            close(fig)
            
            %% Choose a list of samples near the selected one
            dj=ceil(Nx/200);
            nj=10;
            j_list=max(1,j-nj*dj):dj:j+nj*dj;
            
            % Make a video
            F(length(j_list)) = struct('cdata',[],'colormap',[]);
            
            for jj = 1:length(j_list)
                if jj==1
                    mov=figure('units','normalized','outerposition',[0 0 1 1]);
                    [fw_min,fw_max,xmax]=choose_instant(x_deconv,j_list(jj),Fs,r_guess,c1,max(j_list));
                    ax=gca;
                    ax.NextPlot = 'replaceChildren';
                end
                    choose_instant(x_deconv,j_list(jj),Fs,r_guess,c1,xmax,[fw_min fw_max]);
                    F(jj) = getframe(gcf);
            end
            
            commandwindow;
            replay=input('Replay ? (y/n) ','s');
            
            while(replay~='n')
                shg
                movie(mov,F,1,1)
                commandwindow;
                replay=input('Replay ? (y/n) ','s');
            end
           
            commandwindow;
            x_min=str2double(input('\n Choose a time instant: ','s'));

            x_ok=x_deconv(x_min:end);
            disp('End of time selection')

            filt.lims=[x_min,x_max,f_min,f_max];

        %% Warping
        
            N=length(x_ok);
            time=(0:N-1)/Fs;

            [modes,Nmode,choice,filt] = modeSelect(x_ok,Fs,N,r_guess,c1,'samp',filt);

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
                        [~, mode_rtfr,~]=tfrrsp(hilbert(modes_ok(m,:).'),1:N,NFFT,hamming(Nwindow));
                        [tm(m,:),~]=momftfr(mode_rtfr,1,N,time);
                    end

                filt.s_filt=modes_ok;

                %% Last figures
                    % Legend
                    leg=repmat(' ',Nmode,7);
                    figure,
                    colorbar
                    for m=1:Nmode
                        leg_i=['Mode ' num2str(m)];
                        leg(m,1:length(leg_i))=leg_i;
                        subplot(ceil((Nmode+1)/2),ceil(Nmode/2),m)
                        tfr_modes=tfrstft(modes_ok(m,:)',1:N,NFFT,hamming(Nwindow));
                        imagesc(time,freq,abs(tfr_modes))
                        axis xy
                        ylim([f_min f_max])
                        title(leg_i)
                    end

                    % Plot received TFR+dispersion curves
                    tfr=tfrstft(x_dec(x_min:end),1:N,NFFT,hamming(Nwindow));
                    TFR=abs(tfr);
                    FI=imagescFun(time,freq,10*log10(TFR),[f_min f_max],'Received signal and extracted dispersion curves');
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
                        catch
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

