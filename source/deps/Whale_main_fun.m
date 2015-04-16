function [filt]=Whale_main_fun(x,stft_param,filt_param,hdr)

filt=struct('mask',[]); % Save the filtering process to apply it on the other hydrophones

set(0,'DefaultAxesFontName', 'Californian FB')
set(0,'DefaultAxesFontSize', 12)

% Overall parameters that may have to be changed depending on the whale call.

next=zeros(1,4)+1;      % A boolean vector used to go back and forward in the program
Nhydros=1*(hdr.Nchan==3)+15*(hdr.Nchan==15); % DASAR or VLA
[Fs,NFFT]=deal(stft_param.Fs,stft_param.Nfft);
[Nwindow]=deal(stft_param.N_window);
[r_guess,c1]=deal(filt_param.r_guess,filt_param.c1);
flims=filt_param.flims;

%% Signal selection

for Nh=1:Nhydros
    [x(:,Nh),~]=quick_filter(x(:,Nh),Fs,flims(1),flims(2));
end

filt.Fs=Fs;
freq=(0:NFFT-1)/NFFT*Fs;
    
%% Source deconvolution

while next(1)
    if Nhydros==1 % DASAR
        [x_deconv,source] = sourceDeconv_dB(x(:,1),Fs,NFFT,Nwindow,flims,filt_param.Ncontour);
    else % VLA
        [x_deconv,source] = sourceDeconv_dB(x(:,9),Fs,NFFT,Nwindow,flims,filt_param.Ncontour);
    end
        
    next(1:end)=1;
    filt.source=source;
    Nx=length(x_deconv);
    
%% Time selection (for warping)

    while next(2)
        next(2:end)=1;

        disp('Time selection')
        TFR=abs(tfrstft(x_deconv,1:Nx,NFFT,hamming(Nwindow)));
        fig=figure;
        imagescFun((1:Nx)/Fs,Fs*(1:NFFT)/NFFT,10*log10(TFR),'xy');
        ylim(flims)
        title('Signal after source deconvolution: choose the first time instant (for warping)')
        caxis([prctile(prctile(10*log10(TFR),50),50) prctile(prctile(10*log10(TFR),100),100)]);
        clear TFR

        j=ceil(Fs*ginput(1));
        j=j(1);
        close(fig)

 %% Choose a list of samples near the selected one

        dj=ceil(j/2*(1-(49/50)^2)*(1-(c1*j/Fs/r_guess)^2)); % Step (calculated so that the modes shift is about 1/50 of the distance between 2 pekeris cutoff frequencies
        nj=5;                                       % Number of steps in each direction
        j_list=max(1,j-nj*dj):dj:max(1,j+nj*dj);    % List of delay around the selected one

        % Make a video

        Nmodes=3;                           % Number of modes to search and localize
        modes_value=zeros(2,Nmodes,2*nj+1); % 2 estimators of the warping quality

        % Set the scale for the axis for the whole video
        fig=figure;
        [fwlims,clims,xmax,vals]=choose_instant(x_deconv,j_list(1),Fs,r_guess,c1,Nmodes,max(j_list),1);
        modes_value(:,:,1)=vals;      
        close(fig)

        % Starting the video
        mov=figure('units','normalized','outerposition',[0 0 1 1]);

        for jj = 1:length(j_list)
            [~,~,~,vals]=choose_instant(x_deconv,j_list(jj),Fs,r_guess,c1,Nmodes,xmax,fwlims,clims);
            modes_value(:,:,jj)=vals;
            title(['delay : ' num2str((j_list(jj)-j)/Fs) ' s'])
            F(jj) = getframe(mov); % Movie  
        end

        f1_opt=squeeze(prod(squeeze(modes_value(1,:,:)),1));
        f2_opt=squeeze(prod(squeeze(modes_value(2,:,:)),1));

        figure,plot((j_list-j)/Fs,f1_opt/max(f1_opt))
        hold on
        plot((j_list-j)/Fs,1-f2_opt/max(f2_opt))
        xlabel('Delay (s)')
        legend('Modes energy','Modes sparsity')
        title('Warping efficiency')

        [~,opt_delay]=max(f1_opt);

        disp(['The delay that optimize the warping is :' num2str((j_list(opt_delay)-j)/Fs)]) 

        commandwindow;
        replay=input('\n Replay ? (y/n) ','s');

        while(replay~='n')
            shg
            movie(mov,F,1,1) % Replay the movie once with 1 image/sec
            commandwindow;
            replay=input('Replay ? (y/n) ','s');
        end

        commandwindow;
        saveMov=input('\n Do you want to save the movie ? (y/n) ','s');

        if saveMov == 'y' % Saving the movie
            writerObj = VideoWriter('warping.avi');
            writerObj.FrameRate=5;
            open(writerObj);
            for jj=1:length(j_list)
                writeVideo(writerObj,F(jj));
            end
            close(writerObj);
        end

        clear F x_multi modes_value writerObj
        
        % We let the user choose the first time instant
        x_min=floor(j+str2double(input('\n Choose a time instant: ','s'))*Fs);
        
        % Keep the filtering process to apply it to the other hydrophones
        filt.xlims=[x_min,length(x_deconv)];
        filt.flims=flims;
        filt.fwlims=fwlims;
        
        x_ok=x_deconv(x_min:end);
        disp('End of time selection')

        % Apply for all the hydrophones
        writerObj = VideoWriter('warping_Nhydros.avi'); % Show all the warping in a movie
        writerObj.FrameRate=5;
        open(writerObj);

        if Nhydros==15

        for NN=1:Nhydros
            x_multi=Whale_multi_hydro(filt,x(:,NN));

            [s_w, Fe_w]=warp_temp_exa(x_multi,Fs,r_guess,c1);     % s_w: warped signal, Fe_w: new warping frequency
            M=length(s_w);

            t_w=(0:M-1)/Fe_w;                           % Warped time
            f_w=(0:M-1)*Fe_w/M;                         % Warped frequencies
            Nww=1+2*ceil(M/8);                          % Number of points for the tfr of the warped signal

            RTF=abs(tfrstft(s_w,1:M,M,hamming(Nww))).^2; % RTF is squared to better see the modes
            M1=floor(M/2);
            RTF(M1:end,:)=[];                   % Delete negative frequencies
            
            if NN==1; colorb=[0 1.2*max(RTF(:))]; end
            
            fig=figure;
            imagescFun(t_w,f_w,RTF,'ij');
           
            for mm=1:Nmodes
            % Pekeris cutoff frequencies
            [c1,c2]=deal(1439,1673);
            D=55;
            pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);

            colorMode=['r','w','c','y','g','m'];
            hold on
            plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
            end
            
            xlim([0 xmax])
            ylim(fwlims)
            caxis(colorb)
            colorbar
            title(['Hydro n°' num2str(NN) ' (' num2str(hdr.geom.rd(NN)) 'm)'])
            
            F(jj) = getframe(fig);
            writeVideo(writerObj,F(jj));
            close(fig)
            
            % Plot each warping in 4 2x2 subplots
            
            switch NN
                case{1,2,3,4}
                    subplot(2,2,NN)
                case{5}
                    figure,
                    subplot(2,2,1)
                case{6,7,8}
                    subplot(2,2,NN-4)
                case{9}
                    figure,
                    subplot(2,2,1)
                case{10,11,12}
                    subplot(2,2,NN-8)
                case{13}
                    figure,
                    subplot(2,2,1)
                case{14,15}
                    subplot(2,2,NN-12)
                otherwise
                    disp('error : wrong number of hydrophones')
            end
            
        imagescFun(t_w,f_w,RTF,'ij');
        ylim(fwlims)
        caxis(colorb)
        colorbar
        
        xlim([0 xmax])
        hold on
        title(['Hydro n°' num2str(NN) ' (' num2str(hdr.geom.rd(NN)) 'm)'])

        for mm=1:Nmodes
            % Pekeris cutoff frequencies
            [c1,c2]=deal(1439,1673);
            D=55;
            pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);

            colorMode=['r','w','c','y','g','m'];
            plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
        end
        
        end
        close(writerObj);
        
        end

    %% Warping

        N=length(x_ok);
        time=(0:N-1)/Fs;

        [modes,Nmode,choice,filt] = modeSelect(x_ok,Fs,N,r_guess,c1,'samp',filt);

            if choice<3
                next(choice+1:end)=0;
            end
            next(1:choice)=1;
            
            while next(3)
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
                
                for m=1:Nmode
                    leg_i=['Mode ' num2str(m)];
                    leg(m,1:length(leg_i))=leg_i;
                    subplot(ceil(sqrt(Nmode)),ceil(sqrt(Nmode)),m)
                    tfr_modes=tfrstft(modes_ok(m,:)',1:N,NFFT,hamming(Nwindow));
                    imagescFun(time,freq,abs(tfr_modes),'xy')
                    ylim(flims)
                    title(leg_i)
                end

                % Plot received TFR+dispersion curves
                tfr=tfrstft(x(x_min:end,1),1:N,NFFT,hamming(Nwindow));
                TFR=abs(tfr);
                FI=figure;
                imagescFun(time,freq,10*log10(TFR),'xy');
                ylim(flims)
                caxis([prctile(prctile(10*log10(TFR),55),55) prctile(prctile(10*log10(TFR),100),100)])
                title('Received signal and extracted dispersion curves')
                hold on
                plot(tm,freq,'linewidth',2)
                xlabel('Time [s]')
                ylabel('Frequency [Hz]')
                legend(leg)

                temp=str2double(input(['\nAre you satisfied with the mode selection ? \n' ...
                '1: Back to source deconvolution\n'...
                '2: Back to time selection (for warping)\n'...
                '3: Save and exit\n'...
                '4: Exit\n'],'s'));

                if temp==3
                    save('filtered_modes.mat')
                    try
                        saveas(FI,'FI')
                    catch
                        disp('The figure was closed before it was saved')
                    end
                end

                if temp<3
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

