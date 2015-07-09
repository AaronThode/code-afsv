function [filt]=Whale_main_fun(x,stft_param,filt_param,hdr,handles,param_file,var_temp)

%% Create a folder to save the results
saveFold=[strjoin(strsplit(handles.outputdir,'\'),'/') '/' datestr(handles.tdate_start,30)];

if exist(saveFold,'dir')~=7
    mkdir(saveFold)
end

cd(saveFold)

%% General parameters and variables
filt=struct('hdr',hdr); % Save the filtering process to apply it on the other hydrophones

set(0,'DefaultAxesFontName', 'Californian FB')
set(0,'DefaultAxesFontSize', 12)

hydro_ref=str2double(get(handles.edit_chan,'String'));
try
    hydro_ref=hydro_ref-length(find(hdr.geom.Igood(1:hydro_ref)==0)); % In case there are some bad hydrophones
catch
end

next=zeros(1,4)+1;      % A boolean vector used to go back and forward in the program
Nhydros=1*(strcmp(handles.filetype,'GSI'))+15*(strcmp(handles.filetype,'MDAT')); % DASAR or VLA
[Fs,NFFT]=deal(stft_param.Fs,stft_param.Nfft);
N=length(x(:,1));
Nwindow=stft_param.N_window;
[r_guess,c1,c2]=deal(filt_param.r_guess,filt_param.c1,filt_param.c2);
flims=filt_param.flims; % frequency bounds
if min(flims)==0
    uiwait(msgbox('Problem!  Warping requires a lower filtering limit.  Adjust minimum frequency'));
    filt=[];
    return
end
%% Check for bad hydrophones

if Nhydros~=1
    norm_tmp=x/norm(x(:,hydro_ref));
    hydros_ok=var_temp.Nchan(sqrt(sum(norm_tmp.^2,1))>0.01);
else
    hydros_ok=[1]; % DASAR
end

Nhydros=length(hydros_ok);
filt.Nchan=hydros_ok;

%% Signal selection

for Nh=hydros_ok
    [x(:,Nh),~]=quick_filter(x(:,Nh),Fs,flims(1),flims(2));
end

filt.Fs=Fs;
freq=(0:NFFT-1)/NFFT*Fs;

%% Source deconvolution

while next(1)
    if strcmp(param_file,'') % If the user didn't select a file for deconvolution
        [x_deconv,source_phase,t_points,f_points] = sourceDeconv_dB(x(:,hydro_ref),Fs,NFFT,Nwindow,flims,filt_param.Ncontour);
        
    else %otherwise we take the selected points from the file
        load(param_file,'t_points','f_points')
        t=linspace(t_points(1),t_points(end),round(1+Fs*(t_points(end)-t_points(1))));
        iflaw=interp1(t_points,f_points,t,'linear')';
        Fsource=fft(fmodany(iflaw/Fs),N);
        source_phase=angle(Fsource);
        
        disp('End of source deconvolution')
        
        x_f=fft(x(:,hydro_ref),length(source_phase));
        
        x_deconv=ifft(x_f.*exp(-1i*source_phase),'symmetric');
    end
    
    next(1:end)=1;
    
    %% Time selection (for warping)
    
    while next(2)
        next(2:end)=1;
        
        TFR=10*log10(abs(tfrstft(x_deconv,1:N,NFFT,hamming(Nwindow))));
        fig=figure;
        imagescFun((1:N)/Fs,Fs*(1:NFFT)/NFFT,TFR,'xy');
        ylim([0.8*flims(1) 1.2*flims(2)])
        title('Signal after source deconvolution')
        caxis([prctile(TFR(:),60) prctile(TFR(:),100)]);
        clear TFR
        
        if strcmp(param_file,'') % If the user didn't select a file for deconvolution
            saveDeconv=input('Do you want to save the deconvolution ?  (y/n) ','s');
            if saveDeconv=='y'
                var_temp.t_points=t_points;
                var_temp.f_points=f_points;
                var_temp.minfreq=flims(1);
                var_temp.maxfreq=flims(2);
                save('param_deconv.mat','-struct','var_temp')
            end
        end
        
        disp('Time selection')
        title('Signal after source deconvolution: choose the first time instant (for warping)')
        j=max(1,ceil(Fs*ginput(1))); % First time instant guess
        j=j(1);
        close(fig)
        
        %% Choose a list of samples near the selected one
        frac=1/100;
        dj=max(1,ceil(j/2*(1-(1-frac)^2)*(1-(c1*j/Fs/r_guess)^2))); % Step (calculated so that the modes shift is about 1/50 of the distance between 2 pekeris cutoff frequencies
        nj=5;                                       % Number of steps in each direction
        j_list=max(1,j-nj*dj):dj:max(1,j+nj*dj);    % List of delay around the selected one
        
        % Make a video
        
        Nmodes=3;                           % Number of modes to search and localize
        modes_value=zeros(3,Nmodes,length(j_list)); % 2 estimators of the warping quality
        
        % Set the axis scale for the whole video
        fig=figure;
        [params,~]=choose_instant(x_deconv,struct(),j_list(1),Fs,r_guess,c1,c2,Nmodes);
        
        close(fig)
        
        % Starting the video
        mov=figure('units','normalized','outerposition',[0 0 1 1]);
        
        for jj = 1:length(j_list) % try several time instants
            [~,vals]=choose_instant(x_deconv,params,j_list(jj),Fs,r_guess,c1,c2,Nmodes);
            
            modes_value(:,:,jj)=vals;
            title(['delay : ' num2str((j_list(jj)-j)/Fs) ' s'])
            F(jj) = getframe(mov); % Movie
        end
        
        commandwindow;
        replay=input('\n Replay ? (y/n) ','s');
        
        while(replay~='n')
            shg
            movie(mov,F,1,1) % Replay the movie once with 1 image/sec
            commandwindow;
            replay=input('Replay ? (y/n) ','s');
        end
        
        f1_opt=squeeze(prod(squeeze(modes_value(1,:,:)),1)); % First optimisation function (Mode energy)
        f2_opt=squeeze(prod(squeeze(modes_value(2,:,:)),1)); % Second optimisation function (Slope)
        
        figure,hold on
        plot((j_list-j)/Fs,f1_opt/max(f1_opt))
        plot((j_list-j)/Fs,1-f2_opt/max(f2_opt))
        xlabel('Delay (s)')
        legend('Modes energy','Modes sparsity')
        title('Warping efficiency')
        
        figure,hold on
        for mm=1:Nmodes
            plot((j_list-j)/Fs,squeeze(modes_value(3,mm,:)))
            xlabel('Delay (s)')
            ylabel('Slope')
            title('Warping efficiency')
        end
        legend(cellstr(num2str((1:Nmodes)', 'Mode %i')))
        
        [~,opt_delay]=min(sum(squeeze(modes_value(3,:,:)),1));
        
        disp(['The delay that optimize the slope is :' num2str((j_list(opt_delay)-j)/Fs)])
        
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
        x_min=round(j+str2double(input('\n Choose a time instant: ','s'))*Fs);
        
        %% Save the filtering process to apply it to the other hydrophones
        
        filt.xmin=x_min;
        filt.flims=flims;
        filt.fwlims=params.fwlims;
        
        x_ok=ifft(fft(x_deconv,N).*exp(1i*2*pi*x_min/N*(1:N)'),'symmetric');
        
        %% end test
        
        disp('End of time selection')
        
        % Apply for all the hydrophones
        
        writerObj = VideoWriter('warping_Nhydros.avi'); % Show all the warping in a movie
        writerObj.FrameRate=5;
        open(writerObj);
        
        %if Nhydros>1
        
        for NN=1:Nhydros
            % Source phase removal
            x_multi=ifft(fft(x(:,hydros_ok(NN)),N).*exp(1i*(-source_phase+2*pi*x_min*(1:N)'/N)),N,'symmetric');
            % Warping
            [s_w, Fe_w]=warp_temp_exa(x_multi,Fs,r_guess,c1);     % s_w: warped signal, Fe_w: new warping frequency
            
            if NN==1
                M=length(s_w);
                xw_multi=zeros(M,Nhydros);
            end
            
            xw_multi(:,hydros_ok(NN))=s_w;
            
            t_w=(0:M-1)/Fe_w;        % Warped time
            f_w=(0:M-1)*Fe_w/M;      % Warped frequencies
            
            RTF=abs(tfrstft(s_w,1:M,M,hamming(params.Nww)));
            M1=floor(M/2);
            RTF(M1:end,:)=[];                   % Delete negative frequencies
            
            fig=figure;
            colorbar
            imagescFun(t_w,f_w,RTF,'ij');
            
            for mm=1:Nmodes
                % Pekeris cutoff frequencies
                D=55;
                pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
                colorMode=['r','w','c','y','g','m'];
                hold on
                plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
            end
            
            xlim([0 params.xmax])
            ylim(params.fwlims)
            caxis(params.clims)
            try
                title(['Hydro n°' num2str(hydros_ok(NN)) ' (' num2str(hdr.geom.rd(hydros_ok(NN))) 'm)'])
            end
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
            ylim(params.fwlims)
            caxis(params.clims)
            colorbar
            xlim([0 params.xmax])
            hold on
            try
                title(['Hydro n°' num2str(hydros_ok(NN)) ' (' num2str(hdr.geom.rd(hydros_ok(NN))) 'm)'])
            end
            
            for mm=1:Nmodes
                % Pekeris cutoff frequencies
                D=55;
                pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
                
                colorMode=['r','w','c','y','g','m'];
                plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
            end
            
        end
        close(writerObj);
        
        %end
        
        %% Warping
        
        N=length(x_ok);
        time=(0:N-1)/Fs;
        
        % Mode selection
        [~,Nmode,choice,filt] = mode_select3(x_ok,Fs,N,params.Nww,r_guess,c1,params.clims,filt);
        
        if choice<3
            next(choice+1:end)=0;
        end
        next(1:choice)=1;
        
        while next(3)
            %% Mode extraction & source reconvolution (add the initial source phase that was removed)
            
            if Nmode~=0
                filt.mode_stack=zeros(N,Nmode,Nhydros);
                
                % Apply for each hydrophone
                for NN=hydros_ok
                    for mm=1:Nmode
                        M=length(xw_multi(:,NN));
                        rtf=tfrstft(xw_multi(:,NN),1:M,M,hamming(params.Nww));
                        rtf(floor(M/2):end,:)=[];
                        mode_rtf_warp=squeeze(filt.('mask')(mm,:,:)).*rtf;
                        
                        [mode_temp_warp]=real(sum(mode_rtf_warp,1))*2/M/max(hamming(params.Nww)/norm(hamming(params.Nww)));
                        tmp=iwarp_temp_exa(mode_temp_warp,Fe_w,r_guess,c1,Fs,N);
                        filt.mode_stack(:,mm,NN)=ifft(fft(tmp).*exp(1i*source_phase),N,1,'symmetric');
                        clear mode_rtf_warp mode_temp_warp
                    end
                end
                
                modes_ok=filt.mode_stack(:,:,hydro_ref)';
                
                %% Extraction dispersion curves
                
                tm=zeros(Nmode,NFFT);
                
                for m=1:Nmode
                    [~, mode_rtfr,~]=tfrrsp(hilbert(modes_ok(m,:).'),1:N,NFFT,hamming(Nwindow));
                    [tm(m,:),~]=momftfr(mode_rtfr,1,N,time);
                end
                
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
                tfr=tfrstft(x(x_min:end,hydro_ref),1:N,NFFT,hamming(Nwindow));
                TFR=10*log10(abs(tfr));
                figure,
                imagescFun(time,freq,TFR,'xy');
                ylim(flims)
                caxis([prctile(TFR(:),70) prctile(TFR(:),100)])
                title('Received signal and extracted dispersion curves')
                hold on
                plot(tm,freq,'linewidth',2)
                xlabel('Time [s]')
                ylabel('Frequency [Hz]')
                grid on
                legend(leg)
                
                temp=str2double(input(['\nAre you satisfied with the mode selection ? \n' ...
                    '1: Back to source deconvolution\n'...
                    '2: Back to time selection (for warping)\n'...
                    '3: Save and continue\n'...
                    '4: Continue\n'],'s'));
                
                if temp==3
                    save('filtered_modes.mat','filt','c1','stft_param','filt_param')
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

