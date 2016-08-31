function [modes_ok,filtt,modes_signal]=Frequency_warping(x,stft_param,filt_param,hdr,handles,param_file,var_temp)
%function [modes_ok,filtt,modes_signal]=Frequency_warping(x,stft_param,filt_param,hdr,handles,param_file,var_temp)
%  modes_ok: [Nmodes Ntime] of mode time series without source phase
%  modes_signal:  [Nmodes Ntime] of mode time series with source phase
%  added back
%
% var_temp: structure of variables for warping
Nmodes=5;                           % Number of modes to search and localize
D=51;  %Depth for modeling Pekeris frequencies...
%%% Create a folder to save the results
saveFold=[strjoin(strsplit(handles.outputdir,'\'),'/') '/' datestr(handles.tdate_start,30)];

if exist(saveFold,'dir')~=7
    mkdir(saveFold)
end

cd(saveFold)

%%% General parameters and variables
filtt=struct('hdr',hdr); % Save the filtering process to apply it on the other hydrophones

set(0,'DefaultAxesFontName', 'Californian FB')
set(0,'DefaultAxesFontSize', 12)

hydro_ref=str2double(get(handles.edit_chan,'String'));
try
    hydro_ref=hydro_ref-length(find(hdr.geom.Igood(1:hydro_ref)==0)); % In case there are some bad hydrophones
catch
end

next=zeros(1,4)+1;      % A boolean vector used to go back and forward in the program
Nhydros=max([ 1 1*(strcmp(handles.filetype,'GSI'))+15*(strcmp(handles.filetype,'MDAT'))]); % DASAR or VLA
[Fs,NFFT]=deal(stft_param.Fs,stft_param.Nfft);

%%%%Pad x

Nx=length(x(:,1));
N=2^nextpow2(Nx);
buff=round((N-Nx)/2);
y=zeros(N,size(x,2));
y(buff+(1:length(x(:,1))),:)=(x-mean(x)*ones(Nx,1));
x=y;

Nwindow=stft_param.N_window;
[r_guess,c1,c2]=deal(filt_param.r_guess,filt_param.c1,filt_param.c2);
flims=filt_param.flims; % frequency bounds
if min(flims)==0
    uiwait(msgbox('Problem!  Warping requires a lower filtering limit.  Adjust minimum frequency'));
    filtt=[];
    modes_ok=[];
    return
end
%%% Check for bad hydrophones

if Nhydros~=1
    norm_tmp=x/norm(x(:,hydro_ref));
    hydros_ok=var_temp.Nchan(sqrt(sum(norm_tmp.^2,1))>0.01);
else
    hydros_ok=[1]; % DASAR
end

Nhydros=length(hydros_ok);
filtt.Nchan=hydros_ok;

%%%% Signal selection  %%%

for Nh=hydros_ok
    [x(:,Nh),~]=quick_filter(x(:,Nh),Fs,flims(1),flims(2));
end

filtt.Fs=Fs;
freq=(0:NFFT-1)/NFFT*Fs;

%%%% Source deconvolution %%%%%%%


if strcmp(param_file,'')&&var_temp.deconv_chc==1 % If the user didn't select a file for deconvolution
    [x_deconv,source_phase,t_points,f_points] = sourceDeconv_dB(x(:,hydro_ref),Fs,NFFT,Nwindow,flims,filt_param.Ncontour);
    
elseif (var_temp.deconv_chc==1) %otherwise we take the selected points from the file
    load(param_file,'t_points','f_points')
    t=linspace(t_points(1),t_points(end),round(1+Fs*(t_points(end)-t_points(1))));
    iflaw=interp1(t_points,f_points,t,'linear')';
    Fsource=fft(fmodany(iflaw/Fs),N);
    source_phase=angle(Fsource);
    
    disp('End of source deconvolution')
    
    x_f=fft(x(:,hydro_ref),length(source_phase));
    
    x_deconv=ifft(x_f.*exp(-1i*source_phase),'symmetric');
elseif var_temp.deconv_chc==0
    x_deconv=x;
    source_phase=zeros(N,1);
end

next(1:end)=1;

if strcmp(param_file,'')&&var_temp.deconv_chc==1 % If the user didn't select a file for deconvolution
    saveDeconv=input('Do you want to save the deconvolution ?  (y/n) ','s');
    if saveDeconv=='y'
        var_temp.t_points=t_points;
        var_temp.f_points=f_points;
        var_temp.minfreq=flims(1);
        var_temp.maxfreq=flims(2);
        save('param_deconv.mat','-struct','var_temp')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time selection (for warping)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


next(2:end)=1;

[filtt,params]=get_start_time(x_deconv,NFFT,N,Nwindow,Fs,flims,r_guess,c1,c2,D,Nmodes,filtt,var_temp.beta_transform,var_temp.beta);
delays=filtt.xmin;

%%%%filtt.min are times to start extracting...
disp('End of time selection')

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Warping%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

redo=1;
N=length(x_deconv);
%time=(0:N-1)/Fs;

% Mode selection
while redo
    modes_ok=[];
    mode_index=[];
    Nmode=0;
    for I=1:length(delays)
        x_ok=[x_deconv(delays(I):end); zeros(delays(I)-1,1)];
        if var_temp.beta<0
            %nzero=length(x_deconv)-delays(I);
            %x_ok=[x_deconv(1:delays(I)); zeros(nzero,1)];
            x_ok=x_deconv;
         
        end
        [modes,Nm,filtt] = extract_warped_modes(x_ok,Fs,N,params.Nww,r_guess,c1,params.clims,filtt,  ...
           var_temp.beta_transform,var_temp.beta,params.dtmax,params.dt_rc);
        
        temp=input('Enter mode numbers:');
        if isempty(temp)||Nm==0
            continue
        end
        modes=[zeros(Nm,delays(I)-1) modes];
        modes=modes(:,1:N);
        modes_ok=[modes_ok;modes];
        mode_index=[mode_index temp];
        Nmode=Nmode+Nm;
        
    end  %delays
    
    %N=size(modes,1);
    
    [mode_index,Isort]=sort(mode_index);
    modes_ok=modes_ok(Isort,:);
    
    
    %Add back original deconvolved source spectrum...
    for I=1:size(modes_ok)
        N=size(modes_ok,2);
        x_f=fft(modes_ok(I,:).',N);
        modes_signal(I,:)=ifft(x_f.*exp(1i*source_phase),N,'symmetric').';

    end
    time=(0:N-1)/Fs;

%     filtt.mode_stack=zeros(N,Nmode,Nhydros);
%     
%     % Apply for each hydrophone
%     for NN=hydros_ok
%         for mm=1:Nmode
%             M=length(xw_multi(:,NN));
%             rtf=tfrstft(xw_multi(:,NN),1:M,M,hamming(params.Nww));
%             rtf(floor(M/2):end,:)=[];
%             mode_rtf_warp=squeeze(filtt.('mask')(mm,:,:)).*rtf;
%             
%             [mode_temp_warp]=real(sum(mode_rtf_warp,1))*2/M/max(hamming(params.Nww)/norm(hamming(params.Nww)));
%             tmp=iwarp_temp_exa(mode_temp_warp,Fe_w,r_guess,c1,Fs,N);
%             filtt.mode_stack(:,mm,NN)=ifft(fft(tmp).*exp(1i*source_phase),N,1,'symmetric');
%             clear mode_rtf_warp mode_temp_warp
%         end
%     end
%     
%     modes_ok=filtt.mode_stack(:,:,hydro_ref)';
    %%%%% Extract dispersion curves for impulse response
    
    tm=zeros(Nmode,NFFT);
    for m=1:Nmode
        [~, mode_rtfr,~]=tfrrsp(hilbert(modes_ok(m,:).'),1:N,NFFT,hamming(Nwindow));
        [tm(m,:),~]=momftfr(mode_rtfr,1,N,time);
    end 
    leg=repmat(' ',Nmode,7);
    figure
    
    %%%%%Plot impulse response mode estimates
    for m=1:Nmode
        leg_i=['Mode ' num2str(m)];
        leg(m,1:length(leg_i))=leg_i;
        subplot(ceil(sqrt(Nmode)),ceil(sqrt(Nmode)),m)
        tfr_modes=tfrstft(modes_ok(m,:)',1:N,NFFT,hamming(Nwindow));
        imagescFun(time,freq,20*log10(abs(tfr_modes)),'xy')
        ylim(flims)
        title(leg_i)
    end
    
    %%%Plot impulse mode dispersion curve
    if var_temp.beta>=0
        tfr=tfrstft(x_deconv(1:end,hydro_ref),1:N,NFFT,hamming(Nwindow));
    else
        tfr=tfrstft(sum(modes_ok,1)',1:N,NFFT,hamming(Nwindow));
    end
    
    TFR=20*log10(abs(tfr));
    
    figure
    imagescFun(time,freq,TFR,'xy');
    ylim(flims)
    %caxis([prctile(TFR(:),70) prctile(TFR(:),100)])
    title('Impulse response estimate and extracted dispersion curves')
    hold on
    plot(tm,freq,'linewidth',2)
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    grid on
    legend(leg)
    
    %%%%%Plot orignal signal components
    leg=repmat(' ',Nmode,7);
    figure
    for m=1:Nmode
        leg_i=['Mode ' num2str(m)];
        leg(m,1:length(leg_i))=leg_i;
        subplot(ceil(sqrt(Nmode)),ceil(sqrt(Nmode)),m)
        tfr_modes=tfrstft(modes_signal(m,:)',1:N,NFFT,hamming(Nwindow));
        imagescFun(time,freq,20*log10(abs(tfr_modes)),'xy')
        ylim(flims)
        title(leg_i)
    end
    
    %subplot(ceil(sqrt(Nmode)),ceil(sqrt(Nmode)),m+1)
    figure(300)
    tfr_modes=tfrstft(sum(modes_signal,1)',1:N,NFFT,hamming(Nwindow));
    imagescFun(time,freq,20*log10(abs(tfr_modes)),'xy')
    ylim(flims)
    title('All modes')
 
    tm=zeros(Nmode,NFFT);
    for m=1:Nmode
        [~, mode_rtfr,~]=tfrrsp(hilbert(modes_signal(m,:).'),1:N,NFFT,hamming(Nwindow));
        [tm(m,:),~]=momftfr(mode_rtfr,1,N,time);
    end
    
     %%%Plot original mode dispersion curve
    tfr=tfrstft(x(1:end,hydro_ref),1:N,NFFT,hamming(Nwindow));
    TFR=20*log10(abs(tfr));
    
    figure
    imagescFun(time,freq,TFR,'xy');
     ylim(flims)
    %caxis([prctile(TFR(:),70) prctile(TFR(:),100)])
    title('Received signal and extracted dispersion curves')
    hold on
    plot(tm,freq,'linewidth',2)
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    grid on
    legend(leg)
    redo=input('Enter 1 to redo mode selection:');
    if isempty(redo)
        redo=0;
    end
end %while redo

filtt.r_guess=r_guess;
filtt.c1=c1;
filtt.D=D;

%end

end




function writerObj=make_movie(writerObj,s_w,M,t_w,f_w,params,Nmodes,c1,c2,D)
RTF=abs(tfrstft(s_w,1:M,M,hamming(params.Nww)));
M1=floor(M/2);
RTF(M1:end,:)=[];                   % Delete negative frequencies

fig=figure;
colorbar
imagescFun(t_w,f_w,20*log10((RTF)),'ij');

for mm=1:Nmodes
    % Pekeris cutoff frequencies
    
    pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
    colorMode=['r','w','c','y','g','m'];
    hold on
    plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
end

xlim([0 params.xmax])
ylim(params.fwlims)
%caxis(params.clims)
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

imagescFun(t_w,f_w,20*log10(abs(RTF)),'ij');
ylim(params.fwlims)
%caxis(params.clims)
%colorbar
xlim([0 params.xmax])
hold on
try
    title(['Hydro n°' num2str(hydros_ok(NN)) ' (' num2str(hdr.geom.rd(hydros_ok(NN))) 'm)'])
end

for mm=1:Nmodes
    % Pekeris cutoff frequencies
    D=51;
    pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
    
    colorMode=['r','w','c','y','g','m'];
    plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
end
end