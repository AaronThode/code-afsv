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

%JULIEN
%%% If beta transform is desired, go ahead and time reverse the signal
%%% before selecting the start time
if var_temp.beta_transform&&var_temp.beta<0
    x_deconv=flipud(x_deconv);
end

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
time=(0:N-1)/Fs;

% Mode selection
while redo
    modes_ok=[];
    mode_index=[];
    Nmode=0;
    for I=1:length(delays)
        x_ok=[x_deconv(delays(I):end); zeros(delays(I)-1,1)];
        [modes,Nm,filtt] = extract_warped_modes(x_ok,Fs,N,params.Nww,r_guess,c1,params.clims,filtt,var_temp.beta_transform,var_temp.beta);
        
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
    
    [mode_index,Isort]=sort(mode_index);
    modes_ok=modes_ok(Isort,:);
    
    
    %Add back original deconvolved source spectrum...
    for I=1:size(modes_ok)
        N=size(modes_ok,2);
        x_f=fft(modes_ok(I,:).',N);
        modes_signal(I,:)=ifft(x_f.*exp(1i*source_phase),N,'symmetric').';

        
    end
    
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
    %%%%% Extract dispersion curves for impulse respons
    
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
    tfr=tfrstft(x_deconv(1:end,hydro_ref),1:N,NFFT,hamming(Nwindow));
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



function  [filtt,params]=get_start_time(x_deconv,NFFT,N,Nwindow,Fs,flims,r_guess,c1,c2,D,Nmodes,filtt,beta_transform,beta)

%%filtt.xmin is final shift..

TFR=20*log10(abs(tfrstft(x_deconv,1:N,NFFT,hamming(Nwindow))));
fig=figure;
imagescFun((1:N)/Fs,Fs*(1:NFFT)/NFFT,TFR,'xy');
ylim([0.8*flims(1) 1.2*flims(2)]);
clear TFR

disp('Time selection')

if beta<0
    title('Signal after source deconvolution, time-reversed');   
    title('Signal after source deconvolution, time-reversed: choose two times (t=r/c, and then time where signal seems to begin)')
    tmp=ginput(2);
    j=max(1,ceil(Fs*tmp(1,1))); % First time instant guess
    params.jmax0=max(1,ceil(Fs*tmp(2,1)));  %index of tmax in time-reversed signal
    %Convert into time-reversed units
    params.dj=params.jmax0-j;
    params.jmax=max(1,params.dj);  %Index is now location of tmax for normal signal.
    
else    
    title('Signal after source deconvolution: choose the first time instant (for warping)')
    tmp=ginput(1);
    params.jmax=[];
end
%caxis([prctile(TFR(:),60) prctile(TFR(:),100)]);
j=max(1,ceil(Fs*tmp(1,1))); % First time instant guess
j=j(1);
dt=j/Fs;
close(fig)

%%% Choose a list of samples near the selected one
frac=1/100;
if beta>=0
    dj=max(1,ceil(j/2*(1-(1-frac)^2)*(1-(c1*j/Fs/r_guess)^2))); % Step (calculated so that the modes shift is about 1/50 of the distance between 2 pekeris cutoff frequencies
    nj=20;                                       % Number of steps in each direction
    j_list=max(1,j-nj*dj):dj:max(1,j+nj*dj);    % List of delay around the selected time

else %beta<0
    dj=max(1,ceil(2*frac*abs(params.dj)));
    nj=20;                                       % Number of steps in each direction
    j_list=max(1,j-nj*dj):dj:max(1,j+nj*dj);    % List of delay around the selected time
                                                % Note this is relative to
                                                % time-reversed signal
    params.jmax_list=params.jmax0-j_list;   %tmax should be a constant                                           

end


% Make a video

%modes_value=zeros(1,Nmodes,length(j_list)); % 2 estimators of the warping quality

% Set the axis scale for the whole video
fig=figure;
params.beta=beta;
params.beta_transform=beta_transform;
params.jmax=params.jmax_list(1);
[params,~]=choose_instant(x_deconv,params,j_list(1),Fs,r_guess,c1,c2,D,Nmodes);

close(fig)

% Starting the video
%mov=figure('units','normalized','outerposition',[0 0 1 1]);
figure

for jj = 1:length(j_list) % try several time instants
    try
    if beta<0
        params.jmax=params.jmax_list(jj);
    end
    [~,~,slice]=choose_instant(x_deconv,params,j_list(jj),Fs,r_guess,c1,c2,D,Nmodes);
    
    if jj==1
        slice_sum=zeros(length(slice),length(j_list));
    end
    slice_sum(:,jj)=slice(1:size(slice_sum,1));
    %modes_value(:,:,jj)=vals;
    title(['delay : ' num2str((j_list(jj)-j)/Fs) ' s'])
    pause(0.05);
    %F(jj) = getframe(mov); % Movie
    catch
       continue 
    end
end


figure;

Fe_max=max(params.fwlims);
delayy=(j_list-j)/Fs;
subplot(2,1,1)
slicedB=10*log10(slice_sum);
pcolor(delayy,params.f_w(1:size(slice_sum,1)),slicedB-max(max(slicedB)));
shading flat
colorbar;colormap(jet);
%caxis([-35 -15])
axis('xy');
ylim([0 Fe_max]);

%%%Kurtosis computation
Fe_window=100;
Irow=1;
II=find(params.f_w<=Fe_window);
II_step=floor(0.25*median(II));
[~,Imax]=min(abs(Fe_max-params.f_w));
while max(II)<=Imax
    kurtosis_val(Irow,:)=kurtosis(slice_sum(II,:));
    kurtosis_fw(Irow)=params.f_w(floor(median(II)));
    
    Irow=Irow+1;
    II=II_step+II;
    %max_window=params.f_w(max(II));
end

kurtosis_val(Irow,:)=kurtosis(slice_sum);
kurtosis_fw(Irow)=Fe_max;
subplot(2,1,2);
pcolor(delayy,kurtosis_fw,10*log10(kurtosis_val))
xlabel('Delay (s)')
ylabel(sprintf('Kurtosis: %6.2f Hz width',Fe_window))
grid on
ylim([0 Fe_max])
colorbar


% We let the user choose the first time instant
colorrange=input('Enter caxis for top window:');
if ~isempty(colorrange)
    subplot(2,1,1)
    caxis(colorrange);
end
Ncount=input('Number of picks:');
tmp=ginput(Ncount);
tmp=tmp(:,1);

%j_min=round(j+str2double(input('\n Choose a time instant: ','s'))*Fs);
%if isnan(j_min)
%    j_min=j_list(opt_delay);
%end
figure;
for I=1:Ncount
    
    subplot(ceil(Ncount/2),2,I)
    j_min=round(j+Fs*tmp(I));
    if beta<0
        params.jmax=params.jmax0-j_min;
    end
    [~,vals]=choose_instant(x_deconv,params,j_min,Fs,r_guess,c1,c2,D,Nmodes);
    %caxis('auto');
    axis('xy')
    %ylim([0 250]);
    ylim([0 max(params.fwlims)]);
    caxis([-20 0]);
    title(sprintf('Choice: %i, Offset: %6.2f',I,tmp(I)));
end
gtext(sprintf('Range guess: %6.2f km, water depth: %6.2f m',r_guess/1000,D))

Ichc=input('Enter choices for selection; e.g [1 3]: ');
j_min=round(j+Fs*tmp(Ichc));

%j_min=round(j+x_min*Fs);

%%% Save the filtering process to apply it to the other hydrophones

filtt.xmin=j_min;
filtt.flims=flims;
filtt.fwlims=params.fwlims;
%N=2^nextpow2(length(x_deconv));
%x_ok=ifft(fft(x_deconv,N).*exp(1i*2*pi*j_min/N*(1:N)'),'symmetric');

%x_ok=[x_deconv(j_min:end); zeros(j_min-1,1)];

%%% end test
end  %get_start_time

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