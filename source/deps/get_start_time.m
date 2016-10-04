%function  [filtt,params]=get_start_time(x_deconv,NFFT,N,Nwindow,Fs,flims,r_guess,c1,c2,D,Nmodes,filtt,beta_transform,beta)
function  [filtt,params]=get_start_time(x_deconv,NFFT,N,Nwindow,Fs,flims,r_guess,c1,c2,D,Nmodes,filtt,beta_transform,beta)

%%filtt.xmin is final shift..

TFR=20*log10(abs(tfrstft(x_deconv,1:N,NFFT,hamming(Nwindow))));
fig=figure(5);
imagescFun((1:N)/Fs,Fs*(1:NFFT)/NFFT,TFR,'xy');
ylim([0.8*flims(1) 1.2*flims(2)]);
clear TFR

disp('Time selection')
%%%%%Julien, here one choses two points that allow t=r/c to be estimated.
if beta~=0
    
    nothappy=1;
    while ~isempty(nothappy)
        
        disp('Select two points to estimate r/c estimation, last point is tmax:');
        tmp=ginput(2);
        f=tmp(:,2);
        tt=tmp(:,1);
        rr=(f(2)/f(1)).^((1+beta)/beta);
        Trc=(rr*tt(2)-tt(1))/(rr-1);
        fprintf('Estimated rc time is %8.6f seconds\n',Trc);
        hold on
        % line([1 1]*Trc,[0.8*flims(1) 1.2*flims(2)],'linewidth',2,'color','w')
        line([1 1]*Trc,[0.8*flims(1) 1.2*flims(2)],'linewidth',3,'color','w','linestyle','--')
        line([1 1]*tt(2),[0.8*flims(1) 1.2*flims(2)],'linewidth',2,'color','g')
        
        nothappy=input('Not happy?');
    end
    
    j=max(1,ceil((Fs*Trc))); % r/c estimate
    params.jmax=max(1,ceil(Fs*tt(2)));  %index of tmax
    params.j=j;
    %Convert into time-reversed units
    params.dj=max(1,ceil(diff(tt)*Fs));
    
else
    %     title('Signal after source deconvolution: choose the first time instant (for warping)')
    %     tmp=ginput(1);
    %     params.jmax=[];
    %     j=max(1,ceil(Fs*tmp(1,1))); % First time instant guess
    %
    nothappy=1;
    while ~isempty(nothappy)
        beta=1;
        disp('Select two points to estimate r/c estimation, last point is tmax:');
        tmp=ginput(2);
        f=tmp(:,2);
        tt=tmp(:,1);
        rr=(f(2)/f(1)).^((1+beta)/beta);
        Trc=(rr*tt(2)-tt(1))/(rr-1);
        fprintf('Estimated rc time is %8.6f seconds\n',Trc);
        hold on
        line([1 1]*Trc,[0.8*flims(1) 1.2*flims(2)],'linewidth',3,'color','w','linestyle','--')
        line([1 1]*tt(2),[0.8*flims(1) 1.2*flims(2)],'linewidth',2,'color','g')
        
        nothappy=input('Not happy?');
    end
    
    j=max(1,ceil((Fs*Trc))); % r/c estimate
    params.jmax=max(1,ceil(Fs*tt(2)));  %index of tmax
    params.j=j;
    %Convert into time-reversed units
    params.dj=max(1,ceil(diff(tt)*Fs));
    
end
%%%End r/c and initial start time selection

%caxis([prctile(TFR(:),60) prctile(TFR(:),100)]);
j=j(1);
dt=j/Fs;
%close(fig)

%%% Choose a list of time samples surrounding the selected one
frac=1/100;
if beta>=0
    dj=max(1,ceil(j/2*(1-(1-frac)^2)*(1-(c1*j/Fs/r_guess)^2))); % Step (calculated so that the modes shift is about 1/50 of the distance between 2 pekeris cutoff frequencies
    nj=20;                                       % Number of steps in each direction
    j_list=max(1,j-nj*dj):dj:max(1,j+nj*dj);    % List of delay around the selected time
    
else %beta<0
    
    %%%Julien, I found I need a larger spread in possible steps for
    %%%refractive modes
    dj=max(1,ceil(10*frac*abs(params.dj)));
    nj=10;                                       % Number of steps in each direction
    %j_list=max(1,j:dj:max(1,j+2*nj*dj));    % List of delay around the selected time
    j_list=max(1,j-nj*dj):dj:max(1,j+nj*dj);    % List of delay around the selected time
    
    %params.jmax_list=params.jmax0-j_list;   %tmax should be a constant
    j_list=j_list(j_list>params.jmax&j_list<length(x_deconv));
end
params.beta=beta;
params.beta_transform=beta_transform;
%params.jmax=params.jmax_list(1);


% Set the axis scale for the whole video
fig=figure;
[params,~]=choose_instant(x_deconv,params,j_list(1),Fs,r_guess,c1,c2,D,Nmodes);
%filtt.t_w_min=params.t_w_min;

close(fig)

% Starting the video
%mov=figure('units','normalized','outerposition',[0 0 1 1]);
figure

for jj = 1:length(j_list) % try several time instants
    try
        [~,~,slice]=choose_instant(x_deconv,params,j_list(jj),Fs,r_guess,c1,c2,D,Nmodes);
        
        if jj==1
            slice_sum=zeros(length(slice),length(j_list));
        end
        title(sprintf('delay : %6.2f s', ((j_list(jj)-j)/Fs)));
        pause(0.05);
        
        Np=min(size(slice_sum,1),length(slice));
        slice_sum(1:Np,jj)=slice(1:Np);
        %modes_value(:,:,jj)=vals;
        
        %F(jj) = getframe(mov); % Movie
    catch
        disp('crash')
        continue
    end
end


selects=make_kurtosis;

figure;
for I=1:Ncount
    
    subplot(ceil(Ncount/2),2,I)
    j_min=round(j+Fs*selects(I));
    [params_opt{I},~]=choose_instant(x_deconv,params,j_min,Fs,r_guess,c1,c2,D,Nmodes);
    %caxis('auto');
    axis('xy')
    %ylim([0 250]);
    ylim([0 max(params.fwlims)]);
    caxis([-20 0]);
    title(sprintf('Choice: %i, Offset: %6.2f',I,selects(I)));
end
gtext(sprintf('Range guess: %6.2f km, water depth: %6.2f m',r_guess/1000,D))


Ichc=input('Enter choices for selection; e.g [1 3]: ');
j_min=round(j+Fs*selects(Ichc));
[params,~]=choose_instant(x_deconv,params_opt{Ichc},j_min,Fs,r_guess,c1,c2,D,Nmodes);

%j_min=round(j+x_min*Fs);

%%% Save the filtering process to apply it to the other hydrophones

filtt.xmin=j_min;
filtt.flims=flims;
filtt.fwlims=params.fwlims;


%%%%%%%%subroutine make_kurtosis
    function selects=make_kurtosis
        figure(6);
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
        try
            pcolor(delayy,kurtosis_fw,10*log10(kurtosis_val))
        end
        xlabel('Delay (s)')
        ylabel(sprintf('Kurtosis: %6.2f Hz width',Fe_window))
        %grid on
        ylim([0 Fe_max]);
        colorbar
        
        
        % We let the user choose the first time instant
        %colorrange=input('Enter caxis for top window:');
        colorrange=[-10 0];
        if ~isempty(colorrange)&&length(colorrange)==2
            subplot(2,1,1)
            caxis(colorrange);
        end
        Ncount=input('Number of picks:');
        tmp=ginput(Ncount);
        selects=tmp(:,1);
        
    end
end  %get_start_time
