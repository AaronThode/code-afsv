function	[metrics,hf]	=	shipdet_metrics(PSD_struct, param, ha_ref, hf,want_hough)
%function	[metrics,hf]	=	shipdet_metrics(PSD_struct, param, ha_ref, hf,want_hough)
% Ship detection metrics, aggregate function
% Each metric is calculated in its own subfunction,
% and plotted in its own subplot
%
% Jit Sarkar and Aaron Thode
% MPL/SIO
%
% Function form:
%	hf	=	shipdet_metrics(PSD_struct, ha_ref,hf,want_hough)
%
% Inputs:
%	PSD_struct	:	{struct}	Power spectral density structure,
%								must contain the following	fields:
%		PSD		:	{dbl}(t,f)	The 2d power spectral density matrix, linear (not dB) units
%		T		:	{dbl}(t)	Time vector
%		F		:	{dbl}(f)	Frequency vector
%       Idebug:                 scalar that determines if debug output given...
%
%    param:
%                               threshold=5; %dB threshold between peak and surrounding values at +/-df_search
%                               df_search=10; %+-/Hz to examine around each local maximum
%                               f_min	=	25; %Hz
%                               f_max	=	450; %Hz
%
%	ha_ref		:	{dbl}(1)	Axes handle to Specgram plot in main window
%	hf			:	{dbl}(2)	Figure handle to which plots should be made
%   want_hough:  If exists, do a Hough transform analysis...
% Outputs:
%	hf			:	{dbl}(2)	handle to figure windows
%   metrics     :   {struct}    contains .entropy, .kurtosis .avg.peaks, and .med.peaks fields
%       .entropy.med:  entropy of median linear PSD
%       .entropy.avg:  entroyy of averaged linear PSD
%       .kurtosis.med
%       .kurtosis.avg
%       .kurtosis.avgdB:
%       .avg.peaks:  peak pick applied to average of dB PSD of spectrogram
%               
%       .med.peaks:  peak pick applied to median of dB PSD of spectrogram
%
%       each peaks field contains the following (see peak_picker_Thode.m
%                       for details):
%               isi, for "ISI_tone"
%               adp, for "adaptive"

%	Parameter checking
if ~exist('hf','var')
    hf(1)	=	-1;
    hf(2)	=	-1;
elseif isempty(hf)
    hf(1)=20;
    hf(2)=21;
end

if isempty(param)
    
    threshold=5; %dB
    df_search=10; %+-/Hz to examine
    f_min	=	25;
    f_max	=	450;
else
    threshold=param.threshold; %dB
    df_search=param.df_search; %+-/Hz to examine
    f_min	=	param.f_min;
    f_max	=	param.f_max;
end

for	ii	=	1:length(hf);
    if ishandle(hf(ii))
        clf(hf(ii));
    elseif (hf(ii)>0)
        % 		hf(ii)	=	figure('Name', 'Detection metrics',...
        % 				'Units','normalized','Position',[0.5*(ii-1), 0.0, 0.5, 1.0]);
        hf(ii)	=	figure('Name', 'Detection metrics',...
            'Units','normalized');
    end
end


PSD	=	PSD_struct.PSD;
T	=	PSD_struct.T;
F	=	PSD_struct.F;
%Fs	=	PSD_struct.Fs;
%Nfft=	PSD_struct.Nfft;


Nt	=	length(T);
Nf	=	length(F);
dF	=	mean(diff(F));
%dT	=	mean(diff(T));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Limited frequency range for statistical operations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f_min	=	50;
%f_max	=	450;
[~,if_min]	=	min(abs(F - f_min));
[~,if_max]	=	min(abs(f_max - F));
if_min	=	if_min(1);
if_max	=	if_max(1);

PSD_avg		=	sum(PSD,2)/Nt;
PSD_med     =   median(PSD')';
P_tot		=	sum(PSD_avg)*dF;

PSD_dB		=	abs(10*log10(abs(PSD)));


%%Entropy and kurtosis  over time on both dB and linear data...
Pt	=	sum(PSD(if_min:if_max,:),1)/Nf;
PtdB	=	sum(PSD_dB(if_min:if_max,:),1)/Nf;
[EtdB, EmaxtdB]	=	sentropy(PtdB);
KtdB			=	kurtosis(PtdB);
[Et, Emaxt]	=	sentropy(Pt);
Kt			=	kurtosis(Pt);


%	Entropy and kurtosis, median value across window..
PfdB_avg	=	sum(PSD_dB,2)/Nt;
PfdB_med= median(PSD_dB')';

Pf_avg	=	sum(PSD,2)/Nt;
Pf_med= median(PSD')';

[Ef_med, Emaxf_med]	=	sentropy(Pf_med(if_min:if_max));
Kf_med			=	kurtosis(Pf_med(if_min:if_max));

[Ef_avg, Emaxf_avg]	=	sentropy(Pf_avg(if_min:if_max));
Kf_avg			=	    kurtosis(Pf_avg(if_min:if_max));

[EfdB_med, EmaxfdB_med]	=	sentropy(PfdB_med(if_min:if_max));
KfdB_med			=	kurtosis(PfdB_med(if_min:if_max));

[EfdB_avg, EmaxfdB_avg]	=	sentropy(PfdB_avg(if_min:if_max));
KfdB_avg			=	    kurtosis(PfdB_avg(if_min:if_max));

metrics.entropy.med=Ef_med/Emaxf_med;
metrics.entropy.avg=Ef_avg./Emaxf_avg;
metrics.kurtosis.med=Kf_med;
metrics.kurtosis.avg=Kf_avg;
metrics.kurtosis.avgdB=KfdB_avg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Peak Picking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if hf(1)	~=	-1
    Idebug.fig=[20 21];
    Idebug.title='PfdB_avg';
else
    Idebug=[];
end
[metrics.avg.peaks]=peak_picker_Thode(PfdB_avg,F,df_search,[f_min f_max],threshold,Idebug);
if hf(1)	~=	-1
    Idebug.fig=[22 23];
    Idebug.title='PfdB_med';
else
    Idebug=[];
end
[metrics.med.peaks]=peak_picker_Thode(PfdB_med,F,df_search,[f_min f_max],threshold,Idebug);



%%%%%%%%%%%%%%%%%%
%%	Plotting  %%
%%%%%%%%%%%%%%%%%%
if hf(1) ~= -1
    
    %%Raw spectrogram, units of power spectral density
    figure(hf(1));
    subplot(2,3,1);
    imagesc(T,F/1e3,PSD_dB);
    axis on; axis xy;
    title(['Spectrogram, mean total power = ' num2str(10*log10(P_tot)) ' dB']);
    xlabel('Time   (s)');
    ylabel('Frequency   (kHz)');
    try
    xlim(gca,xlim(ha_ref));
    ylim(gca,ylim(ha_ref));
    caxis(gca,caxis(ha_ref));
    end
    
    %%%%Ridge detector results
    % subplot(3,2,2);
    % imagesc(T,F/1e3,detect);
    % axis on; axis xy;
    % title('Ridge detection results');
    % xlabel('Time   (s)');
    % ylabel('Frequency   (kHz)');
    % xlim(gca,xlim(ha_ref));
    % ylim(gca,ylim(ha_ref));
    % caxis(gca,caxis(ha_ref));
    
    
    %%%Entropy over time.
    subplot(2,3,4);
    plot(T, PtdB);
    title_str=sprintf('Time avg: Entropy (linear) = %6.4f  Kurtosis(linear) = %6.2f,  Time avg: Entropy (dB) = %6.4f  Kurtosis(dB) = %6.2f' , ...
        (Et/Emaxt * 100),Kt, 100*EtdB/EmaxtdB,KtdB);
    title(title_str)
    xlabel('Time   (s)');
    ylabel('Power   (db)');
    
    %%mean spectrum, and entropy and kurtosis of time-average and median spectrum..
    subplot(2,3,2);
    plot(PfdB_avg, F/1e3,10*log10(Pf_avg),F/1e3);
    title_str=sprintf('Freq mean: E/Emax = %6.4f %%  K = %6.2f,  E/Em(dB) = %6.4f %% K(dB) = %6.2f', ...
        (Ef_avg/Emaxf_avg * 100),Kf_avg, 100*EfdB_avg/EmaxfdB_avg,KfdB_avg);
    title(title_str)
    legend('dB avg,','linear avg');grid on
    %title(['Entropy = ' num2str(E/Emax * 100) '%, Kurtosis = ' num2str(K)]);
    ylabel('Freq   (kHz)');
    xlabel('Power   (db)');
    
    subplot(2,3,3);
    plot(PfdB_med,F/1e3,10*log10(Pf_med),F/1e3);
    title_str=sprintf(' E/Emax(linear median)=%6.4f %% K(linear median)=%6.2f' , ...
        100*Ef_med/Emaxf_med,Kf_med);
    title(title_str)
    legend('dB median','linear median');grid on
    %title(['Entropy = ' num2str(E/Emax * 100) '%, Kurtosis = ' num2str(K)]);
    ylabel('Freq   (kHz)');
    xlabel('Power   (db)');
    
    %%Add peaks to graphs...
    figure(hf(1));
    subplot(2,3,2)
    hold on
    plot_peaks(metrics.avg.peaks{1});
    
    
    figure(hf(1));
    subplot(2,3,3)
    hold on
    plot_peaks(metrics.med.peaks{1});
    
end

if exist('want_hough','var')&~isempty(want_hough)
    
    %%	Normalized Image data methods
    I	=	mat2gray(PSD);
    Edge_types	=	{'sobel';'prewitt';'roberts';'log';'canny'};
    ie	=	1;
    figure(hf(2));
    edge_detector(I, T, F, ie);
    
    
    %	Selection list
    hspopup	=	uicontrol(hf(2), 'Style','popupmenu',...
        'String', Edge_types,...
        'BackgroundColor', 'white',...
        'Units','characters','Position',[0, 0, 16, 4],...
        'Callback', @hspopup_callback); %#ok<NASGU>
    
end
    function	hspopup_callback(source, eventdata) %#ok<INUSD>
        ie		=	get(source,'Value');
        figure(hf(2));
        edge_detector(I, T, F, ie);
    end


end


%%	Subfunctions
function	[E, Emax]	=	sentropy(S)
S	=	abs(S);
S	=	S./sum(S);
E	=	-sum(S .* log2(S));

Emax	=	log2(length(S));

if	~isreal(E)
    warning('Entropy value not real');
    keyboard;
end
end

%	Image edge detection
function	edge_detector(I, T, F, ie)

Edge_types	=	{'sobel';'prewitt';'roberts';'log';'canny'};

if	~exist('ie', 'var') || isempty(ie) || ie < 1 || ie > length(Edge_types)
    ie	=	1;
end

subplot(2,2,1);
imagesc(T,F,I); colormap('gray');
axis on; axis xy;
title('Original');


%	Default is the Sobel edge detector
%for ie	=	1:length(Edge_types)
eType	=	Edge_types{ie};

BW{ie} = edge(I, eType, [], 'horizontal');
subplot(2,2,2);
imagesc(T,F,BW{ie}); colormap('gray');
axis on; axis xy;
title(eType);
%end
xlabel('Time   (s)');
ylabel('Frequency   (kHz)');


%	Hough Transform (from example)
%for ie	=	1:length(Edge_types)
eType	=	Edge_types{ie};

disp(['Using edge results from ' eType ' detector...']);
[H,theta,rho] = hough(BW{ie});
subplot(2,2,3);
imagesc(theta, rho, imadjust(mat2gray(H)));
title('Hough Transform');
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis xy, hold on;
colormap(hot);

%	Find peaks
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = theta(P(:,2));
y = rho(P(:,1));
hold on;
plot(x,y,'s','color','black');
hold off;

%	Find lines
lines = houghlines(BW{ie},theta,rho,P,'FillGap',5,'MinLength',7);
subplot(2,2,4);
imagesc(T,F,I); %colormap('gray');
axis on; axis xy;
hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(T(xy(:,1)),F(xy(:,2)),'LineWidth',2,'Color','green');
    
    % Plot beginnings and ends of lines
    plot(T(xy(1,1)),F(xy(1,2)),'x','LineWidth',2,'Color','yellow');
    plot(T(xy(2,1)),F(xy(2,2)),'x','LineWidth',2,'Color','red');
    
    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end

% highlight the longest line segment
plot(T(xy_long(:,1)),F(xy_long(:,2)),'LineWidth',2,'Color','red');
hold off;
title([eType ' based']);
%	pause;
%end

end

function plot_peaks(peaks)
plot(10*log10(peaks.isi.P),peaks.isi.F/1000,'gs');  %green square is ISI_tone estimate...
for I=1:length(peaks.adp.PdB)
    Fplot=[ peaks.adp.F(I)+peaks.adp.bandwidth(I)*[-0.5 0.5]];
    Pplot=peaks.adp.PdB(I)*[1 1];
    plot(peaks.adp.PdB,peaks.adp.F/1000,'b^',Pplot,Fplot/1000,'b-');
    
end
for NN=1:length(peaks.adp.F)
   fprintf('Adp Frequency: %6.2f\n',peaks.adp.F(NN)); 
end

for NN=1:length(peaks.isi.F)
   fprintf('ISI Frequency: %6.2f\n',peaks.isi.F(NN)); 
end

end



