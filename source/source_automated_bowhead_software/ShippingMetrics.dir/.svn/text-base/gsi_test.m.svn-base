%%	GSI test script
clc;
clear all;




%%	Setup file and processing options
file_name	=	'/Volumes/Data/Shell2012_GSI_Data/S312gsif/S312F0/S312F0T20120929T000000.gsi';

t_win	=	60;
t_int	=	10*60;
bwidth	=	5;	% in # of fft bins
f_range	=	[50 450];
Nfft	=	256;
overlap =	0.50;




%%	Calculate statistics
SD_metrics	=	gsi_shipdet_metrics(file_name, t_win, t_int, bwidth, f_range, Nfft, overlap);

Num_peaks	=	zeros(length(SD_metrics.Max_Peaks),1);
for	ii = 1:length(SD_metrics.Max_Peaks)
	Num_peaks(ii)	=	length(SD_metrics.Max_Peaks{ii}.Freq);
end

SD_metrics.T	=	SD_metrics.T/60/60;

%%	Plot values
figure;

subplot(2,4,1);
plot(SD_metrics.T, 10*log10(abs(SD_metrics.P_tot)));
xlabel('Time   (hr)');	xlim([0 24]);	grid on
ylabel('Power  (dB)');
title('Total power');

subplot(2,4,2);
plot(SD_metrics.T, Num_peaks, 'k-');%, 'BarWidth',1.0,'LineStyle','none');
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('# of peaks detected');
title('Signficant peaks');

subplot(2,4,3);
plot(SD_metrics.T, SD_metrics.Time_entropy);
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('Entropy  (%)');
title('Time entropy, as % of max');

subplot(2,4,4);
plot(SD_metrics.T, SD_metrics.Time_kurtosis);
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('Kurtosis');
title('Time kurtosis, > 3 implies signifcant peakiness');

subplot(2,4,5);
plot(SD_metrics.T, SD_metrics.Freq_entropy);
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('Entropy  (%)');
title('Freq entropy, as % of max');

subplot(2,4,6);
plot(SD_metrics.T, SD_metrics.Freq_kurtosis);
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('Kurtosis');
title('Freq kurtosis, > 3 implies signifcant peakiness');

subplot(2,4,7);
plot(SD_metrics.T, SD_metrics.avg_CorrCoef);
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('Correllation Coefficients');
title('Mean Correllation Coefficients');

subplot(2,4,8);
plot(SD_metrics.T, SD_metrics.avg_Conf);
xlabel('Time   (hr)');	xlim([0 24]);	grid on 
ylabel('Confidence');
title('Mean CorrCoeff confidence level, < 0.05 preferred');

