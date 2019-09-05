%function [output,tabs_start,tabs_end,azigram_param]=get_DASAR_ambient_statistics(myfile,tstart,tsample,Nfft,ovlap,azigram_param)
% Inputs: 
%         myfile:  full pathname of GSI file to load
%         tstart:  time in seconds to begin loading data
%         tsample: seconds of data to load (Inf) loads all data
%         Nfft:  FFT size of data
%         ovlap:  fractional overlap between Nfft samples
%         azigram_param:  structure with information specific to azigram
%            including...
%               time_chunk: 600  %%Seconds of data to process as azigram at
%                       a time; saves RAM memory and helps with clock drift
%             percentiles: [0 0.2500 0.5000 0.7500 1]  %percentiles for
%                           noise stats
%                  frange: [50 475]  %Frequency range over which to gather
%                  stats
%                 sec_avg: '0'  %%%Ensures no averaging in azigram
%            f_transition: 300  %Frequency at which active and reactive
%                   intensity are swapped to address instrumental error
%                azi_grid: [1×36 double]
%               ItoE_grid: [1×21 double]
%             KEtoPE_grid: [-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6]
%                PdB_grid: [1×46 double]
%     IntensityPhase_grid: [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90]
%
%  Outputs:
%       tabs_start,tabs_end: datenumbers of start and end time of data
%        output:
%               PdB_percentiles: [6×109×5 double]
%                PdB_hist: [4×6×109×45 double]
%                azi_hist: [4×6×109×35 double]
%                ItoE_hist: [4×6×109×20 double]
%                KEtoPE_hist: [4×6×109×12 double]
%                IntensityPhase_hist: [4×6×109×18 double]
%                tabs: [7.3744e+05 7.3744e+05 7.3744e+05 7.3744e+05 7.3744e+05 7.3744e+05]
%                params: [1×1 struct]  %parrots the input parameters, plus
%                       contains myfile
%                FF: [109×1 double]
function [output,tabs_start,tabs_end,azigram_param]=get_DASAR_ambient_statistics(myfile,tstart,tsample,Nfft,ovlap,azigram_param)

[x,~,head]=readGSIfile(myfile,tstart,tsample,1:3,'seconds','nocalibrate');
%x=x(1+round(tstart*head.Fs):end,:);    
x=calibrate_GSI_signal(x, 'DASARC')';
azigram_param.brefa=head.brefa;
azigram_param.goodName=myfile;
tabs_start=head.tabs_start+datenum(0,0,0,0,0,tstart);
tabs_end=tabs_start+datenum(0,0,0,0,0,tsample);

fprintf('File: %s Time start: %s Time end: %s\n', myfile,datestr(tabs_start,31),datestr(tabs_end,31));
%cd(current_dir)

%%%%Parameters for the automated processing
azigram_param.tabs_start=tabs_start;
azigram_param.debug.image=false;  %%Image azigram and spectrogram for each round

%%%%azigram statistics on time series
output=azigram_noise_statistics_compute(x,head.Fs,Nfft,ovlap,azigram_param);
fprintf('Finished processing...\n');

end