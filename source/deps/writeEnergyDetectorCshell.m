% Melania Guerra and Aaron Thode
% function scriptname=writeEnergyDetectoCshell(params)
%
% Function will write a Cshell script for calculating PSD or SEL
% It permits monitoring of the first channel of a multi-channel file only!
%
% Input: structure params:
%  fields:
%   program_dir: Location of JAVA program
%   file_dir: Directory of input files of interest...
%   exten:  file extension to search for in directory..or complete file name
%   dir_out: string with full or relative output directory path
%   chan:  which channel of multichannel file to use.  Hard-wired to '0'
%       for now.
%   dumpsize: Optional argument for number of samples to read into memory
%   nstart: Optional argument for start sample in file to look for
%       detection.  Set to 0 to examine entire file.  If not input, set to
%       '0'
%   nsamples: Optional argument for number of samples to examine.  Set to 0
%       to examine entire file.  If not input, set to '0'.
%   Nfft:  FFT size
%   Fs: sampling frequency
%   ovlap: number of samples to overlap FFT samples (not a percentage)
%   f_low: Lower frequency range in Hz to search
%   f_high: Upper frequency range in Hz to search
%   eq_time: Equalization time in seconds
%   bandwidth: bandwidth in Hz for a subdetector
%   threshold: dB threshold for a detection
%   bufferTime: time in seconds to bound detection when writing data to
%       .snips file.  If less than zero no snips file written
%   TolTime: How much time in seconds to wait before logging as new
%   detection
%   MinTime: Minimum time in seconds for individual detection to last
%   MinTimeTotal:  Minimum time in seconds for complete detection to last.
%       Useful for collecting a sequence of pulses into one detection (e.g.
%       gray whale calls).  If not input, set to MinTime;
%   MaxTime: Maximum time in seconds for detection to last
%   snips_chc:  if zero, don't write a snips file; if 1, write a
%       single-channel snips file, using channel of detection; if two write out all channels of snips,
%       regardless of which channel was used for detection
%   debug: if not zero, will write out SEL or SNR output
% Examples:
% params.file_dir='/Users/thode/Projects/Sitka_Sperm/Sperm_Sitka_2007/Bprobe_download/July17.dir/031_good.dir/sound';
% %params.exten='*_07171109*mt';
% params.exten='*mt';
% params.dir_out='DetectionFiles.dir';
% params.dumpsize=5*60*4096-1;
% params.Fs=4096;
% %params.nstart=60*params.Fs;
% %params.nsamples=15*params.Fs;
% 
% params.nstart=0;
% params.nsamples=0;
% params.Nfft=128;
% params.ovlap=round(.75*params.Nfft);
% params.f_low=800;
% params.f_high=1800;
% params.eq_time=10;
% params.bandwidth=1000;
% params.threshold=10;
% params.bufferTime=-1;
% params.TolTime=0.01;
% params.MinTime=0;
% params.MaxTime=0.5;
% params.snips_chc=2;  %For multiple channel outputs
% params.debug=0;
%
% % Oct 2007
%  Update Dec 2008

function script_name=writeEnergyDetectorCshell(params)
% file_dir, dir_out, main_dir, exten, Nfft, ovlap, Fs, dumpsize,
%nstart, nsamples, sec_avg, f_high, f_low, eq_time, bandwidth, threshold, bufferTime, TolTime, MinTime, debug)
script_name=sprintf('masterEDET_%i.scr',floor(rand*1000));
%script_name=['masterEDET'.scr' ;
f_det = fopen(script_name,'wt');
fprintf(f_det, '#!/bin/csh \n\n\n');
fprintf(f_det,'set dirout = "%s"    \n',params.dir_out);
fprintf(f_det,'set Npsd = 0        \n');
fprintf(f_det,'set Nsel = 0        \n');
fprintf(f_det,'set Ndet = %i        \n',1);

%At present, set input channel to zero (i.e., no processing of channel 2 of
%  a multichannel file currently allowed.
fprintf(f_det,'set chan = %i        \n',0);

if ~isfield(params,'dumpsize')
    fprintf(f_det,'set dumpsize = 1000000\n');
else
    fprintf(f_det,'set dumpsize = %11.0f    \n',params.dumpsize);
end

if ~isfield(params,'nstart')
    fprintf(f_det,'set nstart = 0\n');
else
    fprintf(f_det,'set nstart = %11.0f     \n',params.nstart);
end


if ~isfield(params,'nsamples')&&~isfield(params,'nsample')
    fprintf(f_det,'set nsamples = 0\n');
else
    fprintf(f_det,'set nsamples = %11.0f       \n',params.nsamples);
end


fprintf(f_det,'set Nfft = %i        \n',params.Nfft);
fprintf(f_det,'set Fs = %i    \n',params.Fs);
fprintf(f_det,'set ovlap = %i        \n',params.ovlap);
fprintf(f_det,'set sec_avg = %i        \n',0);
fprintf(f_det,'set flo_det = %i        \n',params.f_low);
fprintf(f_det,'set fhi_det = %i        \n',params.f_high);
fprintf(f_det,'set eq_time = %i        \n',params.eq_time);
fprintf(f_det,'set bandwidth = %i        \n',params.bandwidth);
fprintf(f_det,'set threshold = %i        \n',params.threshold);
fprintf(f_det,'set bufferTime = %i        \n',params.bufferTime);
fprintf(f_det,'set TolTime = %i        \n',params.TolTime);
fprintf(f_det,'set MinTime = %i        \n',params.MinTime);
if isfield(params,'MinTimeTotal')
    fprintf(f_det,'set MinTimeTotal = %i        \n',params.MinTimeTotal);
else
    fprintf(f_det,'set MinTimeTotal = %i        \n',params.MinTime);
end

fprintf(f_det,'set MaxTime = %i        \n',params.MaxTime);
if isfield(params,'snips_chc')
    fprintf(f_det,'set snips_chc = %i        \n',params.snips_chc);
else
    fprintf(f_det,'set snips_chc = %i        \n',1);
end
fprintf(f_det,'set debug = %i        \n\n',params.debug);

 if ~isfield(params,'program_dir')
     params.program_dir=load_program_location;
 end

command_str = sprintf('java -jar %s $i $dirout $Npsd $Nsel $Ndet $chan $nstart $nsamples $dumpsize $Nfft $Fs $ovlap $flo_det $fhi_det $eq_time $bandwidth $threshold $bufferTime $TolTime $MinTime $MinTimeTotal $MaxTime $snips_chc $debug',params.program_dir);

fprintf(f_det,'foreach i (%s/*%s)       \n',params.file_dir, params.exten);
fprintf(f_det,'echo " Starting processing of " $i\n');
command_str = [command_str ' >> out.txt\n'];
fprintf(f_det,command_str);
fprintf(f_det,'echo " Finished processing of " $i\n');
fprintf(f_det,'end\n');
fclose(f_det);
eval(sprintf('!chmod +x %s',script_name));
%!chmod +x masterEDET.scr


