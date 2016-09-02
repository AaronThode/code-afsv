%function [param,LB,UB]=TOC_params(param_chc,param)

function [param,LB,UB]=TOC_params(param_chc,param)
%%TOC params for test data set optimization
%param_chc='Feb14_2008_0830Site3a';
LB=[];UB=[];
% if strcmpi((computer),'maci'),
%     optdir='/Users/Shared/Projects/Arctic_2007/Optimize';
% elseif strcmpi((computer),'mac')
%     optdir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/Arctic_2007/Optimize';
% end

%Choose neural network
[param.net]=load_neural_network_path(param_chc);
        
load_base_scenario;
    function load_base_scenario
        disp('Shell09_initial: base parameter set loaded');
        param.Fs=1000;
        param.Nfft=256;
        param.ovlap=7/8;
        
        %param.energy.program_dir='~/Documents/workspace/PreProcess_PSD_detection/PreProcess.jar';
        param.energy.file_dir='.';
        param.energy.exten='gsi';
        param.energy.dir_out='.';
        param.energy.Nfft=256;
        param.energy.ovlap=128;
        param.energy.f_low = 25;
        param.energy.f_high = 500;
        param.energy.eq_time = 23.8;
        param.energy.bandwidth = 37;
        param.energy.threshold = 6 ;
        param.energy.bufferTime = 1 ;
        param.energy.TolTime = 0.05;
        param.energy.MinTime = 0.1;
        param.energy.MaxTime = 6;
        param.energy.debug=0;
        
        %
        param.interval_remove.Ndet=10;
        param.interval_remove.Nmiss=3;
        param.interval_remove.ICItol=0.75;
        param.interval_remove.ICI_removal_alg='TimingOnly';
        % param.interval_remove.names={'peak_freq','peak_duration'};
        %param.interval_remove.tol_feature=[ 5+(0.5*param.energy.bandwidth) 0.75*1000];
        
        %%Main run
        param.interval_remove.names={'min_freq','max_freq','total_duration'};
        param.interval_remove.tol_feature=[[1  1.5]*(2*param.energy.bandwidth) 5*0.75];
        
        %param.interval_remove.names={'min_freq','total_duration'};
        %param.interval_remove.tol_feature=[[1  ]*(2*param.energy.bandwidth) 0.75];
        
        %param.interval_remove.names={'min_freq','total_duration'};
        %param.interval_remove.tol_feature=[[1  ]*(2*param.energy.bandwidth) 5*0.75];
        
        param.interval_remove.debug_time=datenum(2009,9,12,1,0,0);
        param.interval_remove.debug_index=20;
        
        % Parameters for median filtering..
        param.median.on=0;
        param.median.size=[0.2 20];
        
        %parameters for gaussian filtering..
        param.filter.on=0;  %If one, asymetric gaussian, if two, symmetric gaussian
        param.filter.size=[0.2 0.2];  %Filter in units of [sec Hz] for spectrgram
        param.filter.sigma=0.5;  %Units of pixel size.
        
        
        %Morphological analysis.  Other parameters are derived from above...
        param.morph.eq_time=param.energy.bufferTime-0.25;  %Time to use to create background ambient noise estimate.
        param.morph.threshold_fudge=0.75;
        param.morph.SNRmin=10;  %min dB level, for thresholding background noise
        param.morph.dynamic_range=3.11;  %dB below maximum a pixel is allowed to have
        param.morph.dynamic_range2=7;  %dB below maximum a pixel is allowed to have, horizontal (regional) maximum
        
        param.morph.gap_f=11;  %For dilation results
        param.morph.gap_t=0.04;  %MIGHT NEED TO CHANGE BACK TO 0.3 For dilation results, might need
        
        param.merge.ovlap=0.25;
        param.merge.max_frequency_separation=50;
        param.merge.gap_t=0.1;
        param.merge.gap_f=20;
        
        param.morph.background.on=1;  %Execute contour linking and processing
        param.morph.background.gap_f=param.merge.gap_f;
        param.morph.background.gap_t=param.merge.gap_t;
        
        param.morph.MinFreq=0;  %Minimum frequency to crop image
        param.morph.MaxFreq=500; %Maximum frequency to crop image.
        
        param.morph.duration.min=0.15;
        param.morph.duration.max=param.energy.MaxTime;
        
        %param.morph.time_band_product.min=9.529000000000001e-01;
        param.morph.time_band_product.min=1;
        param.morph.time_band_product.max=Inf;
        
        param.morph.robust_fmin.min=0;
        param.morph.robust_fmin.max=500;
        
        param.morph.robust_fmax.min=0;
        param.morph.robust_fmax.max=500;
        
        param.morph.robust_bandwidth.max=500;
        param.morph.robust_bandwidth.min=0;
        
        param.morph.local_bandwidth.max=500;
        param.morph.local_bandwidth.min=0;
        
        param.morph.Orientation.max=90;
        param.morph.Orientation.min=-90;
        
        %param.morph.Centroid_freq.max=200;
        %param.morph.Centroid_freq.min=85;
        
        param.morph.Centroid_freq.max=500;
        param.morph.Centroid_freq.min=0;
        
        param.morph.Eccentricity.max=1;
        param.morph.Eccentricity.min=0.0;
        
        param.morph.percent_safety_trim=0.75;
        param.morph.want_contour=[];
        
        
        
        %         param.feature.names={'fmin','fmax','BoundingBox'};
        %         param.feature.index=[1 1 3];
        %         param.feature.weight=[1 1 0.25];
        %         param.feature.tol=[.2 .2 1.5 ];
        %         param.feature.Nshapes=3;
        %         param.feature.compare_function_str='relative_error';
        
        param.feature.names={'robust_fmin','robust_fmax','duration','Orientation','Eccentricity', ...
            'local_bandwidth','robust_bandwidth','median_local_kurtosis','Solidity','time_band_product','Nharmonic','harmonic_spacing','SNR','SEL'};
        param.feature.global_names={'Totalfmin','Totalfmax','Totalduration','Contour_Area','Contour_global_bandwidth','Contour_duration', ...
            'Contour_local_bandwidth','Contour_fmin','Contour_fmax'};  %Global names are values that exist once per call
        
        param.feature.Nsegments=1;
        param.feature.index=ones(1,length(param.feature.names));
        %param.feature.weight=[1 1 0.25 0 0 0 0 0 0 0 0];
        %param.feature.tol=[30 100 2  720 2 1000 1000 1000 Inf Inf Inf];
        
        param.feature.compare_function_str='correlate_image'; %'hausdorff_image';
        % param.feature.compare_function_str='absolute_error';
        
        %param.feature.matched_filter_tol=0.49;
        param.feature.hausdorff.tol=2;
        param.feature.xcorr.tol=1-sqrt(1/3);
        param.localize.min_weight=0.0;
        
        
        param.net.skip_feature={'duration_debug','ctime','ctime_min','ctime_debug','miss'};
        
        param.final_filter.names={'Contour_global_bandwidth','Contour_fmax','Contour_fmin'};
        param.final_filter.criteria=[0 20 5; 300 400 500];
    end

param.label=param_chc;

%%Alter baseline parameters...

if ~isempty(findstr(param_chc,'NeuralNetupdated'))  %%Strip away further keywords e.g., Threshholdadjust. etc.
   param_chc= 'NeuralNetupdated';
end


switch param_chc
    case {'Shell09_initial','Shell07_initial'}
        
    case {'Shell09_DoubleNeuralNetFirstRun'}
         
    case {'NeuralNetupdated'}
        %%Used to generate final results...
        param.interval_remove.on=1;
        param.interval_remove.ICI_removal_alg='BearingAndTiming';  %  If 'BearingAndTiming', use new bearing and timing algorithm for ICI removal. Otherwise timing only
        param.interval_remove.Ndet=20;
        param.interval_remove.Nmiss=12;
        param.interval_remove.ICItol=1.00;  %How close does a candidate have to be to predicted ICI time? (sec)
        
       %param.interval_remove.names={'min_freq',  'bearing'};
        %param.interval_remove.tol_feature=[[1 ]*(2*param.energy.bandwidth)    17];
        
        param.interval_remove.names={'bearing'};  %What features should be used when restricting ICI candidates?
        param.interval_remove.tol_feature=15;  %Note that 10 is cleaner, but the weak detections need higher tolerances
       
        %param.interval_remove.ICI_std= 1.0;  %How close do adjacent ICI
        %       estimates have to be (sec)-Dawn investigated.
        param.interval_remove.ICI_std=0.4;
        param.interval_remove.Nstd=7;      % How many ICIs must pass the 'ICI_std' test?
        
       % param.interval_remove.debug_time=datenum(2008,8,26,11,39,10);%26-Aug-2008 10:11:37; 2008,8,26,3,58,45
       % param.interval_remove.debug_index=20;         
       
    case {'Shell09_new_interval_test'}
        param.interval_remove.Ndet=20;
        param.interval_remove.Nmiss=12;
        param.interval_remove.ICItol=0.75;
        
        param.interval_remove.names={'bearing'};
        param.interval_remove.tol_feature=10;  %Note that 10 is cleaner, but the weak
        
        
        
    case {'Shell08_final_parameters','Shell08_March_optimized'}
        I=1;
        disp('Shell08_March_optimized');
        
        
        %%Main run
        param.interval_remove.names={'min_freq','max_freq','total_duration'};
        param.interval_remove.tol_feature=[[1  1]*(2*param.energy.bandwidth) 0.75];
        
        %param.interval_remove.names={'min_freq','total_duration'};
        %param.interval_remove.tol_feature=[[1  ]*(2*param.energy.bandwidth) 0.75];
        
        param.interval_remove.debug_time=datenum(2008,9,15,8+1,3,1);
        param.interval_remove.debug_time=datenum(2008,8,23,5,1,10);
        %param.interval_remove.debug_time=datenum(2008,8,24,13,34,0);
        %param.interval_remove.debug_time=datenum(2008,8,31,4,13,33);
        %param.interval_remove.debug_time=datenum(2008,8,31,4,30,18);
        
        param.interval_remove.debug_index=5;
        
    case 'Shell08_DawnGrebner'
        param.energy.bandwidth=100;
        param.interval_remove.ICI_removal_alg='BearingAndTiming';  %  If 'BearingAndTiming', use new bearing and timing algorithm for ICI removal. Otherwise timing only
        param.interval_remove.Ndet=20;
        param.interval_remove.Nmiss=12;
        param.interval_remove.ICItol=1.00;  %How close does a candidate have to be to predicted ICI time? (sec)
        
        param.interval_remove.names={'bearing'};  %What features should be used when restricting ICI candidates?
        param.interval_remove.tol_feature=[ 17];  %Note that 10 is cleaner, but the weak
        
        param.interval_remove.ICI_std=0.4;  %How close do adjacent ICI estimates have to be (sec)
        param.interval_remove.Nstd=7;      % How many ICIs must pass the 'ICI_std' test?
        
    case 'Shell08_initial'
        I=1;
        disp('Shell08_initial');
        param.Fs=1000;
        param.Nfft=128;
        param.ovlap=7/8;
        
        %param.energy.program_dir='~/Documents/workspace/PreProcess_PSD_detection/PreProcess.jar';
        param.energy.file_dir='.';
        param.energy.exten='gsi';
        param.energy.dir_out='.';
        param.energy.Nfft=128;
        param.energy.ovlap=64;
        param.energy.f_low = 50;
        param.energy.f_high = 500;
        param.energy.eq_time = 23.8;
        param.energy.bandwidth = 37;
        param.energy.threshold = 6 ;
        param.energy.bufferTime = 1 ;
        param.energy.TolTime = 0.05;
        param.energy.MinTime = 0.1;
        param.energy.MaxTime = 5;
        param.energy.debug=0;
        
        
        %         param.interval_remove.Ndet=7;
        %         param.interval_remove.Nmiss=1;
        %         param.interval_remove.ICItol=.3;
        %
        param.interval_remove.Ndet=8;
        param.interval_remove.Nmiss=2;
        param.interval_remove.ICItol=1;
        param.interval_remove.names={'peak_freq','peak_duration'};
        param.interval_remove.tol_feature=[ 5+(0.5*param.energy.bandwidth) 0.75*1000];
        param.interval_remove.debug_index=5;
        
        % Parameters for median filtering..
        param.median.on=0;
        param.median.size=[0.2 20];
        
        %parameters for gaussian filtering..
        param.filter.on=0;  %If one, asymettic gaussian, if two, symmetric gaussian
        param.filter.size=[0.2 0.2];  %Filter in units of [sec Hz] for spectrgram
        param.filter.sigma=0.5;  %Units of pixel size.
        
        
        %Morphological analysis.  Other parameters are derived from above...
        param.morph.eq_time=param.energy.bufferTime-0.25;  %Time to use to create background ambient noise estimate.
        param.morph.threshold_fudge=0.75;
        param.morph.SNRmin=10;  %min dB level, for thresholding background noise
        param.morph.dynamic_range=3;  %dB below maximum a pixel is allowed to have
        
        param.morph.gap_f=10;  %For dilation results
        %param.morph.gap_t=0.15;  %MIGHT NEED TO CHANGE BACK TO 0.3 For dilation results, might need
        param.morph.gap_t=0.05;  %MIGHT NEED TO CHANGE BACK TO 0.3 For dilation results, might need
        
        param.morph.MinFreq=0;  %Minimum frequency to crop image
        param.morph.MaxFreq=500; %Maximum frequency to crop image.
        
        param.morph.duration.min=0.15;
        param.morph.duration.max=param.energy.MaxTime;
        
        %param.morph.time_band_product.min=9.529000000000001e-01;
        param.morph.time_band_product.min=0.25;
        param.morph.time_band_product.max=Inf;
        
        param.morph.robust_fmin.min=20;
        param.morph.robust_fmin.max=400;
        
        param.morph.robust_fmax.min=50;
        param.morph.robust_fmax.max=450;
        
        param.morph.local_bandwidth.max=150;
        param.morph.local_bandwidth.min=0;
        
        param.morph.Orientation.max=60;
        param.morph.Orientation.min=-60;
        
        %param.morph.Centroid_freq.max=200;
        %param.morph.Centroid_freq.min=85;
        
        param.morph.Centroid_freq.max=500;
        param.morph.Centroid_freq.min=0;
        
        param.morph.Eccentricity.max=1;
        param.morph.Eccentricity.min=0.7;
        
        param.morph.percent_safety_trim=0.75;
        
        param.merge.ovlap=0.25;
        param.merge.max_frequency_separation=50;
        param.merge.gap_t=1.5*param.morph.gap_t;
        param.merge.gap_f=param.morph.gap_f;
        
        
        %         param.feature.names={'fmin','fmax','BoundingBox'};
        %         param.feature.index=[1 1 3];
        %         param.feature.weight=[1 1 0.25];
        %         param.feature.tol=[.2 .2 1.5 ];
        %         param.feature.Nshapes=3;
        %         param.feature.compare_function_str='relative_error';
        
        param.feature.names={'robust_fmin','robust_fmax','BoundingBox','Orientation','Eccentricity', ...
            'local_bandwidth','robust_bandwidth','Solidity','Area','SNR','SEL'};
        param.feature.Nshapes=3;
        param.feature.index=[1 1 3 1 1 1 1 1 1 1 1];
        param.feature.weight=[1 1 0.25 0 0 0 0 0 0 0 0];
        param.feature.tol=[30 100 2  720 2 1000 1000 1000 Inf Inf Inf];
        
        param.feature.compare_function_str='correlate_image'; %'hausdorff_image';
        % param.feature.compare_function_str='absolute_error';
        
        param.feature.matched_filter_tol=0.49;
        param.feature.hausdorff.tol=2;
        param.feature.xcorr.tol=1-sqrt(1/3);
        param.localize.min_weight=0.7;
        
        
        %param.net.name='Nnetwork_Shell08_initial_10hidden_pca_20090119T232220.mat';
        %param.net.name='Nnetwork_Shell08_initial_10hidden_20080821.mat';
        param.net.name='Nnetwork_Shell08_initial_10cells_20080821_20080929.mat';
    case 'Shell08_Final_February_archive'
        I=1;
        disp('Shell08_initial');
        param.Fs=1000;
        param.Nfft=128;
        param.ovlap=7/8;
        
        %param.energy.program_dir='~/Documents/workspace/PreProcess_PSD_detection/PreProcess.jar';
        param.energy.file_dir='.';
        param.energy.exten='gsi';
        param.energy.dir_out='.';
        param.energy.Nfft=256;
        param.energy.ovlap=128;
        param.energy.f_low = 50;
        param.energy.f_high = 500;
        param.energy.eq_time = 23.8;
        param.energy.bandwidth = 37;
        param.energy.threshold = 6 ;
        param.energy.bufferTime = 1 ;
        param.energy.TolTime = 0.05;
        param.energy.MinTime = 0.1;
        param.energy.MaxTime = 5;
        param.energy.debug=0;
        
        
        %         param.interval_remove.Ndet=7;
        %         param.interval_remove.Nmiss=1;
        %         param.interval_remove.ICItol=.3;
        %
        param.interval_remove.Ndet=8;
        param.interval_remove.Nmiss=2;
        param.interval_remove.ICItol=1;
        param.interval_remove.names={'peak_freq','peak_duration'};
        param.interval_remove.tol_feature=[ 5+(0.5*param.energy.bandwidth) 0.75*1000];
        param.interval_remove.debug_index=5;
        
        % Parameters for median filtering..
        param.median.on=0;
        param.median.size=[0.2 20];
        
        %parameters for gaussian filtering..
        param.filter.on=0;
        param.filter.size=[0.2 50];  %Filter in units of [sec Hz] for spectrgram
        param.filter.sigma=0.5;  %Units of pixel size.
        
        param.pulse_percent=0.8750;
        param.pulse_Tclear=0.58;
        param.pulse_median_clear=1;  %If exists, then equalize times after pulse end.  Added Feb. 3, 2008
        
        %Morphological analysis.  Other parameters are derived from above...
        param.morph.eq_time=param.energy.bufferTime;  %Time to use to create background ambient noise estimate.
        param.morph.threshold_fudge=0.75;
        param.morph.SNRmin=param.energy.threshold-1;  %min dB level, for thresholding background noise
        param.morph.dynamic_range=10;  %dB below maximum a pixel is allowed to have
        
        param.morph.gap_f=10;  %For dilation results
        param.morph.gap_t=0.15;  %MIGHT NEED TO CHANGE BACK TO 0.3 For dilation results, might need
        %param.morph.gap_t=0.05;  %MIGHT NEED TO CHANGE BACK TO 0.3 For dilation results, might need
        
        param.morph.MinFreq=0;  %Minimum frequency to crop image
        param.morph.MaxFreq=500; %Maximum frequency to crop image.
        
        param.morph.duration.min=0.2;
        param.morph.duration.max=param.energy.MaxTime;
        
        %param.morph.time_band_product.min=9.529000000000001e-01;
        param.morph.time_band_product.min=1;
        param.morph.time_band_product.max=1000;
        
        param.morph.robust_fmin.min=20;
        param.morph.robust_fmin.max=400;
        
        param.morph.robust_fmax.min=50;
        param.morph.robust_fmax.max=450;
        
        param.morph.local_bandwidth.max=150;
        param.morph.local_bandwidth.min=0;
        
        param.morph.Orientation.max=60;
        param.morph.Orientation.min=-60;
        
        %param.morph.Centroid_freq.max=200;
        %param.morph.Centroid_freq.min=85;
        
        param.morph.Centroid_freq.max=500;
        param.morph.Centroid_freq.min=0;
        
        param.morph.Eccentricity.max=1;
        param.morph.Eccentricity.min=0.1;
        
        param.morph.percent_safety_trim=0.75;
        
        param.merge.ovlap=0.25;
        param.merge.max_frequency_separation=50;
        param.merge.gap_t=1.5*param.morph.gap_t;
        param.merge.gap_f=param.morph.gap_f;
        
        
        %         param.feature.names={'fmin','fmax','BoundingBox'};
        %         param.feature.index=[1 1 3];
        %         param.feature.weight=[1 1 0.25];
        %         param.feature.tol=[.2 .2 1.5 ];
        %         param.feature.Nshapes=3;
        %         param.feature.compare_function_str='relative_error';
        
        param.feature.names={'robust_fmin','robust_fmax','BoundingBox','Orientation','Eccentricity', ...
            'local_bandwidth','robust_bandwidth','Solidity','Area','SNR','SEL'};
        param.feature.Nshapes=3;
        param.feature.index=[1 1 3 1 1 1 1 1 1 1 1];
        param.feature.weight=[1 1 0.25 0 0 0 0 0 0 0 0];
        param.feature.tol=[30 100 2  720 2 1000 1000 1000 Inf Inf Inf];
        
        param.feature.compare_function_str='correlate_image'; %'hausdorff_image';
        % param.feature.compare_function_str='absolute_error';
        
        param.feature.matched_filter_tol=0.49;
        param.feature.hausdorff.tol=2;
        param.feature.xcorr.tol=1-sqrt(1/3);
        param.localize.min_weight=0;
        
        
        %param.net.name='Nnetwork_Shell08_initial_10hidden_pca_20090119T232220.mat';
        %param.net.name='Nnetwork_Shell08_initial_10hidden_20080821.mat';
        param.net.name='Nnetwork_Shell08_initial_10cells_20080821_20080929.mat';
        
        
    case 'NorthStar08_nofilt_airgun',
        param.Fs=1000;
        
        %param.energy.program_dir='~/Documents/workspace/PreProcess_PSD_detection/PreProcess.jar';
        param.energy.file_dir='.';
        param.energy.exten='sio';
        param.energy.dir_out='.';
        param.energy.Nfft=256;
        param.energy.ovlap=128+128/4;
        param.energy.f_low = 50;
        param.energy.f_high = 400;
        param.energy.eq_time = 60;
        param.energy.bandwidth = 30;
        % param.energy.threshold =8;
        param.energy.threshold=8;  %For Liberty
        param.energy.bufferTime = 1 ;
        param.energy.TolTime = 0.5;
        param.energy.MinTime = 0.02;
        param.energy.MaxTime = 10;
        param.energy.debug=0;
        
        %param.energy.nsamples=1*60*60*1000;
        
        %%The liberty files accidentally contained two days worth of
        %%data...
        %param.energy.nsamples=24*3600*param.Fs;
        
        param.interval_remove.Ndet=8;
        param.interval_remove.Nmiss=2;
        param.interval_remove.ICItol=1.5;
        param.interval_remove.names={'peak_freq','peak_duration'};
        param.interval_remove.tol_feature=[ 1.5*param.energy.bandwidth 0.75*1000];
        param.interval_remove.debug_index=200;
        
        
    case {'Liberty08_nofilt_airgun'}
        param.Fs=1000;
        
        %param.energy.program_dir='~/Documents/workspace/PreProcess_PSD_detection/PreProcess.jar';
        param.energy.file_dir='.';
        param.energy.exten='sio';
        param.energy.dir_out='.';
        param.energy.Nfft=256;
        param.energy.ovlap=128+128/4;
        param.energy.f_low = 30;
        param.energy.f_high = 400;
        param.energy.eq_time = 60;
        param.energy.bandwidth = 400;
        % param.energy.threshold =8;
        param.energy.threshold=10;  %For Liberty
        param.energy.bufferTime = 1 ;
        param.energy.TolTime = 0.5;
        param.energy.MinTime = 0.05;
        param.energy.MaxTime = 4.5;
        param.energy.debug=0;
        
        %%The liberty files accidentally contained two days worth of
        %%data...
        param.energy.nsamples=24*3600*param.Fs;
        
        param.interval_remove.Ndet=8;
        param.interval_remove.Nmiss=1;
        param.interval_remove.ICItol=0.1;
        
        
end


end
