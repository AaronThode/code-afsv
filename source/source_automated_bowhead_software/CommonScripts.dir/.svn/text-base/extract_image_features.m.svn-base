% %%%%%%%%%%%%%extract_image_features.m%%%%;
%function [features,final_image,Beq,labeled_raw]=extract_image_features(y,cstart,param,Idebug)
%% Take raw input data  y and start time ctime (in c-time), create an
%% equalized spectrogram, then apply morphological techniques to segment and extract features
%%%Inputs:
%     y: row vector of acoustic data
%     cstart: ctime of start of data
%     Idebug: integer specifying detail of debug output.. if 0, no output;
%           if one, only text debug output.  If 2, demo graphic outuput.
%           If 3, every stage has graphic output
%
%     param:  parameters for processing.  Must contain the following
%       fields
%%    param.Fs=1000;  sampling frequency in Hz
%     param.Nfft=128;  FFT size
%     param.ovlap=7/8; overlap fraction for computing spectrogram
%     param.morph: structure with following subfields...
%           Orientation: [1x1 struct]
%           SNRmin: 10
%           background.gap_f:
%           background.gap_t:
%           Centroid_freq: [1x1 struct]
%           duration: [1x1 struct]
%           dynamic_range: 3.1100
%           dynamic_range2: 7
%           Eccentricity: [1x1 struct]
%           eq_time: 0.7500  %Amount of buffer time before start of
%                   detection in y.  This portion of spectrogram will be removed.
%           equalization: [25x2 double]
%                       %if exists, this is an estimate of the long-term background
%                       power spectral density.  [Nfreq x 2] matrix, first column frequencies in Hz
%                       second column power spectral density in dB re 1uPa^2/Hz
%
%           gap_f: 11
%           gap_t: 0.0400
%           local_bandwidth: [1x1 struct]
%           MaxFreq: 500  %Maximum frequency to process image
%           MinFreq: 0    %Minimum frequency to process image
%           percent_safety_trim: 0.7500
%           robust_fmax: [1x1 struct]
%           robust_fmin: [1x1 struct]
%           threshold_chc: 'local_peaks'
%           time_band_product: [1x1 struct]
%           want_contour:  if 1, generate Tcall and Fcall, the image segments...
%   param.morph.background:
%           gap_f, gap_t: morphological elements for contour segments
%           on:  %if 1, activate contour linking.  Otherwise do not
%   param.merge:  parameters for merging segments (ridges) into single
%               detection
%
%           ovlap: 0.2500
%           max_frequency_separation: 50  % separation (Hz) for harmonics
%           gap_t: 0.1000
%           gap_f: 20
%
%
%     param.median.size=[0.2 50] (sec, Hz); %Used to derive bandwidth of
%       median filter
%     param.filter.size=[0.2 50] (sec,Hz); %LoG filter size
%     param.filter.sigma=0.75; %standard deviation for LoG parameter
%
%
%       Otsu option:
%       param.morph.threshold_fudge=what fraction of automated adaptive threshold should be
%           used for converting to binary image...
%           of the background noise
%


%
% % %%%Output %%%%%%%%
%
% Definitions:
%       Features are divided between ridges and contours.
%       A ridge is a local maximum- a "ridge" along topography
%       A contour is a region that exceeds a certain level- a "contour"
%       More than one nonjoint ridge (segment) can exist in a detection.  If this occurs
%           the fields below become horizontal vectors, arranged by
%           descending Area
%
% feature: following fields have Nsegment columns, sorted by descending Area.
%
%                        final_image: {[129x93 logical]  [129x93 logical]}
%                               sparse(squeeze(final_image(:,:,I)));
%                               Cell matrix of an independent image for each segment
%                         Area: Area of segment in pixels.  See
%                                   also time_band_product.
%                    AreaRatio: Ratio of total foreground area to contour
%                               area.
%
%                  BoundingBox: [4xNsegment ]-[Starttime(sec) Maximum_frequency duration  bandwidth]
%                                               of rectangle boxing segment.
%                     Centroid: [2xNsegment]-[Center time Center frequency]
%
%                Centroid_freq: Value of Centroid(2)
%                       Solidity:  the proportion of the pixels in the
%                               convex hull that are also in the region
%                        ctime: start time of segment in c-time.
%                     duration: Duration of segment in seconds.
%                         fend: mean frequency at end of segment (Hz)
%                         fmax: maximum frequency attained by segment (Hz)
%                         fmin: minimum frequency attained by segment (Hz)
%                       fstart: mean frequency at start of segment (Hz)
%             global_bandwidth: total bandwidth spanned by segment
%              local_bandwidth: Maximum instantaneous bandwidth: dF*max(sum(feature.Image,1))
%        median_local_kurtosis: median instantaneous kurtosis
%                               (fourth_moment/local_bandwidth^4)
%             robust_bandwidth: Median instantaneous bandwidth, weighted
%                               by intensity.  Should be robust to pulse
%                               contamination
%                  robust_fmax: maximum frequency attained by mean
%                               frequency .
%                  robust_fmin: minimum frequency attained by mean
%                  frequency.
%            time_band_product: area of segment in sec-Hz (.Area with
%                               units)
%%%%%%%%%
% Contour fields:  All describe complete contour properties associated with
%                   ridge segments. [1 x Nsegment] long, values identical
%                   for each element.
%
%                 Contour_Area: [1xNsegment] Area of dB contour. Same value
%                                   for all elements, same units as Area
%             Contour_duration: Duration of contour (sec)
%                 Contour_fmax: Maximum frequency attained by contour (Hz)
%                 Contour_fmin: Minimum frequency attained by contour (Hz)
%     Contour_global_bandwidth: Total bandwidth traversed by entire segment: BoundingBox(4)
%      Contour_local_bandwidth: Greatest "instantaneous" bandwidth (Hz)\
%%%%%%%%%
%                 Eccentricity: eccentricity of oval fitted over segment.
%                       Extent: Area divided by area of bounding box
%                        Image: []
%              MajorAxisLength: Major axis length of enclosing ellipse in
%                                   pixels
%              MinorAxisLength: Minor axis length in pixels
%                  Orientation: Angle (deg) between segment and time-axis.
%                                   (-90 to 90 degrees)
%%% fields that have one column, regardless of number of segments.  %%%%
%                          SEL: Sound Exposure level of signal (dB re
%                                   1uPa^2-sec)
%                          SNR: Signal to Noise ratio (dB)
%
%                           dF: Height of a pixel (Hz)
%                           dT: Width of a pixel (sec)
%                        power: integrated power spectral density over
%                               global bandwidth
%                        Tcall, Fcall: traced contours of segments, provided if 
%%%%Total fields describe bounding box surrounding all segments in a
%%%%detection.  Thus only a single column.
%             TotalBoundingBox: [4x1 double]
%                Totalduration: 5.120000000000005e-01
%                    Totalfmax: 3.515625000000000e+02
%                    Totalfmin: 2.265625000000000e+02
%
%                      comment: ',plateau,'-describes how ridges in columns
%                                           were linked, via harmonic,
%                                           splice, or contour(plateau).
%
%
%
%
%       final_image: image of shapes, where each pixel for detection i is
%       assigned a value of 'i'.  All segments of a detection share same
%       pixel value.
%
%       Beq:  grayscale image
%       labeled_raw: raw output of bwlabel


% Apr 21, 2008--added an SEL and SNR feature...
% March 4, 2009--added morphological reconstruction arguments
% April 8, 2009--  Added background feature extraction (peaks option)

function [feature,final_image,Beq,labeled_raw]=extract_image_features(y,cstart,param,Idebug)
%global Idebug_global;

Idebug_global=0;  %Hard-wired dipswitch for plotting all output

%Make fake signal for debugging purposes
if 1==0,
    [y,param]=make_fake_signal(param);
    param.morph=rmfield(param.morph,'equalization');
    
end

%%Initialize
feature=[];Beq=[];labeled_raw=[];
%threshold_fudge=param.morph.threshold_fudge;
final_image=[];Bmean=[];

if isempty(y)
    return
end


mydir=pwd;
Fs=param.Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%compute spectrogram, equalize spectrogram, remove buffer part
%%%%                        of image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Create a spectrogram, strip away buffer time (param.morph.eq_time),and restrict frequency
%%  range, then equalize using either data or external spectrum 'eq'
%
%Beq: equalized dB spectrogram, with first param.morph.eq_time s removed
%Borg:   original dB spectrogram, with first param.morph.eq_time s removed
%Bmean:  Estimated background noise spectrum, with same number of
%       frequencies as B and Beq, but expressed in linear (not dB)
%       units.
%T,F:  Time and frequency vector corresponding to spectrogram image.
%      Time has been cropped to remove the first param.morph.eq_time s
% ctime_new:  The ctime of the start of Beq
[Beq,T,F,ctime_new,Bmean,Borg]=equalize_background;
if size(Beq,2)<2,
    feature=[];
    return
end
dT=T(2)-T(1);dF=F(2)-F(1);%dimensions of image "pixels"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Convert image into binary via thresholding, then extract features %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Option 1: Original thresholding method
%param.morph.threshold_chc='otsu';
if strcmp(param.morph.threshold_chc,'local_peaks')
    Bfilt=filter_image(Beq);
    
    %%%%%extract_ridge_trace for foreground (contour) image.  Iopen goes through
    %%%%%bwareaopen after ridge trace
    [Iopen,Back]=extract_ridge_traces(Bfilt,param.morph.SNRmin,param.morph.local_bandwidth.min, ...
        param.morph.local_bandwidth.max,param.morph.dynamic_range,param.morph.dynamic_range2,ceil(param.morph.time_band_product.min/(dT*dF)), ...
        param.Nfft,param.Fs,F,Idebug);
    
    %Close small gaps in foreground (contour) image, using continous close
    [Iridge,SE,debug_params_more]=close_image(Iopen, param.morph.gap_f,param.morph.gap_t, dF, dT);
    
    %%Close small gaps in background image, where Back>SNRmin, using binary
    %%  closing
    
    
    Icontour=close_image(logical(Back), param.morph.background.gap_f,param.morph.background.gap_t, dF, dT);
    Icontour=flipud(Icontour);
    
    [feature0,labeled_contour,numObjects]=extract_features_from_binary_image(Iridge,T,F,ctime_new,Idebug);
    feature=extract_features_from_grayscale_image(feature0,Beq,labeled_contour,numObjects,T,F,Idebug,param.morph.want_contour);
    
    plot_dual_threshold;
elseif strcmp(param.morph.threshold_chc,'otsu')
    %Convert to gray level-plots, using SNRmin as the cutoff
    
    %Optional Filter image
    Bgray=mat2gray(Beq,[param.morph.SNRmin max(max(Beq))]);
    Bgray=flipud(Bgray);
    Bfilt=filter_image(Bgray);
    
    %%%%%%%%%%%%%%Convert to binary%%%%%%%%%%%%%%%%%%%%
    level=graythresh(Bfilt);
    Ibin=im2bw(Bfilt,threshold_fudge*level);
    
    %%%%%%%%Open image to clear out small objects%%%%%%
    Iopen=bwareaopen(Ibin,ceil(param.morph.time_band_product.min/(dT*dF)));
    [Iclose,SE,debug_params_more]=close_image(Iopen, param.morph.gap_f,param.morph.gap_t, dF, dT);
    [feature0,labeled_contour,numObjects]=extract_features_from_binary_image(Iclose,T,F,ctime_new,Idebug);
    feature=extract_features_from_grayscale_image(feature0,Beq,labeled_raw,numObjects,T,F,Idebug);
    
    
elseif strcmp(param.morph.threshold_chc,'reconstruction')
    Bfilt=(filter_image(Beq));
    [Iopen, Bgray,debug_rec]=threshold_image_by_reconstruction(Bfilt,param,dT,dF,Idebug);
    Iopen=flipud(Iopen);
    Bfilt=flipud(Bgray);
    if sum(sum(Bgray))==0, return;end
    [Iclose,SE,debug_params_more]=close_image(Iopen, param.morph.gap_f,param.morph.gap_t, dF, dT);
    [feature0,labeled_contour,numObjects]=extract_features_from_binary_image(Iclose,T,F,ctime_new,Idebug);
    feature=extract_features_from_grayscale_image(feature0,Bfilt,labeled_contour,numObjects,T,F,Idebug);
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Run feature filters%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(feature),
    feature=[];
    if Idebug>1
        disp('no contour peaks traced');
    end
    return
end

Ipass=crude_filter_feature(feature,param,Idebug);
if isempty(Ipass),
    feature=[];
    final_image=(zeros(size(Bfilt)));
    
    if Idebug>1
        disp('extract_image_features: 220: no features pass');
        if strcmp(param.morph.threshold_chc,'otsu')
            plot_morph_steps_otsu;
        elseif strcmp(param.morph.threshold_chc,'local_peaks')
            plot_morph_steps_peaks('after filter feature');
        elseif strcmp(param.morph.threshold_chc,'reconstruction')
            plot_morph_steps_reconstruct;
        end
    end
    return
end
feature=feature(Ipass);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Create final image, where values of pixels are integers...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_final_image=(zeros(size(Bfilt)));
for KK=1:length(Ipass),
    initial_final_image=initial_final_image+KK*ismember(labeled_contour,Ipass(KK));  %Image now indexed relative to the sorted feature
end

if Idebug>1
    figure;
    subplot(2,1,1)
    imagesc(T,F,Iridge);
    title('before feature filtering')
    
    subplot(2,1,2)
    imagesc(T,F,initial_final_image);
    title('after feature filtering');
    %keyboard
    
end

%%%%Attempt to Merge features that are overlapping (harmonic) or close
%%%%  together--this saves time on the reconstruction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[feature,final_image]=link_features(feature,param,initial_final_image,Idebug);

%Ifinal=Ipass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Add background features to merged detections       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(param.morph.threshold_chc,'local_peaks')& param.morph.background.on==1
    %Icontour=extract_background_image(Beq,SE,param,dF,dT,T,F,Idebug);
    
    [feature_background,labeled_background]=extract_features_from_binary_image(Icontour,T,F,ctime_new,Idebug);
    %[feature_background,labeled_background]=link_features(feature_background,param,labeled_background,Idebug);
    
    %labeled_background(labeled_background==0)=Inf;
    
    %%Cycle through final features and images...
    [feature,final_image]=merge_ridge_with_contour(feature,final_image,Icontour,feature_background,labeled_background,T,F);
    
    
    
end

%Create complete image and global statistics
%BWfinal=ismember(labeled_contour,Ifinal);
%[value,Imax]=max(sum(BWfinal));
%feature(1).vertical_percent=value/size(BWfinal,1);
for I=1:size(final_image,3),
    feature(1).final_image{I}=sparse(squeeze(final_image(:,:,I)));
end
feature(1).dF=dF;
feature(1).dT=dT;
feature(1).T=T;
feature(1).F=F;

%Compute SNR and SEL
for JJ=1:length(feature),
    [feature(JJ).SNR,feature(JJ).SEL,feature(JJ).power]=compute_SNR(JJ);
end

if Idebug>=2||Idebug_global==1,
    if strcmp(param.morph.threshold_chc,'otsu')
        plot_morph_steps_otsu;
    elseif strcmp(param.morph.threshold_chc,'local_peaks')
        plot_morph_steps_peaks('final');
    elseif strcmp(param.morph.threshold_chc,'reconstruction')
        plot_morph_steps_reconstruct;
    end
    %close all;
    
    
    % keyboard;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%InnerFunctions%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%merge_ridge_with_contour.m%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [feature,final_image,Iunique]=merge_ridge_with_contour(feature,final_image,Icontour,feature_background,labeled_background,T,F)
        
        back_index=-1*ones(length(feature),1);
        for II=1:length(feature)  %fore each foreground feature...
            
            Iseed=1;
            seed=squeeze(final_image(:,:,II));
            while Iseed
                try
                    Iobr = imreconstruct(seed, Icontour);
                    Iseed=0;
                catch
                    % disp('reconstruct failed')
                    seed=logical(imerode(seed,SE));
                    if isempty(seed)
                        Iseed=-1;
                    end
                end
            end
            
            if Iseed==0
                
                test_image=Iobr.*labeled_background;
                Itest_unique=unique(test_image(:));  %Produces a vector [0 index1 index2...];
                Nunique=length(Itest_unique)-1;
                
                if Idebug>2
                    Fup=flipud(F);
                    figure(3);set(gcf,'pos',[1342   98         554         981]);
                    subplot(4,1,1)
                    imagesc(T,Fup,squeeze(final_image(:,:,II)));axis('xy');
                    title('original foreground image');
                    subplot(4,1,2)
                    imagesc(T,Fup,5*seed+Icontour);title('seed+background');axis('xy');
                    subplot(4,1,3)
                    try
                        imagesc(T,Fup,Iobr);axis('xy');title(sprintf('reconstructed background, matching index %i',back_index(II)));
                    catch
                        %disp
                    end
                    subplot(4,1,4)
                    imagesc(T,flipud(F),test_image);axis('xy');
                    title(sprintf('There are %i components that need to be merged...',Nunique));
                    
                    if Nunique>1
                        disp('Original base background feature:')
                        feature_background(Itest_unique(2))
                        
                    end
                    pause;
                end
                
                
                %%Now for a complex situation...occaisionally the
                %%"merge_features" result on foreground "ridge" analysis
                %% creates a seed that reconstructs two disjoint background
                %% contours.  Thus these contour properties need to be
                %% merged...
                if Nunique>1  %more than zero and one
                    %%Two separate background blobs are associated with
                    %%  each other.  Merge the background blob and
                    %%  features.
                    %disp('Disjoint background blobs linked via foreground');
                    for Ifeature=1:(Nunique-1)  %Elininate other contour indicies other than the first one
                        
                        feature_background(Itest_unique(2))=merge_features(feature_background(Itest_unique(2)),feature_background(Itest_unique(2+Ifeature)));
                        %Review past ridge/contour links and update
                        %  contour index as needed
                        Ireplace=find(back_index==Itest_unique(Ifeature+2));
                        back_index(Ireplace)=Itest_unique(2);
                        if Idebug>1
                            disp(sprintf('Feature to be merge #%i, contour blob %i',Ifeature,Itest_unique(2+Ifeature)));
                            feature_background(Itest_unique(2+Ifeature))
                            feature_background(Itest_unique(2+Ifeature)).Centroid
                            feature_background(Itest_unique(2+Ifeature)).BoundingBox
                            
                            disp('Merged results')
                            feature_background(Itest_unique(2))
                            feature_background(Itest_unique(2)).Centroid
                            feature_background(Itest_unique(2)).BoundingBox
                            
                            
                        end
                    end
                    back_index(II)=Itest_unique(2);
                    
                else
                    back_index(II)=max(Itest_unique);
                end
                %Flag which background object
                if Idebug>2
                    subplot(4,1,3)
                    title(sprintf('reconstructed background, matching index %i',back_index(II)));
                    pause
                end
                
                
            else
                disp('No seed left');
                back_index(II)=0;
            end
            
        end
        
        
        %%Merge foreground ridges together that share a contour...
        
        [Iunique,Ipass]=unique(back_index,'first');  %Iunique are indicies to background contours
        % Ipass are indicies to foreground features
        names=fieldnames(feature);
        for I=1:length(Ipass)
            Iwant=Ipass(I);
            Isame=find(Iunique(I)==back_index);  %feature indicies that share same back index
            
            for J=2:length(Isame);
                index=(Isame(J));
                %feature(Iwant)=merge_features(feature(Iwant),feature(index));
                feature(Iwant)=clump_features(feature(Iwant),feature(index),names,'plateau');
                final_image(:,:,Iwant)=final_image(:,:,Iwant)+final_image(:,:,index);
                
            end
            
            
        end
        feature=feature(Ipass);
        final_image=final_image(:,:,Ipass);
        
        %Merge background contour features into feature vector
        background_names={'Area','global_bandwidth','duration','local_bandwidth','fmin','fmax'};
        
        for I=1:length(feature)
            [feature(I).TotalBoundingBox,feature(I).Totalfmax,feature(I).Totalfmin,feature(I).Totalduration]=merge_BoundingBox(feature(I).BoundingBox);
            
            %%Now add background features to the mix
            Nparts=length(feature(I).Area);
            for J=1:length(background_names)
                try
                    feature(I).(sprintf('Contour_%s',background_names{J}))=(feature_background(Iunique(I)).(background_names{J}))*ones(1,Nparts);
                catch
                    %keyboard;
                    %disp('No contour segment match');
                    feature(I).(sprintf('Contour_%s',background_names{J})) =zeros(1,Nparts);
                end
            end
            feature(I).AreaRatio=feature(I).Contour_Area./sum(feature(I).Area);
            % keyboard;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%filter_image.m%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Imed=filter_image(Bgray)
        
        if param.filter.on==0
            Imed=Bgray;
            return
        end
        %Optional median filter
        if isfield(param,'median')&&param.median.on==1,
            min_duration=ceil(param.median.size(1)/(dT));
            min_freq=ceil(param.median.size(2)/(dF));
            
            Imed=medfilt2(Bgray,[min_freq min_duration],'symmetric');
        else
            Imed=Bgray;
            
        end
        %%Optional Gaussian filtering
        if isfield(param,'filter')&&param.filter.on==2,    %apply gaussian filter to B
            %disp('applying gaussian filter');
            
            min_duration=ceil(param.filter.size(1)/dT);
            min_freq=ceil(param.filter.size(2)/dF);
            
            Hgauss = fspecial('gaussian',[min_freq min_duration],param.filter.sigma);
            %Hgauss = fspecial('average',min_freq);
            Ifilt=imfilter(Imed,Hgauss,'symmetric');
            if Idebug==2,
                figure(3);
                subplot(3,1,1);
                imagesc(T,F,Imed);
                title('Original spectrogram');
                caxis([40 120]); %colorbar;
                axis('xy');
                
                subplot(3,1,2);
                imagesc(T,F,Ifilt);
                title(sprintf('Gaussian filter with  size %s, sigma %6.2f',mat2str(param.filter.size),param.filter.sigma));
                caxis([0 40]); %colorbar;
                axis('xy');
                
                subplot(3,1,3);
                imagesc(T,F,Ifilt-Imed);
                title(sprintf('Gaussian filter with  size %s, sigma %6.2f',mat2str(param.filter.size),param.filter.sigma));
                caxis([-10 10]); %colorbar;
                axis('xy');
                
            end
            Imed=Ifilt;
        elseif isfield(param,'filter')&&param.filter.on==1,    %apply assymetric gaussian filters
            
            
            %disp('applying gaussian filter');
            
            min_duration=ceil(param.filter.size(1)/dT);
            min_freq=ceil(param.filter.size(2)/dF);
            
            min_duration=10;
            min_freq=10;
            param.filter.sigma=1;
            
            Xfilt=zeros(size(Bgray,1),size(Bgray,2),4);
            for III=1:4
                switch III
                    case {1,3}
                        boxx=3*[1 6];
                    case {2,4}
                        boxx=3*[6 1];
                end
                
                Hgauss{III}=gaussian_asym_filter([min_freq min_duration],param.filter.sigma*boxx,III);
                Xfilt(:,:,III)=imfilter(Imed,Hgauss{III},'symmetric');
            end
            
            
            re=(max(Xfilt(:,:,2:4),[],3));
            coef=[0.5 0.5];
            Ifilt=coef(1)*Imed+coef(2)*(re-Xfilt(:,:,1));
            
            if Idebug==2
                figure(3)
                subplot(5,1,1)
                imagesc(T,F,Imed);
                title('Original spectrogram');
                %caxis([0 40]); %colorbar;
                axis('xy');
                
                for III=1:4
                    subplot(5,1,III+1)
                    imagesc(T,F,squeeze(Xfilt(:,:,III)))
                    axis('xy');
                end
                
                figure(4)
                subplot(4,1,1)
                imagesc(T,F,Imed);colorbar
                title('Original spectrogram');%caxis([0 40])
                %caxis([0 40]); %colorbar;
                axis('ij');
                
                subplot(4,1,2);
                imagesc(T,F,squeeze(Xfilt(:,:,1)));caxis([0 40]);colorbar
                title('vertical filter');
                
                subplot(4,1,3)
                imagesc(T,F,re);caxis([0 40]);colorbar;axis('ij')
                title('re--maximum non-vertical dimension')
                
                subplot(4,1,4)
                imagesc(T,F,Ifilt);caxis([0 40]);colorbar;axis('ij')
                
            end
            
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%equalize_background.m%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Create a spectrogram, strip away buffer time (param.morph.eq_time),and restrict frequency
%%  range, then equalize using either data or external spectrum 'eq'
    function [Beq,T,F,ctime_new,Bmean,B]=equalize_background
        %Beq: equalized dB spectrogram, with first param.morph.eq_time s removed
        %B:   original dB spectrogram, with first param.morph.eq_time s removed
        %Bmean:  Estimated background noise spectrum, with same number of
        %       frequencies as B and Beq, but expressed in linear (not dB)
        %       units.
        %T,F:  Time and frequency vector corresponding to spectrogram image.
        %      Time has been cropped to remove the first param.morph.eq_time s
        % ctime_new:  The ctime of the start of Beq
        
        %%PSD is power spectral density
        [S,F,T,PSD]=spectrogram(y,param.Nfft,round(param.ovlap*param.Nfft),param.Nfft,Fs,'yaxis');
        %imagesc(T,F,10*log10(abs(PSD)));
        
        %%%Crop image...%%%%%%
        [junk,IfreqMin]=min(abs(F-param.morph.MinFreq));
        [junk,IfreqMax]=min(abs(F-param.morph.MaxFreq));
        PSD=PSD(IfreqMin:IfreqMax,:);
        F=F(IfreqMin:IfreqMax);
        B=10*log10(abs(PSD+eps));  %Power spectral density..
        
        %%%%%Find median noise background level
        Imedian_sample=find(T<=param.morph.eq_time);
        if length(Imedian_sample)>1,
            Bmean1=median(B(:,Imedian_sample).').';  %Note that the median is unaffected by a log transformation of the data
        else
            Bmean1=B(:,Imedian_sample);
        end
        
        %If external background spectrum provided in [Nfreq by 2] matrix
        % use these interpolated values over the frequency range of the extrnal
        % spectrum
        if isfield(param.morph,'equalization') & any(isinf(param.morph.equalization(:,2)))==0
            eq_coarse=param.morph.equalization(:,2);
            f_coarse=param.morph.equalization(:,1);
            
            if Idebug==3,
                figure(2)
                subplot(4,1,1);
                imagesc(T,F,B);axis('xy');title(ctime2str(cstart));
                subplot(4,1,2);
                plot(f_coarse,eq_coarse);
                hold off;
            end
            
            
            
            Bmean=Bmean1;
            Iedit=find(F>=min(f_coarse)&(F<=max(f_coarse)));
            Bmean(Iedit)=interp1(f_coarse,eq_coarse,F(Iedit));
            
            if Idebug==3,
                subplot(4,1,3);
                figure(2)
                plot(f_coarse,eq_coarse,'r',F,Bmean,'k--',F,Bmean1,'g');grid on;legend('stored','estimated (+5 dB shift)','short term');ylim([40 90]);
                hold off;
            end
            
        else
            Bmean=Bmean1;
            %[Bout,fout]=crude_decimate(Bmean,F,1000,'median');
            
            if isfield(param.morph,'equalization')&&any(isinf(param.morph.equalization(:,2)))==0
                fprintf('Crash: extract_image_features equalization has inf at %s\n',ctime2str(cstart));
            end
            
        end
        
        %%Remove buffer time from image (but don't alter time series)
        Beq=B-Bmean*ones(1,size(B,2));
        Imax=max(Imedian_sample);
        Beq=Beq(:,Imax:end);
        B=B(:,Imax:end);
        ctime_new=cstart+T(Imax);
        T=T(Imax:end);
        %ctime_new=cstart;
        
        if Idebug==3,
            subplot(4,1,4)
            imagesc(T,F,Beq);axis('xy');caxis([0 30]);
            pause;
        end
        %Convert Bmean into linear units
        Bmean=10.^(Bmean/10);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%make_fake_signal%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [y,param]=make_fake_signal(param)
        Idebug=1;
        param.morph.eq_time=0.45;
        param.Fs=1000;
        param.Nfft=128;
        fstart=50;
        fend=75;
        tduration=0.75;
        tstart=1+param.morph.eq_time;
        total_window_time=5;
        SNRdB=20;
        ovlap=0.75;
        Nfft=2^nextpow2(tduration*param.Fs);
        SNRchc='PSD';
        disp('FM Sweep');
        [y,ysignal,ynoise]=simulate_FMsweep_with_noise(param.Fs,param.Nfft,ovlap,fstart,fend, ...
            total_window_time,tduration,tstart,SNRdB,SNRchc);
        
        disp(' ');disp('Signal properties');
        [SELsignal,power_signal]=extract_SEL_n_PSD(y(round(tstart*param.Fs)+(1:round(tduration*param.Fs))),Nfft,param.Fs,[fstart fend]);
        % title('signal');
        
        disp(' ');disp('Noise properties');
        [SELnoise,power_noise]=extract_SEL_n_PSD(y(1:round(tstart*param.Fs)),Nfft,param.Fs,[fstart fend]);
        %title('noise');
        
        %%%SNR=10*log10(power_signal.narrowband/power_noise.narrowband)
        
        [y2,ysignal,ynoise]=simulate_FMsweep_with_noise(param.Fs,param.Nfft,ovlap,fstart,fend, ...
            total_window_time,tduration,tstart+2*tduration,SNRdB/2,SNRchc);
        y=y+0.5*y2;
        
        % keyboard;
        param.morph.MinFreq=50;param.morph.MaxFreq=300;
        
        
    end


%%%%%%%%%%%%%%%%%%%%compute_SNR.m%%%%%%%%%%%%%%%%%%%%
    function [SNR,SEL,pwr,rms]=compute_SNR(I)
        Iy1=max([1 round((min(feature(I).ctime)-cstart)*param.Fs)]);
        Id=round(feature(I).TotalBoundingBox(3)*param.Fs); %signal duration in samples
        In=round(param.morph.eq_time*param.Fs);  %desired noise duration in samples
        snr_Nfft=2^nextpow2(max([Id In]));
        
        % Then, compute noise power by taking sample from start of time
        % series..
        %[noiseSEL,noise_power]=extract_SEL_n_PSD(y(1:In),snr_Nfft,param.Fs,[feature(I).fmin feature(I).fmax]);
        %noise_rmsdB=10*log10(noise_power.narrowband);
        
        %Second estimate of noise power from Bmean (PSD curve)
        [junk,Imin]=min(abs(F-feature(I).Totalfmin));
        [junk,Imax]=min(abs(F-feature(I).Totalfmax));
        noise_rmsdB=10*log10(dF*trapz(Bmean(Imin:Imax)));
        
        %disp(sprintf('Original noise estimate: %6.2f dB, new noise estimate: %6.2f dB, diff %6.2f', ...
        %		noise_rmsdB,long_term.power,noise_rmsdB-long_term.power));
        %Extract signal power
        try
            [temp.SEL,temp.power]=extract_SEL_n_PSD(y(Iy1+(1:Id)-1),snr_Nfft,param.Fs,[feature(I).Totalfmin feature(I).Totalfmax]);
        catch
            disp('postprocess_morph, error loading extract_SEL_N_PSD');
            keyboard;
        end
        SEL=10*log10(abs(temp.SEL.narrowband));
        
        pwr=temp.power.narrowband;
        rms=10*log10(pwr);
        SNR=rms-noise_rmsdB;
        
        if Idebug==4,
            figure(10);set(gcf,'units','norm','pos',[4.5312e-02   8.3333e-02   2.9167e-01   3.5000e-01]);
            Btemp=ismember(labeled_contour,Ifinal(I));
            
            Ion=find(flipud(Btemp)==1);
            Bfilt2=Bfilt;
            
            Bfilt2(Ion)=-20;
            subplot(2,1,1);imagesc(T,F,Bfilt);caxis([-20 30]);colorbar;axis('xy');
            title('Equalized spectrogram of signal');
            subplot(2,1,2);imagesc(T,F,Bfilt2);caxis([-20 30]);colorbar;axis('xy');
            title(sprintf(' SNR: %6.2f dB, SEL: %6.2f fmin: %6.2f fmax: %6.2f',SNR,SEL,feature(I).Totalfmin,feature(I).Totalfmax));
            
            figure(20);set(gcf,'units','norm','pos',[4.7917e-02   5.0083e-01   2.9167e-01   3.5000e-01]);
            subplot(2,1,1);
            spectrogram(y(1:round(param.morph.eq_time*param.Fs)),param.Nfft,round(param.ovlap*param.Nfft),param.Nfft,Fs,'yaxis')
            line([0 param.morph.eq_time],feature(I).Totalfmin *[1 1],'color','k','linewidth',3);
            line([0 param.morph.eq_time],feature(I).Totalfmax *[1 1],'color','k','linewidth',3);
            
            colorbar;axis('xy');title('spectrogram of noise sample used to compute SNR');
            
            subplot(2,1,2);
            spectrogram(y(Iy1+(1:Id)-1),param.Nfft,round(param.ovlap*param.Nfft),param.Nfft,Fs,'yaxis')
            line([0 feature(I).TotalBoundingBox(3)],feature(I).Totalfmin *[1 1],'color','k','linewidth',3);
            line([0 feature(I).TotalBoundingBox(3)],feature(I).Totalfmax *[1 1],'color','k','linewidth',3);
            
            colorbar;axis('xy');title(sprintf('spectrogram of signal sample used to compute SNR, %6.2f seconds',feature(I).TotalBoundingBox(3)));
            pause;
        end
        
        %%%%%%%%end SNR%%%%%%%%%%%%%%%%%%%
        
    end

    function plot_dual_threshold
        if Idebug>1
            figure(1);
            set(1,'pos',[86   214   1437   1030]);
            subplot(3,2,1)
            imagesc(T,F/1000,Borg);%caxis([0 max(max(Beq))]);axis('xy');
            set(gca,'fontweight','bold','fontsize',14);%xlabel('Time');
            ylabel('kHz');
            title(sprintf('Original spectrogram, time %s',ctime2str(cstart) ));
            poss=get(gca,'pos');poss(3)=0.37;poss(4)=0.25;set(gca,'pos',poss);
            
            subplot(3,2,2)
            imagesc(T,F/1000,Beq);
            %caxis([0 max(max(Beq))]);
            axis('xy');
            set(gca,'fontweight','bold','fontsize',14);%xlabel('Time');
            ylabel('kHz');
            title(sprintf('Original equalized image Beq, time %s',ctime2str(cstart) ));
            poss=get(gca,'pos');poss(3)=0.37;poss(4)=0.25;set(gca,'pos',poss);
            
            subplot(3,2,3)
            imagesc(T,F/1000,flipud(Iopen));axis('xy');%xlabel('Time');
            set(gca,'fontweight','bold','fontsize',14);%xlabel('Time');
            ylabel('kHz');
            set(gca,'fontweight','bold','fontsize',14);
            poss=get(gca,'pos');poss(3)=0.37;poss(4)=0.25;set(gca,'pos',poss);
            
            title(sprintf('Opened ridge-traced image (Iopen), SNRmin %6.2f and dynamic ranges %6.2f and %6.2f',param.morph.SNRmin, param.morph.dynamic_range, ...
                param.morph.dynamic_range2));
            
            subplot(3,2,4)
            imagesc(T,F/1000,Back);axis('xy')
            set(gca,'fontweight','bold','fontsize',14);%xlabel('Time');
            ylabel('kHz');
            poss=get(gca,'pos');poss(3)=0.37;poss(4)=0.25;set(gca,'pos',poss);
            
            title(sprintf('Contour (Back): Region above %6.2f dB SNR',param.morph.SNRmin));
            
            subplot(3,2,5)
            imagesc(T,F/1000,flipud(Iridge));axis('xy')
            set(gca,'fontweight','bold','fontsize',14);xlabel('Time (s)');ylabel('kHz');
            poss=get(gca,'pos');poss(3)=0.37;poss(4)=0.25;set(gca,'pos',poss);
            
            title(sprintf('Ridge traced image (Iridge), closed using dt %6.2f and df %6.2f',param.morph.gap_t, param.morph.gap_f));
            
            
            subplot(3,2,6)
            imagesc(T,F/1000,flipud(Icontour));axis('xy');
            set(gca,'fontweight','bold','fontsize',14);;xlabel('Time (s)');ylabel('kHz');
            poss=get(gca,'pos');poss(3)=0.37;poss(4)=0.25;set(gca,'pos',poss);
            
            title(sprintf('Contour(Icontour) after morph closing'));
            
            text(0,0.7,'plot_dual_threshold');
            
            %keyboard
            
        end
    end

%%%%%%%%%%%%%%%%%%plot_morph_steps_reconstruct.m%%%%%%%%%%%%%%%%%%%%%
    function plot_morph_steps_reconstruct
        figure(11);
        set(gcf,'pos',[ 58         160        1153         905]);
        %subplot(3,2,1);
        %imagesc(T,F,Bgray);axis('ij');title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
        %colormap(flipud(gray));
        subplot(3,3,1);
        imagesc(T,(F),(Beq));axis('xy');title('equalized spectrogram')
        colorbar;
        subplot(3,3,2);
        imagesc(T,(F),(debug_rec.Beq));axis('xy');title('equalized spectrogram')
        colorbar;
        
        
        subplot(3,3,3);
        imagesc(T,(F),debug_rec.Iobr);axis('xy'); colorbar;
        try
            title('Iobrh');
            title(sprintf('Opened image, element size %6.2f',debug_rec.se_radius));
            %title(sprintf('Median filter size %i freq by %i time pixels',min_freq,min_duration));
        catch
            %title('nonlinear peak filter');
        end
        subplot(3,3,4);
        imagesc(T,(F),debug_rec.Ih);axis('xy');
        title(sprintf('h-dome, %6.2f dB dynamic range',param.morph.dynamic_range));
        colorbar;
        
        subplot(3,3,5);
        imagesc(T,(F),debug_rec.Ip);axis('xy');
        title(sprintf('Original opened image minus h-dome'));
        colorbar;
        subplot(3,3,6);
        imagesc(T,(F),debug_rec.Ih2);axis('xy');
        title(sprintf('regions associated with peaks %6.2f dB or above',param.morph.SNRmin));
        colorbar;
        subplot(3,3,7);
        imagesc(T,flipud(F),Iopen);axis('xy');
        title('Binary conversion');
        colorbar;
        subplot(3,3,8);
        imagesc(T,flipud(F),Iclose);axis('xy');title(sprintf('imclose, struct element: %s %6.3f Hz by %6.3f sec', ...
            strel_name,dF*strel_size(1),dT*strel_size(2)));
        colorbar;
        
        try
            subplot(3,3,9);
            %imagesc(T,F,BWfinal);
            imagesc(T,flipud(F),sum(final_image,3));axis('xy');hold on;
            for JJJ=1:length(feature),
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1)+[0 feature(JJJ).Totalduration],[feature(JJJ).Totalfmin feature(JJJ).Totalfmin]);set(hh,'Color',[1 0 1]);
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + [0 feature(JJJ).Totalduration],[feature(JJJ).Totalfmax feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + [0 0],[feature(JJJ).Totalfmin feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + feature(JJJ).Totalduration*[1 1],[feature(JJJ).Totalfmin feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
                
            end
            
            figure(12);set(gcf,'units','norm','pos',[  6.968750000e-01     1.058333e-01     3.03120e-01     7.458333e-01]);
            
            subplot(2,1,1);
            %imshow(Bgray,'xdata',T,'ydata',F);axis('xy');
            imagesc(T,(F),(Beq));axis('xy');%colormap(flipud(gray));%caxis([-20 50]);
            
            title(sprintf('original equalized image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
            format long e
            disp(cstart);
            %image(T,F,B,colormap('jet'));
            xlabel('Time (sec)');ylabel('Frequency (Hz)');
            subplot(2,1,2);
            imagesc(T,flipud(F),sum(final_image,3));axis('xy');title('Final image');%colormap('jet');
            xlimm=xlim;
            %subplot(3,1,3)
            %BWskel = bwmorph(BWfinal,'thin','inf');
            %for I=1:length(feature),
            %   plot(feature(I).Tfit,feature(I).Ffit,'ko');hold on; title('Fitted contours');
            
            %end
            hold off;
            ylim([0 500]);
            xlim(xlimm);
            %imagesc(T,flipud(F),BWskel);axis('xy');title('thin');%colormap('jet');
            close(11:12);
        end
    end

%%%%%%%%%%%%%%%%%%plot_morph_steps_otsu.m%%%%%%%%%%%%%%%%%%%%%
    function plot_morph_steps_otsu
        figure(11);
        set(gcf,'pos',[680   128   603   906]);
        %subplot(3,2,1);
        %imagesc(T,F,Bgray);axis('ij');title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
        %colormap(flipud(gray));
        subplot(3,2,1);
        imagesc(T,(F),flipud(Bgray));axis('xy');title('grayscale image')
        
        subplot(3,2,2);
        imagesc(T,flipud(F),Bfilt);axis('xy');
        try
            title('Bfilt');
            %title(sprintf('Median filter size %i freq by %i time pixels',min_freq,min_duration));
        catch
            %title('nonlinear peak filter');
        end
        subplot(3,2,3);
        imagesc(T,flipud(F),Ibin);axis('xy');
        title('binary conversion');
        
        subplot(3,2,4);
        imagesc(T,flipud(F),Iopen);axis('xy');
        title(sprintf('Morphological opening minimum area %6.2f Hz-s',param.morph.time_band_product.min));
        
        subplot(3,2,5);
        imagesc(T,flipud(F),Iclose);axis('xy');title(sprintf('imclose, struct element: %s %6.3f Hz by %6.3f sec', ...
            debug_params_more.strel_name,dF*debug_params_more.strel_size(1),dT*debug_params_more.strel_size(2)));
        
        %subplot(3,2,5);
        %imagesc(T,flipud(F),Ierode);axis('xy');title(sprintf('ierode, struct element: %s %i freq pixel by %i time pixel', ...
        %   strel_name,dF*strel_size(1),dT*strel_size(2)));
        
        %     hold on;
        %     colors=['b' 'g' 'r' 'c' 'm' 'y'];
        %     for k=1:length(bounds),
        %         boundary=bounds{k};
        %         cIpass=mod(k,length(colors))+1;
        %         plot(boundary(:,2),boundary(:,1),...
        %             colors(cIpass),'Linewidth',2);
        %     end
        %     hold off
        try
            subplot(3,2,6);
            %imagesc(T,F,BWfinal);
            imagesc(T,flipud(F),sum(final_image,3));axis('xy');hold on;
            for JJJ=1:length(feature),
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1)+[0 feature(JJJ).Totalduration],[feature(JJJ).Totalfmin feature(JJJ).Totalfmin]);set(hh,'Color',[1 0 1]);
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + [0 feature(JJJ).Totalduration],[feature(JJJ).Totalfmax feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + [0 0],[feature(JJJ).Totalfmin feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
                hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + feature(JJJ).Totalduration*[1 1],[feature(JJJ).Totalfmin feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
                
            end
            
            figure(12);set(gcf,'units','norm','pos',[  6.968750000e-01     1.058333e-01     3.03120e-01     7.458333e-01]);
            
            subplot(2,1,1);
            %imshow(Bgray,'xdata',T,'ydata',F);axis('xy');
            imagesc(T,(F),flipud(Bgray));axis('xy');colormap(flipud(gray));%caxis([-20 50]);
            
            title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
            format long e
            disp(cstart);
            %image(T,F,B,colormap('jet'));
            xlabel('Time (sec)');ylabel('Frequency (Hz)');
            subplot(2,1,2);
            imagesc(T,flipud(F),sum(final_image,3));axis('xy');title('Final image');%colormap('jet');
            xlimm=xlim;
            %subplot(3,1,3)
            %BWskel = bwmorph(BWfinal,'thin','inf');
            %for I=1:length(feature),
            %   plot(feature(I).Tfit,feature(I).Ffit,'ko');hold on; title('Fitted contours');
            
            %end
            hold off;
            ylim([0 500]);
            xlim(xlimm);
            %imagesc(T,flipud(F),BWskel);axis('xy');title('thin');%colormap('jet');
            keyboard;
            close(11:12);
        end
    end

%%%%%%%%%%%%%%%%%%plot_morph_steps_peaks.m%%%%%%%%%%%%%%%%%%%%%
    function plot_morph_steps_peaks(comment)
        
        
        figure(11);
        set(gcf,'pos',[680   128   603   906]);
        %subplot(3,2,1);
        %imagesc(T,F,Bgray);axis('ij');title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
        %colormap(flipud(gray));
        %         subplot(4,1,1);
        %         imagesc(T,(F),(Beq));axis('xy');title('grayscale image')
        
        subplot(1,5,1);
        %imagesc(T,flipud(F),flipud(Bfilt));axis('xy');
        imagesc(T,(F),Borg);axis('xy');
        set(gca,'fontweight','bold','fontsize',14);
        try
            title('Original spectrogram');
            %title(sprintf('Median filter size %i freq by %i time pixels',min_freq,min_duration));
        catch
            %title('nonlinear peak filter');
        end
        subplot(1,5,2);
        imagesc(T,(F),flipud(Iopen));axis('xy');
        set(gca,'fontweight','bold','fontsize',14);
        
        title(sprintf('Morphological opening minimum area %6.2f Hz-s',param.morph.time_band_product.min));
        
        subplot(1,5,3);
        imagesc(T,(F),flipud(Iridge));axis('xy');
        set(gca,'fontweight','bold','fontsize',14);
        
        title(sprintf('final closed Iridge, struct element: %s %6.3f Hz by %6.3f sec', ...
            debug_params_more.strel_name,dF*debug_params_more.strel_size(1),dT*debug_params_more.strel_size(2)));
        
        %         subplot(3,2,5);
        %         imagesc(T,flipud(F),Iclose);axis('xy');title(sprintf('imclose, struct element: %s %6.3f Hz by %6.3f sec', ...
        %             debug_params_more.strel_name,dF*debug_params_more.strel_size(1),dT*debug_params_more.strel_size(2)));
        
        %subplot(3,2,5);
        %imagesc(T,flipud(F),Ierode);axis('xy');title(sprintf('ierode, struct element: %s %i freq pixel by %i time pixel', ...
        %   strel_name,dF*strel_size(1),dT*strel_size(2)));
        
        %     hold on;
        %     colors=['b' 'g' 'r' 'c' 'm' 'y'];
        %     for k=1:length(bounds),
        %         boundary=bounds{k};
        %         cIpass=mod(k,length(colors))+1;
        %         plot(boundary(:,2),boundary(:,1),...
        %             colors(cIpass),'Linewidth',2);
        %     end
        %     hold off
        
        subplot(1,5,4);
        %imagesc(T,F,BWfinal);
        title('after feature filtering');
        imagesc(T,(F),flipud(initial_final_image));axis('xy');
        set(gca,'fontweight','bold','fontsize',14);
        
        subplot(1,5,5);
        %imagesc(T,F,BWfinal);
        for Iim=1:size(final_image,3)
            temp(:,:,Iim)=Iim*final_image(:,:,Iim);
        end
        imagesc(T,(F),flipud(sum(temp,3)));axis('xy');
        set(gca,'fontweight','bold','fontsize',14);
        
        hold on;
        for JJJ=1:length(feature),
            hh=line(min(T)+feature(JJJ).TotalBoundingBox(1)+[0 feature(JJJ).Totalduration],[feature(JJJ).Totalfmin feature(JJJ).Totalfmin]);set(hh,'Color',[1 0 1]);
            hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + [0 feature(JJJ).Totalduration],[feature(JJJ).Totalfmax feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
            hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + [0 0],[feature(JJJ).Totalfmin feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
            hh=line(min(T)+feature(JJJ).TotalBoundingBox(1) + feature(JJJ).Totalduration*[1 1],[feature(JJJ).Totalfmin feature(JJJ).Totalfmax]);set(hh,'Color',[1 0 1]);
            
        end
        hold off;
        title('feature filtered and merged');
        for LL=1:length(feature)
            disp(feature(LL))
            
        end
        try
            figure(12);set(gcf,'units','norm','pos',[  6.968750000e-01     1.058333e-01     3.03120e-01     7.458333e-01]);
            
            subplot(2,1,1);
            %imshow(Bgray,'xdata',T,'ydata',F);axis('xy');
            imagesc(T,(F),Bfilt);axis('xy');colormap(flipud(gray));%caxis([-20 50]);
            set(gca,'fontweight','bold','fontsize',14);
            
            title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
            format long e
            disp(cstart);
            %image(T,F,B,colormap('jet'));
            xlabel('Time (sec)');ylabel('Frequency (Hz)');
            subplot(2,1,2);
            imagesc(T,flipud(F),sum(final_image,3));axis('xy');title('Final image');%colormap('jet');
            set(gca,'fontweight','bold','fontsize',14);
            
            hold on;
            xlimm=xlim;
            %subplot(3,1,3)
            %BWskel = bwmorph(BWfinal,'thin','inf');
            %for I=1:length(feature),
            %   plot(feature(I).Tfit,feature(I).Ffit,'ko');hold on; title('Fitted contours');
            
            %end
            hold off;
            %ylim([0 500]);
            xlim(xlimm);
            %imagesc(T,flipud(F),BWskel);axis('xy');title('thin');%colormap('jet');
             keyboard;
            save Figure_demo T F Borg Beq Iopen param Back Iridge Icontour initial_final_image final_image feature
            
        end
    end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Regular subfunctions%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Back_close=extract_background_image(Beq,SE,param,dF,dT,T,F,Idebug)
%%Create background image
%Back=Beq;

%%Option 1, use reconstruction...
%Back(Beq<0)=0;
%Iobr = imreconstruct(Bcontour-(max(max(Back))-param.morph.SNRmin), flipud(Back));
%I2=imregionalmax(Iobr);
%%Option 2, be stupid...
Back=Beq;
Back(Beq>=param.morph.SNRmin)=1;
Back(Beq<param.morph.SNRmin)=0;
Back=logical(Back);

% strel_name='rect';
% %strel_name='ball';
% strel_size=[ceil(param.morph.gap_f/dF) ceil(param.morph.gap_t/dT)];
% SE=strel(strel_name,strel_size);
% %SE=strel(strel_name,strel_size(2),strel_size(1));
%Back_close = flipud(imdilate(Back,SE));
Back_close=close_image(Back, param.morph.background.gap_f,param.morph.background.gap_t, dF, dT);
Back_close=flipud(Back_close);

if Idebug>1
    figure(1)
    set(1,'pos',[60    85   564   987]);
    subplot(4,1,1)
    imagesc(T,F,Beq);caxis([0 max(max(Beq))]);axis('xy');
    set(gca,'fontweight','bold','fontsize',14);
    
    title(sprintf('Original equalized image Beq' ));
    
    % subplot(4,1,2)
    % imagesc(T,F,Bcontour);axis('xy')
    %title(sprintf('Contour traced image, SNRmin %6.2f and dynamic
    %ranges %6.2f and %6.2f',param.morph.SNRmin, param.morph.dynamic_range, ...
    %    param.morph.dynamic_range2));
    subplot(4,1,3)
    imagesc(T,F,Back);axis('xy')
    set(gca,'fontweight','bold','fontsize',14);
    
    title(sprintf('Back: Region above %6.2f dB SNR',param.morph.SNRmin));
    subplot(4,1,4)
    imagesc(T,F,flipud(Back_close));axis('xy');
    set(gca,'fontweight','bold','fontsize',14);
    
    title(sprintf('Back_close  background image'));
    
    
end

end



function [feature,labeled_raw,numObjects]=extract_features_from_binary_image(Iclose,T,F,ctime_new,Idebug)

%%%%function [feature,labeled_raw]=extract_features_from_binary_image(Iclose,T,F,ctime_new,Idebug)
%%%Segment a binary image and use bwlabel to extract features.
%%%  Convert features  to physical units.
%%%  Output variables and example values...
%  feature: structure array with following fields...
%   Area: 61
%   BoundingBox: [4x1 double]
%   Centroid: [2x1 double]
%   Centroid_freq: 1.241034836065574e+02
%   Eccentricity: 9.731259148550847e-01
%   Extent: 3.696969696969697e-01
%   Image: []
%   Orientation: -3.078138279993195e+01
%   Solidity: 7.820512820512820e-01
%   comment: ''
%   ctime: 1.222830698012441e+09
%   duration: 4.800000000000004e-01
%   fmax: 1.445312500000000e+02
%   fmin: 1.054687500000000e+02
%   global_bandwidth: 3.906250000000000e+01
%   local_bandwidth: 1.562500000000000e+01
%   time_band_product: 7.625000000000007e+00
%  labeled_raw:  binary image with blobs indicated by integers
%  numObjects: number of blobs in image

dT=T(2)-T(1);
dF=F(2)-F(1);

if Idebug>0,
    bounds=bwboundaries(Iclose,8,'noholes');
end

%%Gather feature statistics on largest objects
%% First, segment image
[labeled_raw,numObjects]=bwlabel(Iclose,8);


%feature1=regionprops(labeled_raw,'all');
%% Extract key features
feature1=regionprops(labeled_raw,'Area','Centroid','BoundingBox','Orientation','Eccentricity', ...,
    'Extent','Solidity','Image');
%feature1=regionprops(labeled_raw,'Area','Centroid','BoundingBox','Orie
%ntation','Eccentricity','Extent');
%       'Area'              'ConvexHull'    'EulerNumber'
%       'Centroid'          'ConvexImage'   'Extrema'
%       'BoundingBox'       'ConvexArea'    'EquivDiameter'
%       'SubarrayIdx'       'Image'         'Solidity'
%       'MajorAxisLength'   'PixelList'     'Extent'
%       'MinorAxisLength'   'PixelIdxList'  'FilledImage'
%       'Orientation'                       'FilledArea'
%       'Eccentricity'                      'Perimeter'



%%convert units from pixel to physical units...
%% And extract further derived features from data..
%%% Outputs Centroid, BoundingBox, duration, global_bandwith,
%%% time-band-product, fmin, fmax, and local_bandwidth

for I=1:length(feature1),
    
    sf=convert_pixel_units(feature1(I),F,dF,dT);
    
    %% Correct c-time of object...
    sf.ctime=ctime_new+sf.BoundingBox(1);
    sf.comment='';
    feature(I)=sf;
    
end
if isempty(feature1)
    feature=[];
end


end



function feature=extract_features_from_grayscale_image(feature_in,Binput,labeled_raw,numObjects,T,F,Idebug,want_slope)
%%feature=extract_features_from_grayscale_image(feature_in,Binput,labeled_raw,numObjects,T,F,Idebug)
%% Input:
%%  feature_in: original features extracted from binary image
%%  Binput: equalized spectrogram or other intensity image.  Values of
%%      intensity can be less than zero.
%%  labeled_raw:  labeled image (each blob/segment pixel shares same integer)
%%  numObjects: number of segments(blobs) in labeled_raw.
%%  T,F: x and y axes of Bipnut
%   want_slope:  if not empty, then return time and frequency information
%% Features extracted from image:
%%      feature.robust_fmin,robust_fmax,robust_bandwidth,
%%          median_local_kurtosis,fstart,fend
feature=feature_in;

for I=1:numObjects
    
    %mask the equalized spectrogram by segment of interest
    mask=(Binput-min(min(Binput))).*ismember(labeled_raw,I);
    
    if Idebug>2,
        figure;
        subplot(4,1,1);
        imagesc(T,F,flipud(Binput));axis('xy');
        title('Original equalized image');
        subplot(4,1,2);
        imagesc(T,F,flipud(labeled_raw));axis('xy');
        title('Original labeled_raw');
        subplot(4,1,3);
        imagesc(T,F,flipud(ismember(labeled_raw,I)));axis('xy');
        title('Original template');
        subplot(4,1,4);
        imagesc(T,F,flipud(mask));axis('xy');
        title('Segmented selection');
        keyboard
        
    end
    
    try
        field_name={'robust_fmin','robust_fmax','robust_bandwidth','median_local_kurtosis','fstart','fend'};
        for J=1:length(field_name)
            feature(I).(field_name{J})=[];
        end
        
        %T_fit{1}=[];F_fit{1}=[];
        Fup=flipud(F);
        
        
        %%Trim mask to cover only times it exists
        Itrim=find(sum(abs(mask))>0);
        T1=T(Itrim);
        mask=mask(:,Itrim);
        
        if isempty(mask)
            return
        end
        
        mask_freq=mask.*(Fup*ones(1,size(mask,2))); %Weight by frequency
        mean_F=sum(mask_freq)./sum(mask);
        
        
        %Local bandwidth estimate via integral of (f-mean(f))^2
        mask3=mask.*((flipud(F)*ones(1,length(mean_F))-ones(length(F),1)*mean_F).^2);
        local_bandwidth=2*sqrt(sum(mask3)./sum(mask));
        Igood=find(~isnan(local_bandwidth)&local_bandwidth>0);
        feature(I).robust_bandwidth=median(local_bandwidth(Igood));
        
        mask4=mask.*((F*ones(1,length(mean_F))-ones(length(F),1)*mean_F).^4);
        fourth_moment=(sum(mask3)./sum(mask));
        kurtosis=fourth_moment./((local_bandwidth/2).^4);
        feature(I).median_local_kurtosis=median(kurtosis(~isnan(kurtosis)));
        
        
        Igood=find(~isnan(mean_F));
        feature(I).robust_fmin=min(mean_F(Igood));
        feature(I).robust_fmax=max(mean_F(Igood));
        feature(I).fstart=mean_F(Igood(1));
        feature(I).fend=mean_F(Igood(end));
        if feature(I).fmin<0
            keyboard;
        end
        
        %%Measure slope and curvature..
        if exist('want_slope') && ~isempty(want_slope)
            feature(I).Tcall{1}=T1(Igood)';
            feature(I).Fcall{1}=mean_F(Igood)';
            
            
            %             if length(trc.F_call)>4,
            %                 %Fit curve
            %                 [trc.afit,ss,mu] = polyfit(feature(I).Tcall,feature(I).Tcall,3);
            %                 T_fit{1}=T1(T1>min(feature(I).Tcall)&T1<max(feature(I).Tcall));
            %                 F_fit{1} = polyval(trc.afit,T_fit{1},ss,mu);
            %
            %                 fslope=polyval([0 3*trc.afit(1) 2*trc.afit(2)
            %                 trc.afit(3)],T_fit{1},ss,mu);
            %                 fcurve=polyval([0 0 6*trc.afit(1)
            %                 2*trc.afit(2)],T_fit{1},ss,mu);
            %                 Icenter=floor(length(T_fit{1})/2);
            %
            %                % fslope=fslope(Icenter);
            %                 %fcurv=fcurve(Icenter);
            %
            %
            %             else
            %                 fslope=0;
            %                 fcurv=0;
            %                 F_fit{1}=zeros(size(T1));
            %
            %             end
        end
        %
        % %     %Measure slope at center time.
        %      fmin=min(trc.F_call);
        %      fmax=max(trc.F_call);
        %
        %     fstart=trc.F_call(1);
        %     fend=trc.F_call(end);
        %
        %     fpower=Fup(Imax(Iglobal));
        
        
        if Idebug>2
            subplot(2,1,1);
            hold on
            plot(T1,mean_F,'wx',T1,mean_F+local_bandwidth/2,'w',T1,mean_F-local_bandwidth/2,'w')
            
            hold off;
            %     subplot(3,1,3);
            %     imagesc(T1,flipud(F),(mask));axis('xy');hold on
            %     plot(trc.T_call,trc.F_call,'wo',T_fit{1},F_fit{1},'w',T_fit{1},F_fit{1}+local_bandwidth/2,'y',T_fit{1},F_fit{1}-local_bandwidth/2,'y','linewidth',1);hold off
            %     title('Segment subjected to dynamic range trim and curve fit');
            %     keyboard;
            %     close
        end
        
    catch
        keyboard;
    end
    %[sorted_bandwidth,Ipass]=sort(local_bandwidth);
    %Ipass=Ipass(1:round(param.morph.percent_safety_trim*length(Ipass)));
    %local_bandwidth=median(local_bandwidth(Ipass));
    
    feature(I).comment='';
end
end



%%%%%%%%%%%%%%%%%%%%%convert_pixel_units%%%%%%%%%%%%%%%%%%%%%
%%% Convert blob pixel measurements into physical units.
%%% Outputs Centroid, Centrod_freq, BoundingBox, duration, global_bandwith,
%%% time-band-product, fmin, fmax, and local_bandwidth
%%  function fout=convert_pixel_units(feature,F,dF,dT)

function fout=convert_pixel_units(feature,F,dF,dT)

fout=feature;

%%Determine frequency pixels of bounding box
Ifreq(1)=floor(feature.BoundingBox(2));
Ifreq(2)=Ifreq(1)+ceil(feature.BoundingBox(4));

%%Note that image may be cropped, and lowest frequencies are at TOP
%%of image...
fout.Centroid(2)= min(F)+(length(F)-feature.Centroid(2))*dF;  %Centroid frequency
fout.Centroid(1)= feature.Centroid(1)*dT-dT/2;  %Centroid time.
fout.Centroid=fout.Centroid';

fout.Centroid_freq=fout.Centroid(2);
fout.BoundingBox(1)=0.5*dT+dT*(feature.BoundingBox(1)-1);  %Start time of bounding box
fout.BoundingBox(2)= F(length(F)-Ifreq(1));  %Max frequency
fout.BoundingBox(3)= feature.BoundingBox(3)*dT;  %Duration of segment
fout.BoundingBox(4)= (feature.BoundingBox(4)-1)*dF;  %Bandwidth of segment
fout.BoundingBox=fout.BoundingBox';

fout.duration=fout.BoundingBox(3);
fout.global_bandwidth=fout.BoundingBox(4);
fout.time_band_product=feature.Area*dT*dF;

%%Compute estimate of local bandwidth...
fout.fmax=fout.BoundingBox(2);
fout.fmin=fout.fmax-fout.BoundingBox(4);
fout.local_bandwidth=dF*median(sum(feature.Image,1));  %Added the sum(,1) term in case image only one row.


fout.Image=[];

%Compute perimeter/area ratio..
%fout.Perimeter_ratio=feature.Perimeter./feature.Area;


end

%%%%%%%%%%%%%%%%%%close_image.m%%%%%%%%%%%%%%%%%%%%
%% Morphologically close a binary image (build links) between blobs
% function [Iclose,SE,debug_params_close]=close_image(Iopen, gap_f,gap_t, dF, dT)
function [Iclose,SE,debug_params_close]=close_image(Iopen, gap_f,gap_t, dF, dT)
strel_name='rect';
%strel_name='ball';
strel_size=[ceil(gap_f/dF) ceil(gap_t/dT)];
SE=strel(strel_name,strel_size);
%SE=strel(strel_name,strel_size(2),strel_size(1));
Iclose = imclose(Iopen,SE);
%Iclose=Iopen;
debug_params_close.strel_name=strel_name;
debug_params_close.strel_size=strel_size;
debug_params_close.Iclose=Iclose;
end

%%%%%%%%%%[feature,final_image_cell]=link_features(feature,param,input_image,Idebug)
%%% Link objects in image that overlap or merge, and add
%%% 'Totalfmax','Totalfmin','Totalduration', and 'TotalBoundingBox' to
%%% the feature set.


function [feature,final_image_cell]=link_features(feature,param,input_image,Idebug)

%If only one feature, update and retrun...
if length(feature)==1,
    feature(1).Totalfmax=feature(1).fmax;
    feature(1).Totalfmin=feature(1).fmin;
    feature(1).Totalduration=feature(1).duration;
    feature(1).TotalBoundingBox=feature(1).BoundingBox;
    %feature(1).BoundingBox=feature(1).BoundingBox';
    %feature(1).Centroid=feature(1).Centroid';
    
    feature(1).Nharmonic=0;
    feature(1).harmonic_spacing=0;
    feature(1).comment='';
    
    final_image_cell=input_image;
    return
end

for I=1:length(feature)
    feature(I).Nharmonic=0;
    feature(I).harmonic_spacing=0;
end
%Sort so largest time-bandwidth product comes first
[junk,Isort]=sort([feature.fmin],'ascend');
feature=feature(Isort);

%Useful intermediate variables
Bbox=reshape([feature.BoundingBox],4,length(feature))';
tstart=Bbox(:,1);
duration=Bbox(:,3);
tend=duration+tstart;
flow=[feature.fmin];
fhi=[feature.fmax];
fstart=[feature.fstart];
fend=[feature.fend];

names=fieldnames(feature);
Ibad=[];
Nimages=0;

input_image_org=input_image;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Splice components together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I=1:length(tstart)
    if feature(I).ctime<0
        continue
    end
    
   
        Idebug=0;
    
    if Idebug>3
        final_image0=ismember(input_image,Isort(I));
    
        figure(4);clf
        subplot(3,1,1)
        imagesc(input_image_org)
        title('original input image')
        subplot(3,1,2)
        imagesc(final_image0);
        title(sprintf('I=%i',I));
        subplot(3,1,3);
        imagesc(input_image);
        title(sprintf('current image before splicing: I=%i',I));
        
       
        
    end
    
    
    for J= (I+1):length(tstart)
        
        merge_flag=0;
        if feature(J).ctime<0  %If this has already been merged
            continue
        end
       
        test_spliceb=(tstart(I)-tend(J))<=param.merge.gap_t&(tstart(I)-tend(J))>0;
        if test_spliceb
            [junk,Imin1]= min(abs(fstart(I)-fend(J)));
            test_spliceb= (junk<=param.merge.gap_f);
        end
        
        test_splicea=(tstart(J)-tend(I))<=param.merge.gap_t&(tstart(J)-tend(I))>0;
        if test_splicea
            [junk,Imin1]= min(abs(fstart(J)-fend(I)));
            test_splicea= (junk<=param.merge.gap_f);
        end
         temp_img=ismember(input_image,Isort(J));
           
        if test_spliceb||test_splicea
            
            %Idebug_global=1;
            merge_flag=1;
            feature(I)=merge_features(feature(I),feature(J),Imin1);
            Ibad=[Ibad J];  %Flag indicies to remove from feature object...
            %Update image of call
            %final_image=final_image+ismember(input_image,Isort(J));  %Final image should have same indicies for pieces.
            
            %%Input image should be altered to eliminate the distinction for later harmonic merging...
            Ichange=find(input_image==Isort(J));
            input_image(Ichange)=Isort(I);
            %feature(I).SEL=sum(10.^(feature(I).SEL/10)+10.^(feature(J).SEL/10));
        %else
           % temp_img=[];
            
        end
        
        if merge_flag>0  %Eliminate the higher frequecy segment
            feature(J).ctime=-1;  %Label the 'J' feature to be destroyed
            
            %Finally, update bounding box of 'I' so that future
            %segments might overlap combined call, even if it
            %doesn't overlap the original 'I' segment
            [tstart(I),Imin]=min([tstart(I) tstart(J)]);
            ftemp=[fstart(I) fstart(J)];
            fstart(I)=ftemp(Imin);
            [tend(I),Imax]=max([tend(I) tend(J)]);
            ftemp=[fend(I) fend(J)];
            fend(I)=ftemp(Imax);
            
            flow(I)=min([flow(I) flow(J)]);
            fhi(I)=max([fhi(I) fhi(J)]);
            
            
        end
        %         if I==2
        if Idebug>3
            figure(5);
            subplot(2,1,1)
            imagesc(temp_img);
            title(sprintf('J=%i',J));
            subplot(2,1,2)
            imagesc(input_image);
            title(sprintf('splice loop: I=%i, J=%i merge_flag=%i',I,J,merge_flag));
            pause;
           
            close(5);
        end
    end %J
    
  
    
end  %I


%%Link components together
for I=1:(length(tstart)),
    if feature(I).ctime<0
        continue
    end
    
    final_image0=ismember(input_image,Isort(I));
    final_image=final_image0;
    Nharmonic=0;
    harmonic_spacing=0;
    if Idebug>3
        figure
        subplot(2,1,1)
        imagesc(input_image)
        subplot(2,1,2);
        imagesc(final_image);
        title(sprintf('link: I=%i',I));
        pause;
        %close
        
    end
    
    
    for J= (I+1):length(tstart)
        
        merge_flag=0;
        if feature(J).ctime<0  %If this has already been merged
            continue
        end
        % disp(J);
        %Overlap test--for feature with harmonics use harmonic with
        %    maximum frequency.
        %    Note that feature J will never have harmonics or split
        %    ends,and thus feature(J).(field) will always be a
        %    single element...
        tovlap=min([tend(I) tend(J)])-max([tstart(I) tstart(J)]);
        ovlap=tovlap./min(duration([I J]));  %Negative if no overlap
        
        fovlap=min([fhi(I) fhi(J)])-max([flow(I) flow(J)]);
        fovlap=fovlap./min([fhi(J)-flow(J) fhi(I)-flow(I)]);
        %two splice tests...revision on 4/16/2009 allows splice
        %with any harmonic of feature(I), not just maximum
        %   frequency
        test_ovlap1=(ovlap>=param.merge.ovlap & abs(fhi(I)-flow(J))<param.merge.max_frequency_separation);
        test_ovlap2=(ovlap>=param.merge.ovlap & fovlap>param.merge.ovlap);
        %test_splice1=abs(tstart(I)-tend(J))<=param.merge.gap_t & abs(max(feature(I).fstart)-(feature(J).fend))<=param.merge.gap_f;
        %test_splice2=abs(tstart(J)-tend(I))<=param.merge.gap_t & abs((feature(J).fstart)-max(feature(I).fend))<=param.merge.gap_f;
        
      
        
        if test_ovlap1||test_ovlap2  %harmonics
            Nharmonic=Nharmonic+1;
            harmonic_spacing=harmonic_spacing+abs(feature(I).Centroid(2)-feature(J).Centroid(2));
            merge_flag=1;
            %If overlap, merge into lower frequency bound
            %%Copy features into 'I' feature
            %names=fieldnames(feature);
            feature(I)=clump_features(feature(I),feature(J),names,'harmonic');
            
            Ibad=[Ibad J];
            %Update image of call
            final_image=final_image+(1+length(Ibad))*ismember(input_image,Isort(J));
            %% %Check if segments can be spliced together..
       
        end
        
        if merge_flag>0  %Eliminate the higher frequecy segment
            feature(J).ctime=-1;  %Label the 'J' feature to be destroyed
            
            %Finally, update bounding box of 'I' so that future
            %segments might overlap combined call, even if it
            %doesn't overlap the original 'I' segment
            tstart(I)=min([tstart(I) tstart(J)]);
            tend(I)=max([tend(I) tend(J)]);
            flow(I)=min([flow(I) flow(J)]);
            fhi(I)=max([fhi(I) fhi(J)]);
            
        end
        %         if I==2
        if Idebug>3
            figure(5);
            imagesc(final_image);
            title(sprintf(' link loop: I=%i, J=%i',I,J));
            pause;
            close(5);
        end
    end %J
    
    Nimages=Nimages+1;
    final_image_cell(:,:,Nimages)=final_image;
    %Idebug=4;
    if Idebug>3
        figure
        subplot(3,1,1)
        imagesc(final_image0);
        title('tested segment')
        subplot(3,1,2)
        imagesc(input_image)
        title('original spliced input image')
        subplot(3,1,3);
        imagesc(final_image);
        title('current modified image');
        pause;
        close
        
    end
    
    if Nharmonic>0
        feature(I).Nharmonic=Nharmonic;
        feature(I).harmonic_spacing=harmonic_spacing/Nharmonic;
    end
    
end  %I



%%Elminate features logged in 'Ibad'
Igood=setdiff(1:length(tstart),Ibad);
feature=feature(Igood);

%Add a field for total Bounding box of feature
for I=1:length(feature),
    [feature(I).TotalBoundingBox,feature(I).Totalfmax,feature(I).Totalfmin,feature(I).Totalduration]=merge_BoundingBox(feature(I).BoundingBox);
end
%%Fix harmonic feature parameters

end

function feature_out=merge_features(feature1,feature2,Iharm)
%Merge two blobs into one.  Only feature1 can have harmonics, and if
%Iharm exists, merge feature2 into Iharm'th component of feature1

if ~exist('Iharm')
    Iharm=1;
end
feature_out=feature1;
feature_out.BoundingBox(:,Iharm)=merge_BoundingBox([feature1.BoundingBox(:,Iharm) feature2.BoundingBox]);
feature_out.local_bandwidth(Iharm)=median([feature1.local_bandwidth(Iharm) feature2.local_bandwidth]);
feature_out.global_bandwidth(Iharm)=feature_out.BoundingBox(4,Iharm);

feature_out.fmin(Iharm)=min([feature1.fmin(Iharm) feature2.fmin]);
feature_out.fmax(Iharm)=max([feature1.fmax(Iharm) feature2.fmax]);
feature_out.duration(Iharm)=feature_out.BoundingBox(3,Iharm);
feature_out.Area(Iharm)=sum([feature1.Area(Iharm) feature2.Area]);
feature_out.Centroid(:,Iharm)=[feature_out.BoundingBox(1,Iharm)+0.5*feature_out.BoundingBox(3,Iharm); 0.5*feature_out.fmin(Iharm)+0.5*feature_out.fmax(Iharm)];
feature_out.Centroid_freq(Iharm)=feature_out.Centroid(2);
feature_out.time_band_product(Iharm)=feature1.time_band_product(Iharm)+feature2.time_band_product;

%feature_out.fpower=max([feature1.fpower feature2.fpower]);
feature_out.ctime(Iharm)=min([feature1.ctime(Iharm) feature2.ctime]);
feature_out.comment=[feature1.comment sprintf(',splice%i',Iharm)];

if isfield(feature_out,'robust_fmin')
    feature_out.robust_fmin(Iharm)=min([feature1.robust_fmin(Iharm) feature2.robust_fmin]);
    feature_out.robust_fmax(Iharm)=max([feature1.robust_fmax(Iharm) feature2.robust_fmax]);
    feature_out.robust_bandwidth(Iharm)=median([feature1.robust_bandwidth(Iharm) feature2.robust_bandwidth]);
    feature_out.median_local_kurtosis(Iharm)=median([feature1.median_local_kurtosis(Iharm) feature2.median_local_kurtosis]);
    if feature_out.ctime(Iharm)<=feature2.ctime,
        feature_out.fend(Iharm)=feature2.fend;
    else
        feature_out.fstart(Iharm)=feature2.fstart;
    end
    if isinf(feature_out.robust_fmax(Iharm))
        keyboard;
    end
end
if isfield(feature_out,'Fcall')
    tabs=[feature1.Tcall{1}; feature2.Tcall{1}];
    ff=[feature1.Fcall{1}; feature2.Fcall{1}];
    [tabs,Isort]=sort(tabs);
    feature_out.Tcall{1}=tabs;
    feature_out.Fcall{1}=ff(Isort);
    
end

end

%%%%%%%clump features... assign to common index but keep individual
%%%%%%%features--used for harmonics and features that share a common
%%%%%%%background contour...
function feature_out=clump_features(feature1,feature2,names,comment)
for K=1:length(names)
    
    value=feature2.(names{K});
    %if size(value,2)==1
    feature_out.(names{K})=[feature1.(names{K}) value];
    
    if strcmp(names{K},'Fcall')
       % feature_out.Fcall{1}=feature1.Fcall{1};
       % feature_out.Fcall{2}=feature2.Fcall{1};
        feature_out.Fcall={feature1.Fcall{:} feature2.Fcall{1}};
    end
    if strcmp(names{K},'Tcall')
        %feature_out.Tcall{1}=feature1.Tcall{1};
        %feature_out.Tcall{2}=feature2.Tcall{1};
        feature_out.Tcall={feature1.Tcall{:} feature2.Tcall{1}};
    end
end
feature_out.comment=[feature1.comment ',' comment ',' feature2.comment];
end

%%merge Bounding Boxes
%% Bounding Box is a matrix with [4 Nblob] dimension.
%% Reduce to a single [4 1] vector...
%% Elements are [tstart maxfreq duration bandwidth]'
function [TotalBoundingBox,Totalfmax,Totalfmin,Totalduration]=merge_BoundingBox(BoundingBox)

Totalfmax=max(BoundingBox(2,:));
Totalfmin=min(BoundingBox(2,:)-BoundingBox(4,:));
tstart=(BoundingBox(1,:));
tend=max(tstart+BoundingBox(3,:));
Totalduration=tend-min(tstart);
TotalBoundingBox(1:3,1)=[min(BoundingBox(1,:)) Totalfmax Totalduration]';
TotalBoundingBox(4,1)=Totalfmax-Totalfmin;



end


%function Ipass=crude_filter_feature(feature,param,Idebug)
% Remove features based on certain feature criteria...
function Ipass=crude_filter_feature(feature,param,Idebug)

filter_names={'duration','Centroid_freq','time_band_product','robust_fmin','robust_bandwidth','robust_fmax','Orientation','Eccentricity'};

if Idebug>1,
    disp(sprintf(' %i candidates to test',length(feature)));
end
Ipass=1:length(feature);
for I=1:length(filter_names),
    if ~isfield(feature(Ipass(1)),filter_names{I})
        continue
    end
    try
        Ipass1=find([feature(Ipass).(filter_names{I})]<=param.morph.(filter_names{I}).max & [feature(Ipass).(filter_names{I})]>=param.morph.(filter_names{I}).min);
    catch
        keyboard
    end
    if Idebug>1&findstr(filter_names{I},'robust_bandwidth');
        disp(filter_names{I});
        feature.(filter_names{I})
        param.morph.(filter_names{I})
        
    end
    if isempty(Ipass1),
        if Idebug>1
            disp(sprintf('morph: No candidates pass %s test,: min: %10.6f max: %10.6f', ...
                filter_names{I},param.morph.(filter_names{I}).min,param.morph.(filter_names{I}).max));
            disp(sprintf('Actual values for features: %10.6f \n',feature(Ipass).(filter_names{I})));
            %figure;imagesc(labeled_raw);pause;close;
            
        end
        Ipass=[];
        
        return;
        
    end
    if Idebug>1,
        disp(sprintf('%i candidates pass %s test',length(Ipass1),filter_names{I}));
    end
    Ipass=Ipass(Ipass1);
    
end

end

%%%function Ibin=threshold_image_by_reconstruction(Beq,param,dT,dF);
%% Convert equalized spectrogram into a binary image, reagardless of SNR
%% level.
function [Ibin,Beq,debug_rec]=threshold_image_by_reconstruction(Beq,param,dT,dF,Idebug)
debug_rec=[];
conn=8;
Beq(Beq<(param.morph.SNRmin-param.morph.dynamic_range))=0;
se_area=(ceil(param.morph.time_band_product.min/(dT*dF)))/pi;
debug_rec.se_radius=ceil(sqrt(se_area));
%debug_rec.se_radius=1;
se = strel('diamond',debug_rec.se_radius);

%Remove objects too small in image and flatten peaks so a 2*disksize+1 plateau exists.
% Opening by reconstruction
Ie = imerode(Beq, se);
Iobr = imreconstruct(Ie, Beq);

Ih= imhmax(Iobr,param.morph.dynamic_range,conn);  %Remove top dynamic_range dB from peaks--thus only peaks this tall survive
Ip=(Iobr-Ih);  %Subtract plateaus from peaks to create peaks that are at most dynamic_range high--this step
%permits weak and high SNR signals to be extracted
%simultaneously, as long as peaks in question have
%dynamic_range dB

%Reject 'islands' with max SNR less than SNRmin_peak

% Create a second 'h-dome' that retains only peaks that are SNRmin_peak
%  high--note this is processing original image instead of reconstructed
%  image, beacuse the process of opening reduces SNR size.

Ih2=imhmax(Beq,param.morph.SNRmin,conn);
Ip2=Ip.*sign(Ih2);

%Note: can't do this because I need to preserve data points between SNRmin
% abd SNRmin-dynamic_range
%Ih2=(Beq<param.morph.SNRmin);
%Ip2=Ip;
%Ip2(Ih2)=0;

%Now convert to binary.
Ibin=imregionalmax(sign(Ip2));


if Idebug>=2
    debug_rec.Beq=Beq;
    debug_rec.Iobr=Iobr;
    debug_rec.Ih=Ih;
    debug_rec.Ip=Ip;
    debug_rec.Ih2=Ih2;
end

if 1==0
    figure
    % Icol=10;
    plot(Beq(:,Icol),'ko-');hold on
    plot(Iobr(:,Icol),'g');
    plot(Ih(:,Icol),'k--');
    plot(Ip2(:,Icol),'r');
    legend('original','opened','h-dome','final');
    
end

end

%%function
%%[adjust_threshold,Ibroad_old,Bgray]=adjust_threshold_level(adjust_threshold,feature,param,Ibroad_old,labeled_raw,Beq,Bgray,Idebug)
%%  Replaced by threshold_image_by_reconstruction.m
%% given a set of features that exceed solidity and bandwidth restrictions,
%% attempt to readjust image to reduce "saturation" of threshold.
%% Goal is to break apart harmonics and improve high-intensity call
%% conversion

