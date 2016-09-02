%function [stats,BWfinal,Bmean]=extract_image_features_simple(y,cstart,param,Idebug,Bmean)
%  Simple morphological processing demo
%%%%%%%%%%%%%extract_image_features_simple.m%%%%;
%% take raw input data and y (start time ctime in c-time) and attempt to find a contour in spectrogram
%% using parameter fields in param:
%%    param.Fs=1000;
%     param.Nfft=128;
%     param.ovlap=7/8;
%
%     param.morph.eq_time=0.75;  %Time to use to create background ambient noise estimate.
%     param.morph.threshold_fudge=what fraction of automated adaptive threshold should be
%           used for converting to binary image...
%           of the background noise
%     param.morph.SNRmin=10;  %min dB level, for thresholding background noise
%
%     param.median.size=[0.2 50] (sec, Hz); %Used to derive bandwidth of
%       median filter 
%     param.filter.size=[0.2 50] (sec,Hz); %LoG filter size
%     param.filter.sigma=0.75; %standard deviation for LoG parameter
%
%  Clear small objects
%     param.morph.TimeBandProduct=0.2*50; %Sets minimum spectrogram pixel
%           area to count as an 'object'; used by bwareaopen
%
%  Dilate to close gaps
%     param.morph.gap_f=4*param.Fs/param.Nfft;
%     param.morph.gap_t=0.1;  %maximum time to permit a gap in fm signal in sec

%  Filter out undesirable features, like too short a time or too large a
%  bandwidth...
%     param.morph.minDuration=0.3;  %minimum length in seconds a contour will have to be
%     param.morph.maxBW=20;
%
%     cstart is a c-time of the start of the input data y
%     Beq is an inital equalization spectrum.  When first called the
%       routine will output a value that can be reused.
%     call_nopulses: removes pulses, but does not apply contour tracer...
%function [stats,BWfinal,Bmean]=extract_image_features_simple(y,cstart,param,Idebug,Bmean)

% Apr 21--added an SEL and SNR feature...

function [stats,BWfinal,Bmean]=extract_image_features_simple(y,cstart,param,Idebug,Bmean)
threshold_fudge=param.morph.threshold_fudge;

stats=[];BWfinal=[];Bmean=[];

if isempty(y),
    return
end

if nargin<4,
    Idebug=0;
    Bmean=[];
elseif nargin<5,
    Bmean=[];

end

mydir=pwd;

Fs=param.Fs;

[S,F,T,PSD]=spectrogram(y,param.Nfft,round(param.ovlap*param.Nfft),param.Nfft,Fs,'yaxis');
dT=T(2)-T(1);dF=F(2)-F(1);

B=10*log10(abs(PSD+eps));  %Power spectral density..


%Equalize using median value of data before param.eq_time
%Beq=equalize_background(B,T,param);
Imedian_sample=find(T<param.morph.eq_time);
if length(Imedian_sample)>1,
    Bmean=median(B(:,Imedian_sample).').';  %Note that the median is unaffected by a log transformation of the data
    %noise_PSD=dF*mean(sum(PSD(:,Imedian_sample)).').';

else
    Bmean=B(:,Imedian_sample);
    %noise_PSD=PSD(:,Imedian_sample);

end
Beq=B-Bmean*ones(1,size(B,2));
%Imax=max(Imedian_sample);
%Beq=Beq(:,Imax:end);
%B=B(:,Imax:end);
%PSD=PSD(:,Imax:end);

%T=T(max(Imedian_sample):end);

%Convert to gray level-plots, using SNRmin as the cutoff.
Bgray=mat2gray(Beq,[param.morph.SNRmin 40]);
Bgray=flipud(Bgray);

%Optional median filter
if param.median.on==1,
    min_duration=ceil(param.median.size(1)/(dT));
    min_freq=ceil(param.median.size(2)/(dF));

    Imed=medfilt2(Bgray,[min_freq min_duration],'symmetric');
else
    Imed=Bgray;

end
%%Optional Gaussian filtering
if param.filter.on==1,    %apply gaussian filter to B
    %disp('applying gaussian filter');

    min_duration=ceil(param.filter.size(1)/dT);
    min_freq=ceil(param.filter.size(2)/dF);

    Hgauss = fspecial('gaussian',[min_freq min_duration],param.filter.sigma);
    Ifilt=imfilter(Imed,Hgauss,'symmetric');
    if Idebug>0,
        figure(3);
        subplot(3,1,1);
        imagesc(T,F,Imed);
        title('Original spectrogram');
        %caxis([param.SNRmin 40]); %colorbar;
        axis('xy');

        subplot(3,1,2);
        imagesc(T,F,Ifilt);
        title(sprintf('Gaussian filter with  size %s, sigma %6.2f',mat2str(param.filter.size),param.filter.sigma));
        %caxis([param.SNRmin 40]); %colorbar;
        axis('xy');

    end
    Imed=Ifilt;
end

%Convert to binary
level=graythresh(Imed);
Ibin=im2bw(Imed,threshold_fudge*level);

%Open image to clear out small objects
%Iopen=imopen(Ibin,strel(strel_name,strel_size));

Iopen=bwareaopen(Ibin,ceil(param.morph.TimeBandProduct/(dT*dF)));

%%Dilate image to fill in gaps
strel_name='rect';
strel_size=[ceil(0.5*param.morph.gap_f/dF) ceil(0.5*param.morph.gap_t/dT)];
SE=strel(strel_name,strel_size);
%Dilate in frequency and time
%Iopen=imdilate(Iopen,strel('line',strel_size(2),0));
%Iopen=imdilate(Iopen,strel('line',strel_size(1),90));
Imdilate=imdilate(Iopen,SE);

%%Erode to remove effects of dilation
Ierode = imerode(Imdilate,SE);

if Idebug>0,
    bounds=bwboundaries(Imdilate,8,'noholes');
end

%%Gather feature statistics on largest objects
[labeled,numObjects]=bwlabel(Ierode,8);
%stats1=regionprops(labeled,'all');
stats1=regionprops(labeled,'Area','Centroid','BoundingBox','Orientation','Eccentricity','Image');

%       'Area'              'ConvexHull'    'EulerNumber'
%       'Centroid'          'ConvexImage'   'Extrema'
%       'BoundingBox'       'ConvexArea'    'EquivDiameter'
%       'SubarrayIdx'       'Image'         'Solidity'
%       'MajorAxisLength'   'PixelList'     'Extent'
%       'MinorAxisLength'   'PixelIdxList'  'FilledImage'
%       'Orientation'                       'FilledArea'
%       'Eccentricity'                      'Perimeter'
if isempty(stats1),
    return
end
%Find all areas that meet minimum criteria
duration=[stats1.BoundingBox];
duration=dT*duration(3:4:end);  %Duration of each component
minArea=ceil(param.morph.TimeBandProduct/(dT*dF));

%bandwidth=[stats1.BoundingBox];
%bandwidth=dF*bandwidth(4:4:end);
bandwidth=[];
for IJ=1:length(stats1),
    bandwidth(IJ)=dF*max(sum(stats1(IJ).Image));
end

idx=find(duration>=param.morph.minDuration &bandwidth<=param.morph.maxBW&minArea<=[stats1.Area]);

%Do I reject image if large airgun signal is present?

if isempty(idx),
    return;
end
%idx=find([stats1.Area]>min_freq*min_duration);
[junk,Isort]=sort([stats1(idx).Area],'descend');
Iobject=idx(Isort);
stats=stats1(Iobject);

%Final image contains only largest objects
BWfinal=ismember(labeled,idx);
labeled_tmp=flipud(labeled);


if length(stats)>0,
    [value,Imax]=max(sum(BWfinal));
    stats(1).ctime=cstart+(dT*Imax-(param.energy.bufferTime-param.morph.eq_time));

    stats(1).vertical_percent=value/size(BWfinal,1);

    stats(1).dF=dF;
    stats(1).dT=dT;

    if length(Isort)>1
        stats(1).Area2=sum([stats(1:2).Area]);
    else
        stats(1).Area2=stats(1).Area;
    end
    %

    %%Unit conversions from pixel to physical units...
    Beq2=Beq;
    for I=1:length(stats),
        stats(I).Centroid(2)= (length(F)-stats(I).Centroid(2))*dF;  %%Center frequency..
        stats(I).Centroid(1)= stats(I).Centroid(1)*dT-param.energy.bufferTime;
        stats(I).BoundingBox(3)= stats(I).BoundingBox(3)*dT;  %Duration of segment
        stats(I).BoundingBox(4)= stats(I).BoundingBox(4)*dF;  %Bandwidth of segment
        stats(I).LocalBandwidth=dF*max(sum(stats(I).Image));
        MaxFrequency=(length(F)-stats(I).BoundingBox(2))*dF;

        Ifreq=MaxFrequency-[stats(I).BoundingBox(4) 0];
        Ifreq=[max([1 floor(Ifreq(1)/dF)]) ceil(Ifreq(2)/dF)];
        Ifreq=Ifreq(1):Ifreq(2);
        %keyboard;
        noise_rms=mean(PSD(Ifreq,Imedian_sample).');
        noise_rms=trapz(noise_rms)*dF;
        noise_rmsdB=10*log10(noise_rms);
        %%Produce SEL and SNR estimates
        Inum=Iobject(I);
        Ion=find(labeled_tmp(:)==Inum);
        stats(I).SEL=10*log10((1-param.ovlap)*sum(PSD(Ion)));
        Npixel=round(stats(I).BoundingBox(3)/dT);
        signal_rms=sum(PSD(Ion))*dF/Npixel;
        signal_rmsdB=10*log10(signal_rms);
        stats(I).SNR=signal_rmsdB-noise_rmsdB;

        if Idebug>0,
            figure(10);
            Beq2(Ion)=-20;
            subplot(2,1,1);imagesc(T,F,Beq);caxis([-20 30]);colorbar;axis('xy');
            subplot(2,1,2);imagesc(T,F,Beq2);caxis([-20 30]);colorbar;axis('xy');
            title(sprintf(' SNR: %6.2f, SEL: %6.2f dB',stats(I).SNR,stats(I).SEL));
            %keyboard
            pause;
        end
    end


end
%%What type of signal is it?  Downsweep?  Upsweep?

%%Code to pack BWskel into two vectors for classification...
if Idebug>0,

    figure(1);
    set(gcf,'pos',[680   128   603   906]);
    %subplot(3,2,1);
    %imagesc(T,F,Bgray);axis('ij');title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
    colormap(flipud(gray));
    subplot(3,2,1);
    imagesc(T,flipud(F),Imed);axis('xy');
    try,
    title(sprintf('Median filter size %i freq by %i time pixels',min_freq,min_duration));
    end
    subplot(3,2,2);
    imagesc(T,flipud(F),Ibin);axis('xy');
    title('binary conversion using Otsus method');
    subplot(3,2,3);
    imagesc(T,flipud(F),Iopen);axis('xy');
    title(sprintf('Morphological opening minimum area %6.2f Hz-s',param.morph.TimeBandProduct));

    subplot(3,2,4);
    imagesc(T,flipud(F),Imdilate);axis('xy');title(sprintf('imsopen, struct element: %s %i freq pixel by %i time pixel', ...
        strel_name,dF*strel_size(1),dT*strel_size(2)));

     subplot(3,2,5);
    imagesc(T,flipud(F),Ierode);axis('xy');title(sprintf('ierode, struct element: %s %i freq pixel by %i time pixel', ...
        strel_name,dF*strel_size(1),dT*strel_size(2)));

    hold on;
    colors=['b' 'g' 'r' 'c' 'm' 'y'];
    for k=1:length(bounds),
        boundary=bounds{k};
        cidx=mod(k,length(colors))+1;
        plot(boundary(:,2),boundary(:,1),...
            colors(cidx),'Linewidth',2);
    end
    hold off

    subplot(3,2,6);
    %imagesc(T,F,BWfinal);
    imagesc(T,flipud(F),BWfinal);axis('xy');hold off;
    
    
%     try,
% 
%         hold on
%         plot(stats(1).ConvexHull(:,1),stats(1).ConvexHull(:,2),'b');
%         hold off
%         title(sprintf('final processed image, Bounding Box of max bandwidth: %6.2f Hz Log e: %6.5f, Orientation %5.2f Centroid Freq %6.2f',stats(1).BoundingBox(4),-log10(stats(1).Eccentricity), ...
%             stats(1).Orientation,stats(1).Centroid(2)));
%         stats(1)
% 
%     end

    figure(2);
    subplot(3,1,1);
    imagesc(B,'xdata',T,'ydata',F);axis('xy');colormap(flipud(gray));caxis([-20 50]);
    title(sprintf('original grayscale image: %s',datestr(datenum(1970,1,1,0,0,cstart))));
    format long e
    disp(cstart);
    %image(T,F,B,colormap('jet'));
    xlabel('Time (sec)');ylabel('Frequency (Hz)');
    subplot(3,1,2);
    imagesc(T,flipud(F),BWfinal);axis('xy');title('Final image');%colormap('jet');
    subplot(3,1,3)
    BWskel = bwmorph(BWfinal,'thin','inf');
    imagesc(T,flipud(F),BWskel);axis('xy');title('thin');%colormap('jet');

    keyboard;
    figure(1);clf;
    
end

