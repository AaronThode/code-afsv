%%%%extract_ridge_traces.m%%%%%%%%%%
% function [Iopen, Back, detect]=select_peaks(Beq,SNRmin,minimum_bandwidth,maximum_bandwidth,dynamic_range_v,dynamic_range_h,min_area_pixels,Nfft,Fs,F,Idebug)
%  Given an image, output an image of ridges
%  Input:
%       Beq: equalized and potentially filtered spectrogram image in SNR dB units
%       SNRmin:  minimum dB SNR required for an image pixel to be considered
%       minimum_bandwidth, maximum_bandwidth: minimum and maximum bandwidth
%           of ridge to be considered in Hz.
%       dynamic_range_v: sets cutoff point for ridge in vertical direction;
%           ridge edge lies (dynamic_range_v) dB below local ridge maximum
%       dynamic_range_h: ensures that all points in ridge lie within
%           dynamic_range_h dB of local maximum.  Currently deactivitated
%       min_area_pixels:  Minimum contour area in pixels
%       Nfft,Fs:  FFt size used to generate Beq, and sampling rate in Hz.
%       F: vector of frequencies (Hz) associated with Beq
%       Idebug: a scalar, if greater than 2 plot debug output
%   Output:
%       detect: matrix size of Beq, with only ridges remaining (SNR levels
%           preserved)
%       Iopen:  Binary version of detect, followed by morphological opening
%           to remove small objects.  Logical MATLAB format
%       Back: binary image of Beq, where Beq>SNRmin=1, otherwise 0.
%
function [Iopen,Back,detect]=extract_ridge_traces(Beq,SNRmin,minimum_bandwidth,maximum_bandwidth,dynamic_range_v,dynamic_range_h, ...
    min_area_pixels,Nfft,Fs,F,Idebug)

if ~exist('Idebug')
    Idebug=0;
    F=0;
end

%%Step 1--locate local minima along frequency axis with SNRmin dB minimum
%%height
detect=zeros(size(Beq));
Nt=size(Beq,2);
Iminimum_bandwidth=round(0.5*minimum_bandwidth*(Nfft/Fs));
Imax_bandwidth=round(0.5*maximum_bandwidth*(Nfft/Fs));

%%Step 2--define the two thresholds, a minimum SNR, plus a contour trace
Back=Beq;
Back(Beq>=SNRmin)=1;
Back(Beq<SNRmin)=0;

Back2=Beq;
Back2(Beq>=SNRmin+dynamic_range_v)=1;
Back2(Beq<SNRmin+dynamic_range_v)=0;



Icol=(1:size(Beq,1))';  %indicies of rows
for I=1:Nt  %for every image column...
    
    %%%Find peaks from FM sweep
    %%  peaks must exceed minSNR and be far enough from top and
    %% bottom edges.
    %It=find(Beq(:,I)>SNRmin);%Must exceed cutoff
    It=Icol.*Back2(:,I);  %examine only pixels with a minimum (SNRmin+contour value)
    
    %%%Peak frequency must be greater than minimum bandwidth; less than
    %%  max frequency minus minimum bandwidth
    It01=find(It<length(F)-Iminimum_bandwidth-1&It>1+Iminimum_bandwidth);
    It=It(It01);
    
    %%Peak frequency must, indeed, be a local maximum
    It2=find(Beq(It+1,I)<Beq(It,I)&Beq(It-1,I)<Beq(It,I)); %must be a peak
    Ipeaks=It(It2);
    
    SNRpeaks=Beq(Ipeaks,I);
    
    %%For each local maxima at a given time, find upper and lower boundary.
    %% Note that the edges of the contour can be greater than
    %% dynamic_range_v below the peak--creating a bandwidth fatter than that produced by
    %   Back;  thus a "robust bandwidth" is a
    %% better measure...
    myflag=0;
    for J=1:length(SNRpeaks)
        index=Ipeaks(J)+(0:Imax_bandwidth);
        index=index(index<=length(F));
        %%Find upper bound on contour
        Iupper_bound=min(find(Beq(Ipeaks(J),I)-Beq(index,I)>dynamic_range_v));
        I1=index(1:Iupper_bound);
        %%Confirm that no local maxima exists between upper bound and
        %%current local maxima
        monotonic1=~any(Beq(I1,I)>Beq(Ipeaks(J),I));
        
        index2=sort(Ipeaks(J)-(0:Imax_bandwidth));
        index2=index2(index2>0);
        Ilower_bound=max(find(Beq(Ipeaks(J),I)-Beq(index2,I)>dynamic_range_v));
        I2=index2(Ilower_bound:end);
        monotonic2=~any(Beq(I2,I)>Beq(Ipeaks(J),I));
        myflag=0;
        if monotonic1&&monotonic2&&~isempty(Iupper_bound)&&~isempty(Ilower_bound);
            detect(I1,I)=Beq(I1,I);
            detect(I2,I)=Beq(I2,I);
            myflag=1;
        end
        
    end
    
    %%%Optional debug plot of peak-picking
    if Idebug>2&&I>1
        figure(3);
        plot(F/1000,Beq(:,I));ylim([-10 10]);grid on;
        title(sprintf('Time index: %i', I));
        ylabel('Equalized SNR (dB)');xlabel('Freq (kHz)');
        hold on; plot(F(Ipeaks)/1000,Beq(Ipeaks,I),'gx');
        try
            plot(F(I1)/1000,Beq(I1,I),'r');
            plot(F(I2)/1000,Beq(I2,I),'g');
        end
        title(sprintf('green x: local maxima; red: upper frequency ridge; green: lower frequency ridge, detected?: %i',myflag));
        pause;hold off;
    end
end


if Idebug>1
    figure(4)
    set(get(gcf,'Number'),'pos',[118   103   554   980]);
    subplot(2,2,1)
    imagesc([],F,Beq);title('extract_ridge_traces: Beq');axis('xy')
    caxis([0 70]); colorbar('southoutside');
    subplot(2,2,2)
    imagesc([],F,(detect));title('after dynamic_range_v');caxis([0 70]);axis('xy')
    colorbar('southoutside');
    
end



%%Convert 'detect' to binary and remove small objects.
Ibin=sign(abs(detect));
Iopen=bwareaopen(Ibin,min_area_pixels);
%Back=bwareaopen(Back,min_area_pixels);
Iopen=flipud(Iopen);
detect=flipud(detect);
if Idebug>1
    figure(4)
    subplot(2,2,3)
    
    imagesc([],F/1000,(Iopen));title('after opening');
end

%%Step 3... restrict SNR values further when compared on a global scale;
%% trim each blob in image so that it only has points within
%% dynamic_range_h of the maximuum value
if 1==0  %Set 1==0 to retain original function
    [labeled_raw,numObjects]=bwlabel(Iopen,8);
    for I=1:numObjects
        index=find(labeled_raw==I);
        %     figure(5);subplot(3,1,1);
        %     imagesc(labeled_raw==I);
        %     subplot(3,1,2);
        %     detect0=detect;
        %     detect0(labeled_raw~=I)=0;
        %     imagesc(detect0);colorbar('southoutside')
        maxval=max(detect(index));
        Ikill= detect(index)<=maxval-dynamic_range_h;
        Iopen(index(Ikill))=0;
        
        %     detect0(index(Ikill))=0;
        %     subplot(3,1,3);imagesc(detect0);colorbar('southoutside')
        %     pause;
        
    end
end


if Idebug>1
    subplot(2,2,4);
    if isempty(dynamic_range_h)
        imagesc([],F/1000,Ibin);title('Initial binary image');axis('xy')
    else
        imagesc([],F/1000,flipud(Iopen));title('After dynamic range horizontal');axis('xy')
        
    end
end


end



