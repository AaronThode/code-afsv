%%%%%%%%%image_threshold.m%%%%%%%%%%%%
% Demonstrates basic idea of morphological reconstruction
% function image_threshold(I,SNRmin_peak,dynamic_range,disksize,conn,Icol)
% SNRmin_peak: minimum SNR of peak of contour section required
% dynamic_range: How many dB down from peak to preserver 'N-dB bandwidth'
% disksize: radius of disk used for initial opening.
% Icol: debug variable.  If > 0, plot cross sections of imagesc at column
% Icol
function image_threshold(I,SNRmin_peak,dynamic_range,disksize,conn,Icol)

se = strel('diamond',disksize); 

Ibad=find(I<(SNRmin_peak-dynamic_range));
I(Ibad)=0;

Iobr=Iopen_by_reconstruct(I,se);  %Remove objects too small in image and flatten peaks so a 2*disksize+1 plateau exists.

Ih= imhmax(Iobr,dynamic_range,conn);  %Remove top dynamic_range dB from peaks--thus only peaks this tall survive
Ip=(Iobr-Ih);  %Subtract plateaus from peaks to create peaks that are at most dynamic_range high--this step
               %permits weak and high SNR signals to be extracted
               %simultaneously, as long as peaks in question have
               %dynamic_range dB

%Reject 'islands' with max SNR less than SNRmin_peak

% Create a second 'h-dome' that retains only peaks that are SNRmin_peak
%  high--note this is processing original image instead of reconstructed
%  image, beacuse the process of opening reduces SNR size.
Ih2=imhmax(Iobr,SNRmin_peak,conn);
Ip2=Ip.*sign(Ih2);

%Now convert to binary.
Ibin=imregionalmax(sign(Ip2));

%Now attempt to join pieces together.
Iclose=Imclose(Ibin,strel('rect',[3 0.2/0.02]));


%Ip2=Imopen_grayscale(Ip,SNRmin_peak);

figure;
subplot(4,2,1);
imagesc(I);title(sprintf('equalized spectrogram thresholded to %6.2f',SNRmin_peak-dynamic_range));
subplot(4,2,2);
imagesc(Iobr);title(sprintf('Opened image, element size %6.2f',disksize));
subplot(4,2,3);
imagesc(Ih);title(sprintf('h-dome, %6.2f dB dynamic range',dynamic_range));
subplot(4,2,4);
imagesc(Ip);title(sprintf('Original image minus Ih'));
subplot(4,2,5);
imagesc(Ih2);title(sprintf('h-dome, %6.2f dB min peak',SNRmin_peak));

subplot(4,2,6);
imagesc(Ibin);title(sprintf('Final threshold' ));
subplot(4,2,7);
imagesc(Iclose);title(sprintf('BWareaopen, %i pixels',10 ));


if Icol>0
   figure
  % Icol=10;
   plot(I(:,Icol),'ko-');hold on
   plot(Iobr(:,Icol),'g');
   plot(Ih(:,Icol),'k--');
   plot(Ip2(:,Icol),'r');
   legend('original','opened','h-dome','final');
    
end

end

function Iobr=Iopen_by_reconstruct(I,se)

Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);

end

function If=Imopen_grayscale(I,SNRmin_peak)
neighbor=[1 1 1; 1 1 1; 1 1 1];
se=strel(neighbor,SNRmin_peak*ones(size(neighbor)));
Ie=imerode(I,se);
Ibad=find(Ie<0);
Ie(Ibad)=0;
If=imdilate(Ie,se);
end
