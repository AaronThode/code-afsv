function [dist,shiftt] = hausdorff_translate(img1, img2)
% HAUSDORFF_TRANSLATE Computes the Hausdorff distance between two images
%  that are gradually translated relative to each other.
%
% function dist = hausdorff_translate(img1, img2)
%
% This function computes the Hausdorff distance between two sliding BW images
% using the distance transform
%
% Params
% ------
% IN:
% img1 = The first image.
% img2 = The second image.
% OUT:
% dist = The distance.
%
% Pre
% ---
% - The images must be of the same size and shape.
% - The images must be logicals.
%
% Post
% ----
% - The distance image is returned.
%
% SeeAlso
% -------
% bwdist
%
% Examples
% --------
% Computing the similarity between two images:
% >> sim = hausdorff2(img1, img2)

% Checking the input images types
%if not(isa(img1,'logical') && isa(img2,'logical'))
%error('The input images must be logicals (binary)!');
%end

Ilim1=find(sum(img1)>0);
Ilim2=find(sum(img2)>0);

%Image 2 is doing the sliding, where move so it is just right of image 1
%object?

minIlim2=min(Ilim2);
maxIlim2=max(Ilim2);
maxgoal=max(Ilim1)+2;
Ibound(1)=maxgoal-minIlim2;

%How far do I need to move image 2 so it is just left of image 1?
mingoal=min(Ilim1)-2;
Ibound(2)=mingoal-maxIlim2;

Ibound=sort(Ibound);

%Place limits on shift, can't wrap image 2 around.

Ibound(1)=max([Ibound(1) -minIlim2+1]);
Ibound(2)=min([Ibound(2) size(img2,2)-maxIlim2-1]);



% Computing the two reciprocal distances:
shiftt=Ibound(1):Ibound(2);
if isempty(shiftt)
    shiftt=0;
end
Ns=length(shiftt);
%img2=double(full(img2));

% Computing the distance transforms of one  image:
%Igood=find(img2>0);
distTrans = bwdist((img1));

dist=Inf*ones(1,Ns);
for I=1:Ns
    
    Z=circshift(img2,[0 shiftt(I)]);
    dist(1,I) = max(max(distTrans.*Z));
    
    %Z=circshift(distTrans,[0 -shiftt(I)]);
    %dist(2,I) = max(max(Z.*img2));
   
   
end

% tmp=max(diff(dist));
% if tmp>0
%     keyboard;
% end
% dist=dist(1,:);
%disp('');

% figure;
% subplot(3,1,1)
% imagesc(img1);
% title('hausdorff_translate');
% subplot(3,1,2)
% imagesc(img2)
% subplot(3,1,3)
% plot(shiftt,dist);grid on
% %ylim([0 20])
% keyboard;

    

%dist=max(dist);