function dist = hausdorff2(img1, img2)
% HAUSDORFF2 Computes the Hausdorff distance between two images
%
% function dist = hausdorff2(img1, img2)
%
% This function computes the Hausdorff distance between two BW images
% using the distance transform.
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
if not(isa(img1,'logical') && isa(img2,'logical'))
error('The input images must be logicals (binary)!');
end

% Computing the distancetransforms of the two images:
distTrans = bwdist(full(img1));

% Computing the two reciprocal distances:
dist = max(max(double(distTrans).*double(full(img2)))); 