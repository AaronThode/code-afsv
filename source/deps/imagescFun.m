% imagesc function with additional options (full screen, colorbar, 
% color axis [10th percentile -> 100th percentile], bounds and axis xy)

function [] = imagescFun(x,y,z,axOrientation)
imagesc(x,y,z)
colormap(jet)
if strcmp(axOrientation,'xy')
    axis xy
else
    axis ij
end
end