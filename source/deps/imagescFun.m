% imagesc function with additional options (full screen, colorbar, 
% color axis [10th percentile -> 100th percentile], bounds and axis xy)

function [] = imagescFun(x,y,z,axOrientation)
%%Assume z is in dB
z=z-max(max(z));
imagesc(x,y,z)
colormap(jet)
colorbar
caxis([-30 0]);
if strcmp(axOrientation,'xy')
    axis xy
else
    axis ij
end
end