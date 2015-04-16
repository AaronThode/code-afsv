% imagesc function with additional options (full screen, colorbar, 
% color axis [10th percentile -> 100th percentile], bounds and axis xy)

function [fig] = imagescFun(x,y,z,bounds,titleName)

fig=figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(x,y,z)
try
    caxis([prctile(prctile(z,15),15) prctile(prctile(z,100),100)])
catch
    
end
if ~isempty(bounds)
    ylim(bounds)
end
title(titleName)
colorbar
axis xy

end