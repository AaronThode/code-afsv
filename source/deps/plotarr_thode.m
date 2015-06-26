function plotarr_thode(ARRFIL)
close all

Narrmx = 4;
[ Arr, Pos ] = read_arrivals_asc( ARRFIL, Narrmx );

for I=1:Narrmx
    %tmp=squeeze(Arr.SrcAngle(:,I,:))';
    tmp=squeeze(Arr.delay(:,I,:))';
    
    imagesc(Pos.r.range,Pos.r.depth,tmp)
    grid on;colorbar;
    title(sprintf('Arrival: %i',I))
    pause
end