function [y,x] = modes_rect(RTF,thresh_dB,i_pek,di)

[i_min,i_max]=deal(max(1,i_pek-ceil(di/2)),min(i_pek+ceil(di/2),length(RTF(:,1))));
RTF=RTF(ceil(i_min/2):floor(i_max/2),:);

f_RTF=10*log10(RTF);

[M,I]=max(f_RTF(:));

[ic,jc]=ind2sub(size(RTF),I);

jm=jc;
while(f_RTF(ic,jm)>M-thresh_dB && jm>1)
    jm=jm-1;
end

L=jc-jm;
H=i_max-i_min+1;

x=[i_min:i_max linspace(i_max,i_max,2*L+1) i_max:-1:i_min linspace(i_min,i_min,2*L+1)];
y=[linspace(jc-L,jc-L,H) jc-L:jc+L linspace(jc+L,jc+L,H) jc+L:-1:jc-L];
end

