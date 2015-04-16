function [y,x] = modes_rect(RTF,thresh_dB,i_pek,di)
Nt=length(RTF(1,:));
Nf=length(RTF(:,1));

[i_min,i_max]=deal(max(1,i_pek-di),min(Nf,i_pek+di));
RTF=RTF(i_min:i_max,:);

f_RTF=10*log10(RTF);

[M,I]=max(f_RTF(:));

[ic,jc]=ind2sub(size(RTF),I);

[im,iM]=deal(ic,ic);
[jm,jM]=deal(jc,jc);

while f_RTF(im,jc)>M-thresh_dB(1) && im>1
    im=im-1;
end

while f_RTF(iM,jc)>M-thresh_dB(1) && iM<i_max-i_min+1
    iM=iM+1;
end

while f_RTF(ic,jm)>M-thresh_dB(2) && jm>1
    jm=jm-1;
end

while f_RTF(ic,jM)>M-thresh_dB(2) && jM>Nt
    jM=jM+1;
end

L=ceil((jM-jm)/2);
H=ceil((iM-im)/2);

[i_inf,i_sup]=deal(max(1,i_min+ic-H),min(Nf,ic+H+i_min));
[j_inf,j_sup]=deal(max(1,jc-L),min(Nt,jc+L));

H=i_sup-i_inf+1;
L=j_sup-j_inf+1;

x=[i_inf:i_sup linspace(i_sup,i_sup,L) i_sup:-1:i_inf linspace(i_inf,i_inf,L)];
y=[linspace(j_inf,j_inf,H) j_inf:j_sup linspace(j_sup,j_sup,H) j_sup:-1:j_inf];
end

