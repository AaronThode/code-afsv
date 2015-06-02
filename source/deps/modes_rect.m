function [y,x,zone] = modes_rect(RTF,thresh_dB,zone)

[Nf,Nt]=size(RTF);

%% Define the area of the mode we're looking for

if ~isfield(zone,'fmax') % First time looking for modes for this delay
    [~,tmp]=max(RTF(:));
    [zone.fmax,~]=ind2sub(size(RTF),tmp);
    
    if zone.fmax<zone.pek_cutoff*0.5+zone.dpek/2 % Max = mode 1
        zone.mode1=[zone.fmax-zone.dpek/2 zone.fmax+zone.dpek/2];
        zone.mode2=[zone.fmax+zone.dpek/2 zone.fmax+3*zone.dpek/2];
        zone.mode3=[zone.fmax+3*zone.dpek/2 zone.fmax+5*zone.dpek/2];
    elseif zone.fmax>zone.pek_cutoff*2.5-zone.dpek/2 % Max = mode 3
        zone.mode1=[zone.fmax-5*zone.dpek/2 zone.fmax-3*zone.dpek/2];
        zone.mode2=[zone.fmax-3*zone.dpek/2 zone.fmax-zone.dpek/2];
        zone.mode3=[zone.fmax-zone.dpek/2 zone.fmax+zone.dpek/2];
    else % Max = mode 2
        zone.mode1=[zone.fmax-3*zone.dpek/2 zone.fmax-zone.dpek/2];
        zone.mode2=[zone.fmax-zone.dpek/2 zone.fmax+zone.dpek/2];
        zone.mode3=[zone.fmax+zone.dpek/2 zone.fmax+3*zone.dpek/2];
    end
    
    zone.mode1=min(max(1,round(zone.mode1)),Nf);
    zone.mode2=min(max(1,round(zone.mode2)),Nf);
    zone.mode3=min(max(1,round(zone.mode3)),Nf);
    
    zone.next=1; % Do the 1st mode first
end

if zone.next==1
    [i_min,i_max]=deal(zone.mode1(1),zone.mode1(2));
    f_RTF=10*log10(RTF(zone.mode1(1):zone.mode1(2),:));
    
elseif zone.next==2
    [i_min,i_max]=deal(zone.mode2(1),zone.mode2(2));
    f_RTF=10*log10(RTF(zone.mode2(1):zone.mode2(2),:));
else
    [i_min,i_max]=deal(zone.mode3(1),zone.mode3(2));
    f_RTF=10*log10(RTF(zone.mode3(1):zone.mode3(2),:));
end

zone.next=zone.next+1; % Next mode

%% Define a bounding box
[M,I]=max(f_RTF(:));

[ic,jc]=ind2sub(size(f_RTF),I);

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
