function [Ksout, Ns, EE_sort, VV] = extractKsbest_contour(app, x, ovlap, Nfft, chann, frange, fr, fbad, Fs, M, keep_zeros, nowin)
%%%%%%%%%%%%%%%%%%%%%%%extractKsbest_contour.m%%%%%%%%%%%
% [Ksout,Ns]=extractKsbest_contour(x,ovlap,Nfft,goodel,frange,fr,fbad,Fs,M,nowin);
%  Aaron Thode
%  April 2, 2004
%  Generates averaged cross-spectral outputs
%    from an FM contour
% INPUT:
%x=array of data, rows are time, columns are channels
%  ovlap=overlap of time samples
%  Nfft-number of points used in FFT
%  chann-vector containing element indicies: referenced to *bottom* element
%  frange-Vector of frequencies: first is initial start of contour,
%   second is the start of second harmonic, etc.
%  fr-Search space in terms of +-freq
%  fbad-frequencies of constant interference, etc.
% M time window
% nowin-if exists, don't window the data before using fft.used for source signatureestimates
% keep_zeros: if exists, keep frequency bins that have no samples.  Useful for plotting...
% OUTPUT:
%    Ksout: structure array containing
%     Kstot CSDM size(Nel,Nel,Nfreq)
%     freq: frequencies corresponding to Ks in third dimension
%     fcontour: frequencies detected, (Is,If)--can use to make contours..
%     fcount: Number of times each frequency has been averaged..
%     fpower: Power in each bin
%The replica must also be conjugated to work properly with Kraken.
%	Don't know why.
% Feb 9-remove bad elements before summing power
% April 1, 2004-normalize by Fs*Nfft to put units as power spectral density
EE_sort=[];VV=[];
figure;
if ~exist('nowin', 'var'),
    nowin=0;
elseif  nowin==1
    nowin=1;
else
    nowin=0;
end

Nel=length(chann);
if ~exist('M', 'var')||M<0, M=Nfft;end

%Compute number of time snapshots
Ns=floor(((size(x,1)/M)-ovlap)/(1-ovlap));
disp(['Ns: ' int2str(Ns)]);
%if Ns==0,Ns=1;end

f=linspace(0,Fs,Nfft+1);
df=diff(f);df=df(1);
nbins=ceil(fr/df);
bins=-nbins:nbins;

fcount=zeros(size(f));
fpower=fcount;
Kstot=zeros(Nel,Nel,length(f));
%Kstot=zeros(Nel,Nel,length(findex));

%Select appropriate frequency bin
for I=1:length(frange),
    [junk,findex(I)]=min(abs(f-frange(I)));
    frange(I)=f(findex(I));
end

for I=1:length(fbad)
    [tmp,findexjunk(I)]=min(abs(f-fbad(I)));
    fbad(I)=f(findexjunk(I));
end

if Ns<0,
    disp('Signal too short for one shapnot, will center pad for FFT:');
    %pause;
    Ns=1;Nx=size(x,1);
    x0=zeros(Nfft,size(x,2));
    index=floor(Nfft/2-(Nx/2));
    index=(index:(index+Nx-1));
    x0(index,:)=x;
    x=x0;
    clear x0
    M=size(x,1);
end
if nowin==0,
    win=kaiser(M,2.5);
else
    win=ones(M,1);
end

%Determine the frequency with greatest average power near your bin!
Pt=[];
t=[];
for I=0:(Ns-1),
    index=round(I*M*(1-ovlap)+1);
    t(I+1)=index(1)/Fs;
    xindex=(index:(index+M-1));
    xh=x(xindex,chann);
    for Ic=1:size(xh,2),
        xh(:,Ic)=xh(:,Ic)-mean(xh(:,Ic));
        xh(:,Ic)=xh(:,Ic).*win;
    end
    Xh=fft(xh,Nfft);
    Pwr=abs(Xh).^2;
    Pt=sum(Pwr,2);
    if ~isempty(fbad)
        Pt(findexjunk)=0; %Remove bad freqencies.
        Pt(findexjunk+1)=0;
        Pt(findexjunk-1)=0;
    end
    %Pt=cat(2,Pt,Pwr);
    %end
    %Pt=sum(Pt,2);
    for If=1:length(findex),
        subplot(length(findex),1,If);

        %plot(f(findex(If)+bins),Pt(findex(If)+bins));
        [junk,fi]=max(Pt(findex(If)+bins));
        findex(If)=findex(If)+bins(fi);
        fcount(findex(If))=fcount(findex(If))+1;
        fcontour(I+1,If)=f(findex(If));
        fpower(findex(If))=fpower(findex(If))+Pt(findex(If));
        disp(f(findex(If)));
        pgoal=Xh(findex(If),:);
        %Make pgoal a vertical array
        %pgoal=conj(pgoal);   %Test to see if conjugation is the problem
        %   PREVIOUS STATEMENT COMMENTED OUT BY A THODE MARCH 15, 2004
        %     HE ALSO CHANGED write_covmat.m to remove conjugation as well

        pgoal=pgoal.';        %Rotate so vector is vertical, like KRAKEN
        Kstot(:,:,findex(If))=Kstot(:,:,findex(If))+pgoal*pgoal';
        %Ksout.pgoal
    end
    %title(int2str(Ns));
    %pause(0.25);
    disp('');
end

%%Collapse Kstot to non-zero components if desired
Igood=find(fcount>0);

if exist('keep_zeros', 'var') %%Keep all frequency bins, even if no power...
    %Igood=[min(Igood):max(Igood)];
    Igood=1:length(fcount);
end
fcount=fcount(Igood);
fpower=fpower(Igood);
freq=f(Igood);
Kstot=Kstot(:,:,Igood);

%Normalize by sample size, if necessary
Iavg=find(fcount>1);
for I=1:length(Iavg),
    Kstot(:,:,Iavg(I))=Kstot(:,:,Iavg(I))/fcount(Iavg(I));
    fpower(Iavg(I))=fpower(Iavg(I))/fcount(Iavg(I));
end

%Normalize the |X(f)|^2 term to make units power spectral density (power
%per Hz)
Kstot=Kstot/(Fs*Nfft);
%keyboard;

if Ns>1
    for J=1:length(freq)
        Ks=squeeze(Kstot(:,:,J));
        [V,D,FLAG]=eigs(Ks,4,'LM');
        %keyboard;
        D=real(diag(D));
        %%sqrt(D(1)*V(:,1)) gives scaled eigenvector
        %D_sort=sort(D);
        %D=flipud(D_sort)
        pgoal(:,J)=V(:,1)*sqrt(abs(D(1)));
        SN(1,J)=10*log10(abs(D(1)/D(2)));
        if fcount(J)>0
            disp(['Est. S/N for ' num2str(freq(J)) ' is ' num2str(SN(1,J)) 'dB']);
        end
    end
end

Ksout.Kstot=Kstot;
Ksout.freq=freq;
Ksout.fcount=fcount;
Ksout.fpower=fpower;
Ksout.fcontour=fcontour;
Ksout.tcontour=t;
Ksout.pgoal=pgoal;
Ksout.SN=SN;

close;
end
