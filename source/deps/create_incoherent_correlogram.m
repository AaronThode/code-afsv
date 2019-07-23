function [mean_corr,tindex,TT_plot,pwr,pwr_tot]= create_incoherent_correlogram(TT,FF,B,param,flo,fhi)
%function [mean_corr,tindex,TT_plot,pwr,pwr_tot]= create_incoherent_correlogram(TT,FF,B,param,flo,fhi)
%% TT,FF, B:  B [Nfreq Ntime] is a linear spectrogram matrix, TT and FF are time and frequency axes.
% param.ici_range;  %Autocorrelation lag to examine
% param.time_sample;  %sec,  amount of spectrogram time columns to process for autocorrelation..
% param.ovlap;  %How much to shift the autocorrelation window between samples
% param.teager;  %Apply teager-kaiser operation to spectrogram before processing...
% flo, fhi:   min and max frequency in Hz

ici_range=param.ici_range;  %Autocorrelation lag to examine
time_sample=param.time_sample;  %sec,  amount of spectrogram time columns to process for autocorrelation..
ovlap=param.ovlap;  %fraction of how much to shift the autocorrelation window between samples (e.g. 0.5 means 50% overlap)
teager=param.teager;  %Apply teager-kaiser operation to spectrogram before processing...

% A single event may thus persist for up to time_sample/dX times, 1/(1-ovlap) times



%freq=flo:bandwidth:fhigh;
B=10*log10(B+eps);
if teager
    B=B(:,2:end-1).^2-(B(:,1:end-2).*B(:,3:end));
end
dT= TT(2)- TT(1);
dF= FF(2)- FF(1);
Ncol=length(find( TT-TT(1)<=time_sample));  %Number of columns to process at a time:
Iindex=1:ceil((1-ovlap)*Ncol):length( TT);
tindex=dT*Iindex;  %Time axis of absolute time
dX=tindex(2)-tindex(1);  %X axis of new correlation image
Ntime=length(Iindex)-1;
%yscale=dT*(0:(Ncol-1));

Igood=find( FF>flo& FF<fhi);
maxlag = Ncol - 1;  %maximum autocorrelation lag
laggs=-maxlag:maxlag;
I_range=find(laggs*dT>ici_range(1)&laggs*dT<ici_range(2));
TT_plot=laggs(I_range)*dT;
mean_corr=zeros(length(I_range),Ntime);
pwr=zeros(Ntime,Ncol);
pwr_tot=zeros(Ntime);
xcovv=zeros(length(Igood),length(I_range));

for I=1:Ntime-5
    if rem(I,1000)==0,disp(sprintf('%6.2f percent done',100*I/Ntime));end
    
    try
        A= B(Igood,Iindex(I)+(0:(Ncol-1)));  
        
        [N,M] = size(A);  %Autocorrelate columns
        X=fft(A'-ones(M,1)*mean(A'), 2^nextpow2(2*M-1));
        tmp=ifft(abs(X).^2);
        
        tmp1=[tmp(end-maxlag+1:end,:);tmp(1:maxlag+1,:)];  %Rearrange and trim to maxlag
        tmp=tmp1./(ones(size(tmp1,1),1)*tmp1(maxlag+1,:));  %Normalize by value at zero lag
        
        xcovv=tmp(I_range,:);  %xcovv has dimentions of [ laggs, frequency bin]
        %end
        mean_corr(:,I)=median(xcovv')';  %Get median value across frequency
        pwr(I,:)=sum(A)/length(Igood);
        pwr_tot(I)=max(pwr(I,:));
    catch
        
        mean_corr(:,I)=[];
    end
    %     figure(2);
    %     subplot(2,1,1)
    %     imagesc([], FF(Igood)/1000,A);
    %     subplot(2,1,2)
    %     imagesc(TT_plot, FF(Igood)/1000,xcovv{I});
    %     caxis([0 1]);
    %     colorbar
    %     pause(0.25);
    
    
end  %Ntime
end
