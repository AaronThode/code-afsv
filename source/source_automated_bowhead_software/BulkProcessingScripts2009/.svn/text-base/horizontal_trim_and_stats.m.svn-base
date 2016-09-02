function [detect,stats]=horizontal_trim(detect_in,dynamic_range_h,F,Idebug,Iplot)
detect=detect_in;
stats=[];
[labeled_raw,numObjects]=bwlabel(detect_in,8);
mindet=min(min(detect_in));
for I=1:numObjects
    index=find(labeled_raw(:)==I);
    
    mask=(detect_in-mindet).*ismember(labeled_raw,I);
   
    %%Trim mask
    Itrim=find(sum(abs(mask))>0);
    %T1=T(Itrim);
    mask=mask(:,Itrim);
    
    %%Statistics of distribution...
    mask_freq=mask.*(F*ones(1,size(mask,2))); %Weight by frequency
    mean_F=sum(mask_freq)./sum(mask);


    %Local bandwidth estimate
    mask3=mask.*((F*ones(1,length(mean_F))-ones(length(F),1)*mean_F).^2);
    local_bandwidth=2*sqrt(sum(mask3)./sum(mask));
    
    mask4=mask.*((F*ones(1,length(mean_F))-ones(length(F),1)*mean_F).^4);
    fourth_moment=(sum(mask3)./sum(mask));
    local_kurtosis=fourth_moment./((local_bandwidth/2).^4);
    stats(I).mean_F=mean_F;
    stats(I).local_bandwidth=local_bandwidth;
    stats(I).local_kurtosis=local_kurtosis;
    stats(I).mean_SNR=max(mask+mindet);
    %Final trim..dumps mask
    maxval=max(detect_in(index));
    Ilow=find(detect_in(index)<=maxval-dynamic_range_h);
    detect(index(Ilow))=0;

end



if Idebug>1
    figure(Iplot);
    subplot(3,1,2)
    imagesc(detect_in)
    title('After vertical thresholding')
    subplot(3,1,3)
    imagesc(detect);
    title('After horizontal thresholding');

end