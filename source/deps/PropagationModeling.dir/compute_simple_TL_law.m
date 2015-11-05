%function [prop_law,prop_law_mean,TLmeandB,TL_est]=compute_simple_TL_law(rplot,TLdB,freq,nfreq)
%
%  Input: TLdB: [range by nfreq];
% Output:
%   prop_law:  alpha value for each frequency column of TLdB
%   prop_law_mean:  broadband alpha value all frequencies weighted equally
%   TLmeandB: model transmission loss (broadband) vs. rplot
%   TL_est: polyval curve fit of best alpha value to rplot

function [prop_law,prop_law_mean,TLmeandB,TL_est]=compute_simple_TL_law(rplot,TLdB,freq,nfreq)

plot(log10(rplot),TLdB);

TL=10.^(-TLdB/20);
if nfreq>1
    
    TLmean=mean(TL,2);
else
    TLmean=TL;
end
TLmeandB=-20*log10(abs(TLmean));

xlabel('log(R)');
prop_law=zeros(1,nfreq);
for III=1:nfreq
    if nfreq==1
        mm = polyfit(log10(rplot),TLdB,1);
        
    else
        
        mm = polyfit(log10(rplot'),TLdB(:,III),1);
    end
    TL_est=polyval(mm,log10(rplot));
    hold on
    plot(log10(rplot),TL_est,'r--','linewidth',3);
    %disp(sprintf('Simple law, %6.2f Hz: TL=%6.2f *log(R)',freq(III),mm(1)));
    %pause
    prop_law(III)=mm(1);
end
if nfreq>1
 
    mm = polyfit(log10(rplot'),TLmeandB,1);
end
TL_est=polyval(mm,log10(rplot));
hold on
plot(log10(rplot),TL_est,'k','linewidth',5);
disp(sprintf('Simple law, arithmetic mean: TL=%6.2f *log(R)',mm(1)));
%pause
prop_law_mean=mm(1);
hold off