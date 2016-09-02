%%quick_review_links.m%%%
close all
cred_limit1=0;
cred_limit2=0;

Nhits=sum(sign(loc_index.I));
link_frac=loc_index.credibility_score./(Nhits.*(Nhits-1));

Igood1=find(sum(loc_index.I)>0);
Igood=find(loc_index.credibility_score>=cred_limit1&loc_index.credibility_score<=cred_limit2);
disp(sprintf('%i out of %i pass',length(Igood),length(Igood1)));
for J=Igood
    clf
    for I=1:7
        subplot(1,7,I)
        II=loc_index.I(I,Igood(J));
        try
            IM=station(I).Image{II};
            Ilim=find(sum(IM)>0);
            imagesc(IM(:,Ilim));
        end
    end
    pause
end


for I=1:7
    xvec=loc_index.I(I,:);
    xvec=xvec(1,xvec>0);
    
    disp(sprintf('Out of %i elements %i are unique',length(xvec),length(unique(xvec))));
end