function [frac_duplicate,Ntotal]=strip_duplicates(ctimes,tol)
ctimes_out=[];
Ndets.org=size(ctimes,1);
Nunits=size(ctimes,2);

ctimes_unique1=unique(ctimes,'rows');
Ndets.one=size(ctimes_unique1,1);


for J=1:Nunits
    ctimes_nonzero=ctimes(ctimes(:,J)>0,J);
    [single_unit_unique,Iun,Irept]=unique(ctimes_nonzero);
    frac_duplicate(J)=1-length(Iun)/length(ctimes_nonzero);
    Ntotal(J)=length(Iun);
    disp(sprintf('%i unique ctimes in unit %i, out of %i total non-zero ctimes, %6.2f percent redundancy',length(Iun),J,length(ctimes_nonzero),100*frac_duplicate(J)));
    Nsame=zeros(1,length(Iun));
    for K=1:length(Iun)
        Isame=find(Irept==K);
        Nsame(K)=length(Isame);
       % disp(sprintf('ctime: %s  appears %i times in unit %i',ctime2str(single_unit_unique(K)),Nsame(K),J));
        %if Nsame(K)>1
        %    keyboard;
        %end
    end
    
    Ndup=histc(Nsame,1:10);
    Ncat=find(Ndup(2:end)>0);
    for LL=1:length(Ncat)
       disp(sprintf('Unit %i has %i total detections, %i unique detections, and %i detections repeated %i times',J,length(ctimes_nonzero),Ndup(1),Ndup(1+Ncat(LL)),Ncat(LL))); 
    end
    
    disp('/n');

end

end