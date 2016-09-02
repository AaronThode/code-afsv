function auto_corrected=update_reviewed_results_missed(auto_corrected,manual_ind,Iloc,decision,comment,Isite)
%Create a log of manual analysis of missed manual calls
persistent fname
%%Check whether we have to load a previous file
Iptr=1;
update_count=3;
if isempty(auto_corrected)|isempty(fname)
    Ichc1=menu('Select an option:','Load existing file','Create new file');
    if Ichc1==1
        fnames=dir('MissedManualReview*mat');
        for I=1:length(fnames)
            ftest{I}=fnames(I).name;
        end
        if length(fnames)>0
            Ichc2=menu('Select a file:',ftest);
            fname=ftest{Ichc2};
            data=load(fname);
            auto_corrected=data.auto_corrected;
            Iptr=size(auto_corrected.ctime,1)+1;
            
        else
            ctime=max(manual_ind.ctime(:));
            fname=sprintf('MissedManualReview_Site%i_%s.mat',Isite,datestr(ctime2str(ctime),30));
            helpdlg(sprintf('No MissedManualReview exists, creating %s',fname),'Writing dialog');
            Iptr=1;
        end
        
        
        
    else
        ctime=max(manual_ind.ctime(:));
        fname=sprintf('MissedManualReview_Site%i_%s.mat',Isite,datestr(ctime2str(ctime),30));
    end
else
    
    Iptr=size(auto_corrected.ctime,1)+1;
end
ctime=max(manual_ind.ctime(:));


auto_corrected.ctime(Iptr,:)=manual_ind.ctime(Iloc,:);
auto_corrected.flo(Iptr,:)=manual_ind.flo(Iloc,:);
auto_corrected.fhi(Iptr,:)=manual_ind.fhi(Iloc,:);
auto_corrected.duration(Iptr,:)=manual_ind.duration(Iloc,:);
auto_corrected.bearing(Iptr,:)=manual_ind.bearing(Iloc,:);
auto_corrected.decision{Iptr}=decision;
auto_corrected.comment{Iptr}=comment;
auto_corrected.location_reference=Iloc;

if ~isempty(findstr(lower(decision),'uncertain'))
    myfig=gcf;
    figure(3);
    orient tall
    print('-djpeg',sprintf('%s_%s_%i.jpg',fname,datestr(datenum(1970,1,1,0,0,ctime),30),Iptr));
    figure(myfig);
end

Iptr=Iptr+1
if rem(Iptr,update_count)==0
    helpdlg(sprintf('Writing results to %s',fname),'Writing dialog');
    save(fname,'auto_corrected');
    
end



