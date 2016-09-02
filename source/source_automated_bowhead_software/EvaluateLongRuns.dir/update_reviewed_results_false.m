function auto_corrected=update_reviewed_results_false(auto_corrected,current_location,Ilocation,decision,comment,Isite)
persistent fname
%%Check whether we have to load a previous file
Iptr=1;
update_count=3;
if isempty(auto_corrected)|isempty(fname)
    Ichc1=menu('Select an option:','Load existing file','Create new file');
    if Ichc1==1
        fnames=dir('FalseAlarmReview*mat');
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
            ctime=max(current_location.ctime_min);
            fname=sprintf('FalseAlarmReview_Site%i_%s.mat',Isite,datestr(ctime2str(ctime),30));
            helpdlg(sprintf('No FalseAlarmReview exists, creating %s',fname),'Writing dialog');
            Iptr=1;
        end
        
       
        
    else
        ctime=max(current_location.ctime_min);
        fname=sprintf('FalseAlarmReview_Site%i_%s.mat',Isite,datestr(ctime2str(ctime),30));
    end
else
        
    Iptr=size(auto_corrected.ctime,1)+1;
end
 ctime=max(current_location.ctime_min);
   


auto_corrected.ctime(Iptr,:)=current_location.ctime_min';
auto_corrected.flo(Iptr,:)=current_location.Totalfmin';
auto_corrected.fhi(Iptr,:)=current_location.Totalfmax';
auto_corrected.duration(Iptr,:)=current_location.Totalduration';
auto_corrected.bearing(Iptr,:)=current_location.bearing';
auto_corrected.decision{Iptr}=decision;
auto_corrected.comment{Iptr}=comment;
auto_corrected.location_reference=Ilocation;

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
    
    
    
