path(path,'../CommonScripts.dir')
path(path,'../BulkProcessingScripts')
run_options.quick_load=1;
Isite=2;
tol=0.2;

for Isite=1:5
    TSV_dir=sprintf('/Volumes/macmussel1/Arctic_2008/Processed/Site_0%i/Shell08_Final_February.dir/WEST_TSV_files.dir/',Isite);
    fnames=dir([TSV_dir '/*.tsv']);
    strr='ABCDEFGHIJK';
    for I=1:10,
        goodDASAR{I}=sprintf('S%i08%s0',Isite,strr(I));
    end
    DASAR_coords=load_DASAR_coords('Shell08',Isite,goodDASAR);
   
    diary redundancy.txt
    for I=1:length(fnames),
        
        % for I=15:15
        Idot=findstr('.tsv',fnames(I).name)-1;
        save_name=[fnames(I).name(1:Idot) '.mat'];

        disp(fnames(I).name);
        if run_options.quick_load==0
            [ind,localized]=read_tsv([TSV_dir '/' fnames(I).name],0,Inf,goodDASAR);
            save(save_name,'ind','localized');
        else
            load(save_name)

        end

        %keyboard;
        [frac_duplicate(I,:),Ntotal(I,:)]=strip_duplicates(ind.ctime,tol);
        
        %Conduct error check...
        Isurv=[];
        for Iloc=1:length(localized.ctev)
            Ikeep=find(ind.ctime(Iloc,:)>0);
            dt_true=ind.ctime(Iloc,Ikeep)-ind.ctime(Iloc,Ikeep(1));
            VM=[localized.utmx(Iloc) localized.utmy(Iloc)];
            [dt_err,dt_est]=dt_error_check(VM,DASAR_coords(Ikeep,:),dt_true.');
            Ikeep=Ikeep(find(abs(dt_err)<0.75|dt_err==0));
            if length(Ikeep)>1
               Isurv=[Isurv Iloc]; 
            end
        end
        
        [frac_duplicate2(I,:),Ntotal2(I,:)]=strip_duplicates(ind.ctime(Isurv,:),tol);
        Iorg_total(I)=length(localized.ctev);
        Isurv_total(I)=length(Isurv);
    end
    
    save(sprintf('redundancy_study_%i',Isite),'frac_duplicate','Ntotal','frac_duplicate2','Ntotal2','goodDASAR','Isurv_total')
end



for Isite=1:5
     strr='ABCDEFGHIJK';
    for I=1:10,
        goodDASAR{I}=sprintf('S%i08%s0',Isite,strr(I));
    end

    load(sprintf('redundancy_study_%i',Isite));
    figure
    subplot(3,2,1)
    imagesc(frac_duplicate);colorbar
    subplot(3,2,2)
    imagesc(frac_duplicate2);colorbar
   
    subplot(3,2,3)
    imagesc(frac_duplicate.*Ntotal);
    colorbar
    subplot(3,2,4)
    imagesc(frac_duplicate2.*Ntotal2);
    colorbar
    
     subplot(3,2,5)
    imagesc(Ntotal);
    colorbar
    subplot(3,2,6)
    imagesc(Ntotal2);
    colorbar
end



