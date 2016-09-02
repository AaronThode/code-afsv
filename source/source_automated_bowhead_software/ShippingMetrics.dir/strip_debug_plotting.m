%%strip_debug_plotting.m%%%
close all
%Remove failed locations (only one bearings)
for JJJ=1:1 %local and combined DASAR peak detections
    
    %Review and cleanse detections...
    if JJJ==1
        Igood=find(sum(locations_ctime'>0)>1);
        Ibad=find(sum(locations_ctime'>0)<2);  %Rejections if enough rejections to fail localization
        
       % Igood=find(sum(locations_ctime'<0)==0);  %Reject localization if *any* failures...
        %Ibad=find(sum(locations_ctime'<0)>0);
        
        
        %Igood=setdiff(1:size(locations_ctime,1),Ibad);
        titstr='single DASAR';
        
    else
        %Review and cleanse detections...
        Igood=find(sum(locations_ctime_combined'>0)>1);
        Ibad=find(sum(locations_ctime_combined'>0)<2);
        titstr='combined';
    end
    disp(titstr);
    %plot original and cleaned results
    time_inc=24;
    run_options.center_dist_limit=50000;
    
    auto.locations{1}=data.locations(Ibad);
    auto.locations_ctime{1}=data.locations_ctime(Ibad,:);
    fprintf('Plotting rejects %s, choice %s\n',titstr, chc);
    
    try
        plot_movie_all_auto(DASAR_coords,[],auto,[],[],I,run_options.center_dist_limit,time_inc,goodFile);
    end
    
    title(sprintf('%i Rejected localizations at Site %i, Rule: %s, %s, File: %s',length(Ibad),I,rule,titstr,fnames.name),'interp','none','fontsize',8);
    orient landscape
    Idot=findstr(fnames.name,'.')-1;
    print('-djpeg',sprintf('%s_%s_rejected.jpg',titstr,fnames.name(1:Idot)));
    
    %Statistics on the frequency differences between detections and shipping metric...
    tmp=locations_dfmin(Ibad,:);
    tmp=tmp(:);
    Ipass= (tmp<=500);
    stats.reject.dfmin=tmp(Ipass);
    
    tmp=locations_dfmax(Ibad,:);
    tmp=tmp(:);
    Ipass= (tmp<=500);
    stats.reject.dfmax=tmp(Ipass);
    
    tmp=locations_fmin(Ibad,:);
    tmp=tmp(:);
    Ipass= (tmp<=500);
    stats.reject.fmin=tmp(Ipass);
    
    tmp=locations_fmax(Ibad,:);
    tmp=tmp(:);
    Ipass= (tmp<=500);
    stats.reject.fmax=tmp(Ipass);
    
    
    tmp=locations_solidity(Ibad,:);
    tmp=tmp(:);
    Ipass= (tmp<=500);
    stats.reject.solidity=tmp(Ipass);
    
%     figure
%     subplot(3,1,1);
%     hist(stats.reject.dfmin,0:1:100);grid on
%     xlabel('dfmin for rejected detections');
%     subplot(3,1,2);
%     hist(stats.reject.dfmax,0:1:100);grid on
%     xlabel('dfmax for rejected detections');
%     subplot(3,1,3);
%     hist(stats.reject.fmin,0:1:100);grid on
%     xlabel('duration for rejected detections');
    
    
    figure
    auto.locations{1}=data.locations(Igood);
    auto.locations_ctime{1}=data.locations_ctime(Igood,:);
    fprintf('Plotting accepted %s, choice %s\n',titstr, chc);
    
    try
        Ifalse_org=plot_movie_all_auto(DASAR_coords,[],auto,[],[],I,run_options.center_dist_limit,time_inc,goodFile);
    end
    title(sprintf('%i Retained localizations at Site %i, Rule: %s, %s File: %s',length(Igood),I,rule,titstr,fnames.name),'interp','none','fontsize',8);
    print('-djpeg',sprintf('%s_%s_retained',titstr,fnames.name(1:Idot)));
    
    %Statistics on the frequency differences between detections and shipping metric...
    
    Ifalse=Igood(Ifalse_org);
    tmp=locations_dfmin(Ifalse,:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.false.dfmin=tmp(Ipass);
    
    tmp=locations_dfmax(Ifalse,:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.false.dfmax=tmp(Ipass);
    
    tmp=locations_fmin(Ifalse,:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.false.fmin=tmp(Ipass);
    
    tmp=locations_fmax(Ifalse,:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.false.fmax=tmp(Ipass);
    
     tmp=locations_solidity(Ifalse,:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.false.solidity=tmp(Ipass);
    
%     figure
%     subplot(3,1,1);
%     hist(stats.false.dfmin,0:1:500);grid on
%     xlabel('dfmin for accepted detections');
%     subplot(3,1,2);
%     hist(stats.false.dfmax,0:1:500);grid on
%     xlabel('dfmax for accepted detections');
%     subplot(3,1,3);
%     hist(stats.false.fmin,0:1:500);grid on
%     xlabel('duration for accepted detections');
%     
    
    Itrue=setdiff(1:length(Igood),Ifalse_org);
    tmp=locations_dfmin(Igood(Itrue),:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.true.dfmin=tmp(Ipass);
    
    tmp=locations_dfmax(Igood(Itrue),:);
    tmp=tmp(:);
    Ipass= tmp<500;
    stats.true.dfmax=tmp(Ipass);
    
    tmp=locations_fmin(Igood(Itrue),:);
    tmp=tmp(:);
    Ipass=find(tmp<500);
    stats.true.fmin=tmp(Ipass);
    
    tmp=locations_fmax(Igood(Itrue),:);
    tmp=tmp(:);
    Ipass=find(tmp<500);
    stats.true.fmax=tmp(Ipass);
    
    tmp=locations_solidity(Igood(Itrue),:);
    tmp=tmp(:);
    Ipass=find(tmp<500);
    stats.true.solidity=tmp(Ipass);
    
    
%     figure
%     subplot(3,1,1);
%     hist(stats.true.dfmin,0:1:500);grid on
%     xlabel('dfmin for accepted detections');
%     subplot(3,1,2);
%     hist(stats.true.dfmax,0:1:500);grid on
%     xlabel('dfmax for accepted detections');
%     subplot(3,1,3);
%     hist(stats.true.fmin,0:1:500);grid on
%     xlabel('fmin for accepted detections');
%     
%     
%     figure
%     subplot(3,1,1);
%     plot(stats.true.dfmin,stats.true.dfmax,'kx');grid on
%     xlabel('dfmin accepted and true');ylabel('dfmax accepted but false');
%     subplot(3,1,2)
%     plot(stats.reject.dfmin,stats.reject.dfmax,'ro');grid on
%     xlabel('dfmin rejected');ylabel('dfmax rejected');
%     subplot(3,1,3);
%     plot(stats.false.dfmin,stats.false.dfmax,'gx');grid on
%     xlabel('dfmin accepted and false');ylabel('dfmax accepted but false');
    
    [N,printname,Ibin]=hist2D([stats.true.dfmin stats.true.dfmax]',1:5:100,1:5:100,{'dfmin','dfmax'});
    title('true');
    
    [N,printname,Ibin]=hist2D([stats.reject.dfmin stats.reject.dfmax]',1:5:100,1:5:100,{'dfmin','dfmax'});
    title('reject');
    
    [N,printname,Ibin]=hist2D([stats.false.dfmin stats.false.dfmax]',1:5:100,1:5:100,{'dfmin','dfmax'});
    title('false');
    
    [N,printname,Ibin]=hist2D([stats.true.dfmin stats.true.fmin]',1:5:100,0:5:500,{'dfmin','fmin'});
    title('true');
    
    [N,printname,Ibin]=hist2D([stats.reject.dfmin stats.reject.fmin]',1:5:100,1:5:500,{'dfmin','fmin'});
    title('reject');
    
    [N,printname,Ibin]=hist2D([stats.false.dfmin stats.false.fmin]',1:5:100,1:5:500,{'dfmin','fmin'});
    title('false');
    
    
    [N,printname,Ibin]=hist2D([stats.true.fmin stats.true.fmax-stats.true.fmin]',1:2:500,1:2:50,{'fmin','fmax-fmin'});
    title('true');
    
    [N,printname,Ibin]=hist2D([stats.reject.fmin stats.reject.fmax-stats.reject.fmin]',1:2:500,1:2:50,{'fmin','fmax-fmin'});
    title('reject');
    
    [N,printname,Ibin]=hist2D([stats.false.fmin stats.false.fmax-stats.false.fmin]',1:2:500,1:2:50,{'fmin','fmax-fmin'});
    title('false');
    
    
    
    [N,printname,Ibin]=hist2D([stats.true.fmin stats.true.solidity]',1:5:200,0:0.01:1,{'fmin','solidity'});
    title('true');
    
    [N,printname,Ibin]=hist2D([stats.reject.fmin stats.reject.solidity]',1:5:200,0:0.01:1,{'fmin','solidity'});
    title('reject');
    
    [N,printname,Ibin]=hist2D([stats.false.fmin stats.false.solidity]',1:5:200,0:0.01:1,{'fmin','solidity'});
    title('false');
    
    
    
end