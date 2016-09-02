function [locations,test]=compute_bearings(locations,station,goodFile,param,algchc,Idebug,Nsamples)

if strcmp(algchc,'sel_ratio')
    buffertime=param.morph.eq_time;
elseif strcmp(algchc,'sel_ratio_FFT')
    buffertime=param.morph.eq_time;
else
    buffertime=0;
end

    test.ang=zeros(length(locations),length(goodFile));
    test.sd=test.ang;
    
    if ~isempty(findstr('T2007',goodFile{1}))
        cal07flag=1;
        brefa_table=calibrate_bearing_Shell2007(param.calibration_dir,goodFile{1});
    else
        cal07flag=0;
    end

for J=1:length(locations),
    if rem(J,50)==0,disp(sprintf('%5.2f percent done',100*J/length(locations)));end
    locations{J}.bearing=zeros(length(goodFile),1);
    locations{J}.kappa=locations{J}.bearing;
    locations{J}.bearing_sd=locations{J}.bearing;
     locations{J}.bearing_rel=locations{J}.bearing;
    for K=1:length(goodFile)

        if ~strcmp(algchc,'spectrogram')

            ctime=locations{J}.ctime_min(K)-buffertime;
            if ctime<=0,
                locations{J}.bearing(K)=NaN;
                locations{J}.bearing_rel(K)=NaN;
                continue
            end
            
            tlen=locations{J}.Totalduration(K)+2*buffertime;
            [y,t,head]=readfile(goodFile{K},ctime,tlen,1:3,'ctime','calibrate');
            
        elseif strcmp(algchc,'spectrogram')
            %
            %             [y,start_ctime,t,head]=extract_signal_from_localization(locations,K,J,goodFile{K},'final',1:3);  %Assumes no buffertime
            %             buffertime=0;
            %         else  %extract image from spectrogram
            head=readgsif_header(goodFile{K});
            if Idebug>0
               disp(sprintf('Location %i, station %i, time %s', J, K, ctime2str(locations{J}.ctime_min(K)))); 
            end
            y=extract_image_from_localization(locations,K,J,station,goodFile{K},param,1:3,Idebug);  %Assumes no buffertime

        end

        if isempty(y),
            locations{J}.bearing(K)=NaN;
            locations{J}.bearing_rel(K)=NaN;
            continue
        end

        %algchc={'sel_ratio','sel_ratio_FFT'};
        %for KK=1:2,
        %thet0=extract_bearings(y,buffertime,param.Nfft,param.Fs,locations{J}.feature(K).robust_fmin,locations{J}.feature(K).robust_fmax,algchc);
        [thet0,kappa,sd]=extract_bearings(y,buffertime,param.Nfft,param.Fs,locations{J}.Totalfmin(K),locations{J}.Totalfmax(K),algchc,Nsamples);

        test.ang(J,K)=thet0;
        try
            test.sd(J,K)=sd;
        end
        
        %Convert measured bearing into true bearing...
        if cal07flag==0
            thet=bnorm(thet0+head.brefa);
        else
            [junk,Icol]=calibrate_bearing_Shell2007(param.calibration_dir,goodFile{K},1);
            thet= bnorm(interp1(0:360,brefa_table(:,Icol),bnorm(thet0)));
        end
        if Idebug>0,
            info_str=sprintf('Location %i, %s, DASAR %i from alg_chc %s relative bearing: %6.2f true bearing: %6.2f +/- %6.2f', ...
                J,ctime2str(locations{J}.ctime_min(K)), K,algchc,thet0,thet,sd);
            disp(info_str);
           % plot_spectrogram(full(y{1}))
           % disp(' ');

        end

        locations{J}.bearing(K)=thet;
        locations{J}.bearing_rel(K)=thet0;
        if ~isempty(kappa)
             locations{J}.kappa(K)=kappa;
             locations{J}.bearing_sd(K)=sd;
        end
    end

    %locations{J}.bearing=(locations{J}.bearing);

end



    function plot_spectrogram(xx)
        figure
        for Ishow=1:3,
            [S,F,T,PP1]=spectrogram(xx(:,Ishow),param.Nfft,round(param.Nfft*param.ovlap),param.Nfft,head.Fs,'yaxis');
            set(gcf,'units','norm','pos',[6.666e-02     4.858333e-01     2.91666e-01     3.5000e-01]);
            subplot(3,1,Ishow);
            imagesc(T,F,10*log10(abs(PP1)));axis('xy');%caxis([0 50]);
        end
        subplot(3,1,1);
        title(info_str);
        pause

    end
end





