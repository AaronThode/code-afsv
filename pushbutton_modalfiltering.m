% --- Executes on button press in pushbutton_modalfiltering.
function pushbutton_modalfiltering_Callback
Bfilt=[];MFP_scenario=[];
mode_filter_chc=2;  %If 1, SMS, if 2, pseudoinverse, if 0, my original crude, inefficient code
Icase_chc=5;
data_source='simulation'; %'simulation' or 'data'
env_chc='fixed'; %'fixed' (use same envrionment for all data), or 'optimized' (for inversions specific to that signal)
for Icase=Icase_chc
    save temp Icase Bfilt MFP_scenario mode_filter_chc Icase_chc data_source env_chc
    clear 
    load temp
    close all
    
    [range_guess,maxdelay,plot_chc,Nfft2,depths,tilt_offset,Maxmodes,MFP_scenario,head,Fs,minfreq,maxfreq]=load_scenario;
   
    if Icase==min(Icase_chc)
        frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
        [N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.05],Fs);
        Bfilt = firpm(N,Fo,Ao,W);
    end
    Nchan=size(data.data.x,1);
    y=zeros(size(data.data.x));
    for I=1:Nchan
        try
            y(I,:)=filtfilt(Bfilt,1,data.data.x(I,:)-mean(data.data.x(I,:)));
        catch
            keyboard
        end
    end
    
    Np=size(y,2);
    Nfft=2^ceil(log10(Np)/log10(2));
    %B=permute(B,[2 3 1]);  %time, element, freq
    B=fft([zeros((Nfft-Np)/2,Nchan); y'.*(hanning(Np)*ones(1,Nchan)); zeros((Nfft-Np)/2,Nchan)],Nfft);
    B=B.'; %(Nel, freq);
    
   
    
    for Imodel=1:length(MFP_scenario)
        
        model=load(MFP_scenario{Imodel});
        disp(sprintf('Loaded %s',MFP_scenario{Imodel}));
        
        maxdelay(4:Maxmodes)=maxdelay(3);
        
        
        Igood=zeros(1,length(head.rd));
        for I=1:length(head.rd)
            [junk,Igood(I)]=min(abs(model.rd-head.rd(I)));
        end
        
        
        Igood_source=zeros(1,length(depths));
        for I=1:length(depths)
            [junk,Igood_source(I)]=min(abs(model.rd-depths(I)));
        end
        
        model.sd=model.rd(Igood_source);
        model.Urd=model.U(Igood,:,:);
        
        model.U=model.U(Igood_source,:,:);
        model.rd=model.rd(Igood);
        
        
        %%%Match frequency bins with frequencies computed in model
        If_min=floor(minfreq*Nfft/Fs);
        If_max=ceil(maxfreq*Nfft/Fs);
        index=If_min:If_max;
        FF=(Fs/Nfft)*index;
        Iwant=zeros(length(index),1);  %Iwant will be indicies in model frequencies that match FFT bins in index
        for If=1:length(index)
            [junk,Iwant(If)]=min(abs(model.freq-FF(If)));
        end
        kk=2*pi*model.freq(Iwant)./1495;
        
        %%Compute expected group delays between mode pairs
        if any(maxdelay<0)
            for Im=1:Maxmodes
                Igg1=find(model.cgroup(Iwant,Im))>0;
                Igg1=Iwant(Igg1);
                for Im2=(Im+1):Maxmodes
                    Igg2=find(model.cgroup(Iwant,Im2))>0;
                    Igg2=Iwant(Igg2);
                    Sij(Im,Im2)=median(abs(1./model.cgroup(Igg1,Im)-1./model.cgroup(Igg2,Im2)));
                end
                
            end
        end
        
        
        windowpwr=sum(hanning(Nfft).^2)/Nfft;
        normm=1./sqrt(Nfft*Fs*windowpwr);
        %normm=1;
        depth_weights=sqrt([model.rd(1) diff(model.rd) ]');
        
        
        %%Gosh, this has been the incorrect formula for a linear tilt for a long time...
        %rd_cum=model.rd./cumsum(model.rd);
        
        rd_cum=-(model.rd-min(model.rd))./(max(model.rd)-min(model.rd));
        
        %%Depth estimator setups...
        magg{Imodel}=zeros(length(tilt_offset),length(range_guess),Maxmodes);
        talign{Imodel}=magg{Imodel};
        talign_allModeOne{Imodel}=zeros(length(tilt_offset),length(range_guess));
        talign_all{Imodel}=zeros(length(tilt_offset),length(range_guess));
        
        %%Cross-correlation has different FFT size than orignal processing, so adjust relevent variables...
        %Nxc=2*round(Fs*maxdelay(1))+1;
        Nxc=2*round(Fs*max(maxdelay))+1;
        Nfft_xc=2^ceil(log10(Nxc)/log10(2));
        If_min1=floor(minfreq*Nfft_xc/Fs);
        If_max1=ceil(maxfreq*Nfft_xc/Fs);
        index2=If_min1:If_max1;
        Iwant2=zeros(length(index2),1);  %Iwant will be indicies in model frequencies that match FFT bins in index
        for If=1:length(index2)
            [junk,Iwant2(If)]=min(abs(model.freq-Fs*index2(If)/Nfft_xc));
        end
        
        %%%Corrections in November...
        %%Create psuedo inverse matrix...
        %U=squeeze(model.Urd(:,Iwant,1)).*(depth_weights*ones(1,length(index)));  %(rd,freq)
        U=squeeze(model.Urd(:,Iwant,1:Maxmodes));
        Nrd=length(model.rd);
        Nf=length(FF);
        VVinv=zeros(Maxmodes,Nrd,Nf);
        
        for J=1:Nf
            %Weight modeled modes by element separation
            U0=squeeze(U(:,J,:)).*(depth_weights*ones(1,Maxmodes));
            
            %%Compute Pseudoinverse
            
            if mode_filter_chc==1
                VVinv(:,:,J)=U0';
            else
                HH=U0'*U0;
                VVinv(:,:,J)=HH\U0';
                
                VVinv(:,:,J)=pinv_autotol(HH,0.1)*U0';
            end
            %In case small eigenvalues need to be removed (unstable inversion)
            %[VV,EE]=eig(HH);  %dd is
            %HHinv=zeros(Maxmodes,Maxmodes);
            %for Im=1:Maxmodes
            %    HHinv=HHinv+VV(:,Maxmodes-Im+1)*VV(:,Maxmodes-Im+1)'./EE(Maxmodes-Im+1,Maxmodes-Im+1);
            %end
            
            %VVinv(:,:,J)=HHinv*U0';
            %EEall(J,:)=diag(EE);
        end
        
        %%Enter frequency-independent time delay
        kg=zeros(Nf,Maxmodes);
        for Im=1:Maxmodes
            kg(:,Im)=real(model.kr{1}(Iwant,Im))-(2*pi*model.freq(Iwant)'./max(model.cgroup(Iwant,1)));
        end
        
        for Itilt=1:length(tilt_offset)  %Each tilt offset will get a separate figure...
            Bout=zeros(Maxmodes,Nfft);
            disp(sprintf('Tilt: %6.2f',tilt_offset(Itilt)));
            tilt_matrix=exp(-1i*tilt_offset(Itilt)*rd_cum'*kk);
            
            pp=B(:,index).*tilt_matrix;  %%Correct sampled fourier data for tilt
            for J=1:Nf
                Bout(:,index(J))=squeeze(VVinv(:,:,J))*pp(:,J);  %%Mode coefficients [Maxmodes,Nfreq];
            end
            
            Bout=Bout.';
            Bout_corrected=Bout;
            
            dvec=zeros(length(range_guess),Maxmodes,length(index2));
            for Ir=1:length(range_guess)
                if rem(Ir,10)==0, disp(sprintf('%6.2f percent ranges done',100*Ir/length(range_guess)));end
                
                %%Add frequency dispersion to phase shift...
                kg_dispersed=exp(1i*range_guess(Ir)*kg);  %Note no sqrt(kr*r) correction
                Bout_corrected(index,:)=Bout(index,:).*kg_dispersed;
                %%%Unchanged from here...
                yout=2*real(ifft(Bout,Nfft));
                yout_corrected=2*real(ifft(Bout_corrected,Nfft));
                
                if plot_chc==1
                    tilt_check(Itilt)=plot_mode_filtering;
                    if length(tilt_offset)>1
                        figure(1)
                        pause(0.1);
                        disp('capturing frame');
                        MM(Itilt)=getframe(gcf);
                        
                    end
                end
                
                for Im=2:Maxmodes
                    %maxdelay(Im)=2*range_guess(Ir).*(1./model.cgroup
                    
                    if any(maxdelay<0)
                        maxdelayt(1,Im)=Sij(1,Im)*range_guess(Ir);
                    else
                        maxdelayt(1,Im)=maxdelay(Im);
                    end
                    [xc,lags]=xcov(yout_corrected(:,1),yout_corrected(:,Im),round(Fs*maxdelayt(1,Im)));  %lags are positive if array tilting towards signal
                    [magg{Imodel}(Itilt,Ir,Im),Ibest]=max(abs(xc));
                    magg{Imodel}(Itilt,Ir,Im)=sign(xc(Ibest))*magg{Imodel}(Itilt,Ir,Im);
                    talign{Imodel}(Itilt,Ir,Im)=lags(Ibest)/Fs;  %talign estimates time difference between modes 1 and Im
                    %%%Add square of time offsets between phase-shifted modal arrivals, weighted by peak of
                    %%%correlation (which should be proportional to modal excitation).
                    %%% To avoid situations where a low peak reduces error, normalize weighting factor by peak
                    %%% power (so relative weighting between modes preserved)
                    
                    if Im==2 %first round
                        magg_norm=abs(magg{Imodel}(Itilt,Ir,Im));  %Mode 12 interference weight
                        
                    end
                    
                    wght{Imodel}(Itilt,Ir,Im)=abs(magg{Imodel}(Itilt,Ir,Im))./magg_norm;
                    talign_allModeOne{Imodel}(Itilt,Ir)=talign_allModeOne{Imodel}(Itilt,Ir)+wght{Imodel}(Itilt,Ir,Im).*talign{Imodel}(Itilt,Ir,Im).^2;
                    talign_all{Imodel}(Itilt,Ir)=talign_all{Imodel}(Itilt,Ir)+wght{Imodel}(Itilt,Ir,Im).*talign{Imodel}(Itilt,Ir,Im).^2;
                    
                    if Im>1
                        for Im2=(Im+1):Maxmodes
                            %fprintf('Computing xcorr of %i and %i\n',Im,Im2);
                            if any(maxdelay<0)
                                maxdelayt(Im,Im2)=Sij(Im,Im2)*range_guess(Ir);
                                
                            else
                                maxdelayt(Im,Im2)=maxdelay(Im);
                            end
                            
                            [xc2,lags2]=xcov(yout_corrected(:,Im2),yout_corrected(:,Im),round(Fs*maxdelayt(Im,Im2)));  %lags are positive if array tilting towards signal
                            [magg2,Ibest2]=max(abs(xc2));
                            wght2{Imodel}(Itilt,Ir,Im,Im2)=abs(sign(xc2(Ibest2))*magg2/magg_norm);
                            talign2=lags2(Ibest2)/Fs;  %talign estimates time difference between modes 1 and Im
                            talign_all{Imodel}(Itilt,Ir)=talign_all{Imodel}(Itilt,Ir)+wght2{Imodel}(Itilt,Ir,Im,Im2).*talign2.^2;
                            
                        end
                    end
                    
                    %if min(maxdelay)~=maxdelay(Im)  %When estimating depth use full xc
                    [xc,lags]=xcov(yout_corrected(:,1),yout_corrected(:,Im),round(Fs*max(maxdelay)));  %lags are positive if array tilting towards signal
                    % end
                    %%Try to keep the function as even as possible (as small a phase as possible)
                    %%Try to multiply Bout directly?  But need to make sure no sidelobes are included...
                    zpad=floor(0.5*(Nfft_xc-length(xc)));
                    xc_pad=fftshift([zeros(zpad,1); hanning(length(xc)).*xc ;zeros(zpad,1)]);
                    xc_FFT=fft(xc_pad,Nfft_xc);
                    dvec(Ir,Im,:)=xc_FFT(index2);
                    
                    if plot_chc==1
                        plot_xcorr;
                    end
                end  %Im
                
                
            end %Ir
            %%%Estimate best range from each tilt computation...
            compute_depth_estimate;
            
        end %Itilt
        %%Plot summary of range estimation time lags...
        plot_range_estimates;
        
        if plot_chc==1
            if exist('MM');
                movie2avi(MM,'tiltTest');
            end
            keyboard
            for I=1:12
                if any(get(0,'child')==I)
                    figure(I)
                    orient tall
                    print(I,'-djpeg',sprintf('ModalBeamform%i',I));
                end
            end
        end
    end %Imodel
    
end  %Icase

%%%%Inner function for loading data
    function [range_guess,maxdelay,plot_chc,Nfft2,depths,tilt_offset,Maxmodes,MFP_scenario,head,Fs,minfreq,maxfreq]=load_scenario
        
        Iparam=2;
        if strcmp(data_source,'simulation')
            switch Icase
                case 1
                    data=load('Depth40m_Range1km.mat');
                    range_source=1000; %Modeled range of source...
                case 2
                    data=load('Depth40m_Range7km.mat');
                    range_source=7000; %Modeled range of source...
                case 3
                    data=load('Depth40m_Range17km.mat');
                    range_source=17000; %Modeled range of source...
                case 4
                    data=load('Depth40m_Range35km.mat');
                    range_source=35000; %Modeled range of source...
                case 5  %Tilt check..
                    data=load('../Tilt_test5m_40mSourceDepth/Depth40m_Range17km.mat');
                    range_source=17000; %Modeled range of source...
                    %data.answer={'5000:100:20000','0.1','0.025','0','2048','1:1:55','-10:2:10','2'};
            end
        else
            
        end
        
        answer=data.answer;
        range_guess=str2double(answer{1});
        maxdelay(1)=str2double(answer{2});
        maxdelay(2)=maxdelay(1);
        maxdelay(3)=str2double(answer{3});
        plot_chc=str2double(answer{4});
        Nfft2=str2double(answer{5});
        depths=str2double(answer{6});
        tilt_offset=str2double(answer{7});
        Maxmodes=str2double(answer{8});
        %%Use if statement if want to use same set of environments for all inversions...
        if strcmp(env_chc,'fixed')&&Icase==min(Icase_chc)
            MFP_scenario=data.MFP_scenario;
        elseif ~strcmp(env_chc,'fixed')
            MFP_scenario=data.MFP_scenario;
        end
        head=data.head;
        Fs=data.Fs;
        tdate_start=data.tdate_start;
        
        %%%Change parameters here
        minfreq=50;maxfreq=150;
        %tilt_offset=[-2 0];
         switch Iparam
            case 1
                tilt_offset=[-6 6 -3 3 0];
                Maxmodes=3;
            case 2  %Debug case
                Maxmodes=3;
                tilt_offset=[-5];  %Debug case...
                range_guess=range_source;
                plot_chc=1;
        end
        
    end
%%%Inner function for plotting range estimates
    function plot_range_estimates
        %%%Plot time alignment vs hypothesized source range and look for zero crossing
        persistent best_range fitt best_range2 fitt2
        
        if Imodel==1
            best_range=[]; fitt=[]; best_range2=[]; fitt2=[];
        end
        plotstr='ods+*';
        plotstr_tilt='gbrck';
        
        if length(range_guess)>1  %Plot only once...
            for Itilt2=1:length(tilt_offset)
                figure(10+Itilt2)
                
                %%Plot individual mode pair combinations
                for Imm=2:Maxmodes
                    subplot(Maxmodes+1,1,Imm-1);
                    %amplitudes=1000*abs(magg{Imodel}(Itilt2,:,Imm));
                    %[ax,h1,h2]=plotyy(range_guess/1000,1000*(squeeze(talign{Imodel}(Itilt2,:,Imm))),range_guess/1000,amplitudes);
                    
                    amplitudes=squeeze(wght{Imodel}(Itilt2,:,Imm));
                    [ax,h1,h2]=plotyy(range_guess/1000,1000*(squeeze(talign{Imodel}(Itilt2,:,Imm))),range_guess/1000,amplitudes);
                    
                    
                    set(h1,'markeredgecolor',[plotstr_tilt(Itilt2) ])
                    set(h1,'marker', plotstr(Imodel));set(ax(2),'xtick',[])
                    set(gca,'fontsize',14,'fontweight','bold');
                    set(ax(2),'fontsize',14,'fontweight','bold');
                    ylabel(ax(2),'Amplitude');
                    grid on
                    xlabel('Range guessed (km)');ylabel('lag (ms)');
                    
                    title(sprintf('legend: %s, Cross-correlation between modes 1 and %i, Tilt: %6.2f Freq range: %6.1f-%6.1f Hz', ...
                        plotstr,Imm, tilt_offset(Itilt2),minfreq,maxfreq));
                    hold on
                    line(range_guess([1 end])/1000,maxdelayt(1,Imm)*[1 1]*1000);
                    line(range_guess([1 end])/1000,-maxdelayt(1,Imm)*[1 1]*1000);
                    %if Imm==2, ylim([-25 25]);end
                    ylim(maxdelayt(1,Imm)*[-1 1]*1000);
                    ylimm=ylim;
                    set(ax(1),'ytick',linspace(ylimm(1),ylimm(2),5));
                end
                
                
                %%Plot MSE error vs. range estimate, mode 1 combos only%%%%%%%%
                [fitt(Itilt2,Imodel),Ir_chc]=min(talign_allModeOne{Imodel}(Itilt2,:));
                best_range(Itilt2,Imodel)=range_guess(Ir_chc)/1000;
                
                subplot(Maxmodes+1,1,Maxmodes);
                plot(range_guess/1000,1000*sqrt(squeeze(talign_allModeOne{Imodel}(Itilt2,:))),[plotstr_tilt(Itilt2) '-' plotstr(Imodel)]);
                set(gca,'fontsize',14,'fontweight','bold');
                grid on;xlabel('Range guessed (km)');ylabel('MSE (ms)');
                title(sprintf('legend: %s, MSE, all mode 1 combos: %6.3g ms at %6.2f km, Tilt: %6.2f Freq range: %6.1f-%6.1f Hz',plotstr, ...
                    1000*sqrt(fitt(Itilt2,Imodel)), range_guess(Ir_chc)/1000, tilt_offset(Itilt2),minfreq,maxfreq));
                hold on
                %line(range_guess([1 end])/1000,maxdelay(Imm)*[1 1]*1000);
                %line(range_guess([1 end])/1000,-maxdelay(Imm)*[1 1]*1000);
                plot(best_range(Itilt2,Imodel),1000*sqrt(fitt(Itilt2,Imodel)),[plotstr_tilt(Itilt2) 's' ],'markerfacecolor',plotstr_tilt(Itilt2));
                %ylim([0 50]);
                
                %%%Plot MSE error vs. range estimate for all mode combinations...
                subplot(Maxmodes+1,1,Maxmodes+1);
                plot(range_guess/1000,1000*sqrt(squeeze(talign_all{Imodel}(Itilt2,:))),[plotstr_tilt(Itilt2) '-' plotstr(Imodel)]);
                set(gca,'fontsize',14,'fontweight','bold');
                grid on
                xlabel('Range guessed (km)');ylabel('MSE (ms)');
                title(sprintf('legend: %s, MSE, all modes, Tilt: %6.2f Freq range: %6.1f-%6.1f Hz',plotstr, tilt_offset(Itilt2),minfreq,maxfreq));
                hold on
                [fitt2(Itilt2,Imodel),Ir_chc2]=min(talign_all{Imodel}(Itilt2,:));
                best_range2(Itilt2,Imodel)=range_guess(Ir_chc2)/1000;
                plot(best_range2(Itilt2,Imodel),1000*sqrt(fitt2(Itilt2,Imodel)),[plotstr_tilt(Itilt2) 's' ],'markerfacecolor',plotstr_tilt(Itilt2));
                
                %line(range_guess([1 end])/1000,maxdelay(Imm)*[1 1]*1000);
                %line(range_guess([1 end])/1000,-maxdelay(Imm)*[1 1]*1000);
                
            end  %Itilt2
            
        end  %if range_guess length is 1
        
        
        if Imodel==length(MFP_scenario)
            %yes=menu('Save and Print result?','Yes','No');
            yes=1;
            if yes==1
                for Itilt2=1:length(tilt_offset)
                    savename_base=sprintf('%s_%ikm_tilt%4.2f_%ito%iHz_%imodels',datestr(tdate_start,30),round(range_source), ...
                        tilt_offset(Itilt2),round(minfreq),round(maxfreq),Imodel);
                    
                    figure(10+Itilt2)
                    orient landscape
                    print(gcf,'-djpeg',sprintf('RangeEstimate_%s.jpg',savename_base));
                    
                    %%%%%Plot best range vs. fit%%%%%%
                    figure(1)
                    
                    subplot(2,1,1)
                    for Imodel2=1:length(MFP_scenario)
                        %subplot(length(MFP_scenario),1,Imodel2);
                        plot(best_range(Itilt2,Imodel2),1000*sqrt(fitt(Itilt2,Imodel2)),[plotstr_tilt(Itilt2) plotstr(Imodel2)], ...
                            'markersize',10,'markerfacecolor',plotstr_tilt(Itilt2));
                        hold on
                    end
                    ylim([0 2]);
                    set(gca,'fontweight','bold','fontsize',14);
                    xlabel('Best range (km)');
                    ylabel('MSE (msec)');
                    title(sprintf('Mode 1 contribution only: Color %s is tilt, %s is inversion model',plotstr_tilt,plotstr));grid on
                    
                    subplot(2,1,2)
                    for Imodel2=1:length(MFP_scenario)
                        %subplot(length(MFP_scenario),1,Imodel2);
                        plot(best_range2(Itilt2,Imodel2),1000*sqrt(fitt2(Itilt2,Imodel2)),[plotstr_tilt(Itilt2) plotstr(Imodel2)]);
                        hold on
                    end
                    set(gca,'fontweight','bold','fontsize',14);
                    
                    xlabel('Best range (km)');
                    ylabel('MSE (msec)');
                    title(sprintf('All mode contributions: Color %s is tilt, %s is inversion model',plotstr_tilt,plotstr));grid on
                    
                end
                
                %%%%%Plot min error vs array tilt...
                figure(2)
                
                subplot(2,1,1)
                for Imodel2=1:length(MFP_scenario)
                    plot(tilt_offset,1000*sqrt(fitt(:,Imodel2)),['k-' plotstr(Imodel2)]);hold on;grid
                    
                end
                xlabel('tilt offset (m)');
                ylabel('mse error (msec)');
                title(sprintf('Mode 1 only: %s is inversion model',plotstr));grid on
                
                subplot(2,1,2)
                for Imodel2=1:length(MFP_scenario)
                    plot(tilt_offset,1000*sqrt(fitt2(:,Imodel2)),['k-' plotstr(Imodel2)]);hold on;grid
                    
                end
                xlabel('tilt offset (m)');
                ylabel('mse error (msec)');
                title(sprintf('All modes: %s is inversion model',plotstr));grid on
                
                
                figure(1)
                orient landscape
                subplot(2,1,1)
                legend('full model','Pekeris model','Averaged Pekeris')
                print(gcf,'-djpeg',sprintf('RangeError1_%s.jpg',savename_base));
                
                figure(2)
                orient landscape
                legend('full model','Pekeris model','Averaged Pekeris')
                print(gcf,'-djpeg',sprintf('RangeError2_%s.jpg',savename_base));
                
                save(sprintf('ModalEstimateRange_%s.mat',savename_base), ...
                    'range_guess','tilt_offset','talign','talign_all','talign_allModeOne','MFP_scenario',...
                    'minfreq','maxfreq','maxdelay','plotstr','plotstr_tilt','Maxmodes', ...
                    'fitt','best_range','fitt2','best_range2','magg','wght','tdate_start','Imodel');
                
                
            end
        end
    end

%%%%%Inner function for depth estimation
    function compute_depth_estimate
        persistent depth_chc
        depth_chc=1;
        
        if Itilt==1&&Imodel==1
            %depth_chc=menu('Select depth estimation choice:','No depth estimation','Depth estimate for all ranges','Depth estimate at best mode one only range estimate');
        end
        %depth_chc=1;
        
        if depth_chc==2
            Ir_chc=1:length(range_guess);
        else
            Ir_chc=1;
            fitt=Inf;
            if length(range_guess)>1
                %dt=diff(sign(squeeze(talign{Imodel}(Itilt,:,2))));dt=[0 dt];
                [fitt,Ir_chc]=min(talign_allModeOne{Imodel}(Itilt,:));
                
            end
            if isempty(Ir_chc)
                Ir_chc=-1;
                fprintf('No best range for Model %i..\n',Imodel);
            else
                fprintf('Best range for Model %i: %6.2f km, error %6.2f msec\n',Imodel,range_guess(Ir_chc)/1000,sqrt(fitt));
            end
            if depth_chc==1
                return
            end
        end
        
        %%%Estimate source depth by relative amplitudes of cross correlations with mode 1 extraction...
        mvec=model.U(:,Iwant2,1:Maxmodes);  %All real values
        sqkr=sqrt(real(model.kr{1}(Iwant2,1:Maxmodes)));
        for Is=1:length(model.sd)
            mvec(Is,:,:)=(squeeze(model.U(Is,Iwant2,1))'*ones(1,Maxmodes)).*squeeze(mvec(Is,:,:))./sqkr;  %U1(zs)*Um(zs)/sqrt(km)
        end
        mvec(isnan(mvec))=0;
        amb{Imodel}=zeros(length(model.sd),length(Iwant2));
        amb_abs{Imodel}=amb{Imodel};
        amb_nowght{Imodel}=amb{Imodel};
        
        for Ir_best=Ir_chc
            
            dd=squeeze(dvec(Ir_best,:,:));
            for I=1:Maxmodes
                freq_wght(I,:)=abs(dd(I,:))./max(abs(dd(I,:)));
            end
            freq_wght=mean(freq_wght);  %Frequency weight for all modes
            for J=1:length(Iwant2)
                for I=1:length(model.sd)
                    mm=real(squeeze(mvec(I,J,:)));
                    Iexist=find(mm~=0);
                    mm(Iexist)=mm(Iexist)/norm(mm(Iexist));
                    %Don't penalize absent modes...
                    dd(Iexist,J)=dd(Iexist,J)/norm(dd(Iexist,J));
                    if Ir_best>0
                        amb{Imodel}(I,J)=freq_wght(J).*abs(dd(Iexist,J)'*mm(Iexist)).^2;  %FT of xcorr data, unnormalized...,
                        amb_abs{Imodel}(I,J)=freq_wght(J).*abs(dd(Iexist,J))'*abs(mm(Iexist));  %FT of xcorr data, unnormalized...,
                        
                        amb_nowght{Imodel}(I,J)=abs(dd(Iexist,J)'*mm(Iexist)).^2; %No frequency weighting
                        %amb2{Imodel}(I,J)=(mm'*dvec1).^2;   %Real/real comparison...
                    else
                        amb{Imodel}(I,J)=NaN;
                        amb_abs{Imodel}=NaN;
                        amb_nowght{Imodel}(I,J)=NaN;
                        %amb2{Imodel}(I,J)=NaN;
                    end
                end
            end
            
            amb{Imodel}=amb{Imodel}/max(max(amb{Imodel}));
            amb_abs{Imodel}=amb_abs{Imodel}/max(max(amb_abs{Imodel}));
            amb_nowght{Imodel}=amb_nowght{Imodel}/max(max(amb_nowght{Imodel}));
            
            %%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            savename_base=plot_depth_estimates;
        end  %Ir_best
        
        
        if Itilt==1&&Imodel==1
            
            yes=menu('Save result?','Yes','No');
            if yes==1
                
                
                save(sprintf('ModalEstimateDepthTilt%i_%s.mat',tilt_offset(Itilt),savename_base), ...
                    'tdate_start','range_guess','Ir_best','Iwant','tilt_offset','Itilt','amb','amb_abs','amb_nowght', ...
                    'model','Imodel','MFP_scenario','minfreq','maxfreq','maxdelay','plotstr','plotstr_tilt','Maxmodes','fitt');
                
            end
        end
        
        %%%%%%%%%Inner function plot_range_estimates%%%%%%%%%%%%
        
        %%%Inner function plot_depth_estimates%%%%%%
        function savename_base=plot_depth_estimates
            %%Plot graphs of depth estimators
            
            savename_base=sprintf('%s_%ikm_tilt%4.2f_%ito%iHz_%imodels',datestr(tdate_start,30),round(range_guess(min(Ir_best))/1000), ...
                tilt_offset(Itilt),round(minfreq),round(maxfreq),Imodel);
            
            plotstr='oxs+*';
            plotstr_tilt='kbgrc';
            
            %%%Plot 2-D depth estimator
            figure(30+10*Imodel+Itilt+Ir_best);
            
            subplot(2,1,1)
            imagesc(model.freq(Iwant),model.sd,10*log10(amb{Imodel}));
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');
            ylabel('Source depth (m)');
            title(sprintf('Depth estimate of call starting at %s, using FT of xc function',datestr(tdate_start)));
            plot_letter_label('a)');
            
            subplot(2,1,2)
            imagesc(model.freq(Iwant),model.sd,10*log10(amb_nowght{Imodel}));
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');
            ylabel('Source depth (m)');
            title('No frequency weighting');
            orient tall
            print(gcf,'-djpeg',sprintf('DepthEstimate2D_%s.jpg',savename_base));
            
            
            %Plot frequency-averaged depth estimator%
            
            figure(100+Itilt+Ir_best);
            plot1=10*log10(sum(amb{Imodel}'))/length(Iwant);
            plot_norm=10*log10(sum(amb_nowght{Imodel}'))/length(Iwant);
            plot_abs=10*log10(sum(amb_abs{Imodel}'))/length(Iwant);
            plot(plot1/max(abs(plot1)),model.sd,[plotstr_tilt(Itilt) '-' plotstr(Imodel)]);hold on
            %plot(plot_norm/max(abs(plot_norm)),model.sd,[plotstr_tilt(Itilt) '-' plotstr(Imodel+1)]);
            plot(plot_abs/max(abs(plot_abs)),model.sd,[plotstr_tilt(Itilt) '-' plotstr(Imodel+1)]);
            
            set(gca,'fontweight','bold','fontweight','bold');grid on
            ylabel('Source depth (m)');xlabel('Normalized depth response');
            try
                title(sprintf('legend: %s, Depth estimate, Range: %6.2f km, Tilt: %6.2f m, Freq range: %6.1f-%6.1f Hz', ...
                    plotstr_tilt,range_guess(Ir_best)/1000,tilt_offset(Itilt),minfreq,maxfreq));
            end
            axis('ij')
            hold on
            legend('frequency weighted','abs value');
            print(gcf,'-djpeg',sprintf('DepthEstimate_%s.jpg',savename_base));
            
            
            
            
            
        end  %Inner function plot_depth_estimates
        
    end  %inner function compute_depth_estimates
%%%%%%%%%%%%%%%%%


%%Inner function for plotting
    function tltchk=plot_mode_filtering
        persistent timmes_signal timmes_noise pwrr yes
        
        if Itilt==1
            %yes=menu('What intermediate data do you want?','Filtered spectrograms','SNR estimate');
            yes=1;
            pwrr=[];
        end
        if yes==1
            
            y1=yout(:,1);
            y2=yout(:,2);
            tltchk=y1'*y2./(norm(y1).*norm(y2));
            fprintf('Xcor for bottom sensor, tilt %6.2f: %6.4f, \n',  ...
                tilt_offset(Itilt),tltchk);
            
            figure(1)
            
            for Imm=1:3
                [Bbeam,FF1,TT] = specgram(yout(:,Imm),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
                figure(1)
                set(gcf,'units','norm','pos',[ 0.1589    0.3345    0.2216    0.6102]);
                
                subplot(3,1,Imm);
                imagesc(TT,FF1,20*log10(abs(Bbeam).*normm));%
                orient landscape
                grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
                axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);
                caxis([-150 -70])
                if Imm==1
                    %title(sprintf('Mode %i, modal beamforming only: vertical tilt: %6.2f deg',Imm, tilt));
                end
                ylabel('Hz');
                xlabel('Time (s)');
                if Imm<3
                    set(gca,'xticklabel',[]);
                    xlabel('');
                end
                set(gca,'units','norm');
                poss=get(gca,'pos');
                poss(4)=.27;
                %set(gca,'pos',poss,'xtick',0:0.25:3);
                
                text(0.1,0.8,sprintf('Mode %i',Imm),'units','norm','color','w','fontsize',14,'fontweight','bold');
                title(sprintf('Modal Filtering, Tilt: %6.2f deg, correlation: %6.4f',tilt_offset(Itilt),tltchk));
            end
            hh=colorbar('south','fontweight','bold','fontsize',14);
            set(hh,'pos',[0.1369    0.0239    0.7462    0.0234]);
            
            
            
            for Imm=1:3
                [Bbeam1,FF1,TT] = specgram(yout_corrected(:,Imm),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
                figure(2)
                subplot(3,1,Imm);
                imagesc(TT,FF1,20*log10(abs(Bbeam1).*normm));%
                orient landscape
                grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
                ylabel('Hz'); xlabel('Time (s)');
                
                axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);grid on
                caxis([-150 -70])
                if Imm==1
                    %title(sprintf('Mode %i, modal beamforming and dispersion correction: vertical tilt: %6.2f deg',Imm, tilt));
                end
                if Imm<3
                    set(gca,'xticklabel',[]);
                    xlabel('');
                    
                end
                set(gca,'units','norm');
                poss=get(gca,'pos');
                poss(4)=.27;
                
                %set(gca,'pos',poss,'xtick',0:0.25:3);
                text(0.25,110,sprintf('Mode %i',Imm),'color','w','fontsize',14,'fontweight','bold');
                title(sprintf('Phase Correction, Range: %6.2f km, Tilt: %6.2f deg',range_guess(Ir)/1000,tilt_offset(Itilt)));
            end
            colorbar('south','fontweight','bold','fontsize',14);
            
            
            
        else
            %%Estimate SNR of signal...
            
            [tmp,lags]=xcov(yout(:,1),y(1,:)');
            [junk,Imax]=max(abs(tmp));
            Imax=lags(Imax);
            
            figure(10);
            
            [Bbeam0,FF1,TT] = specgram(y(1,:),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
            subplot(2,1,1);
            
            imagesc(TT,FF1,20*log10(abs(Bbeam0).*normm));%
            orient landscape
            grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
            axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);
            xlimm=xlim;
            title('Single element');
            
            subplot(2,1,2);
            set(gcf,'units','norm','pos',[ 0.1589    0.3345    0.2216    0.6102]);
            yout=yout(Imax:end,:);
            
            [Bbeam,FF1,TT] = specgram(yout(:,1),Nfft2,Fs,Nfft2,round(0.9*Nfft2)); %B(freq,time,element);
            
            imagesc(TT,FF1,20*log10(abs(Bbeam).*normm));%
            orient landscape
            grid on;set(gca,'fontweight','bold','fontsize',14,'ycolor','k','xcolor','k','xminorgrid','on');
            axis('xy');ylim([minfreq maxfreq]);%caxis([70 120]-30);
            xlim(xlimm)
            title(sprintf('Mode 1 filter, tilt %6.2f',tilt_offset(Itilt)));
            
            if Itilt==1
                disp('Select two times for signal estimate...');
                
                timmes_signal=ginput(2);
            end
            In=round(Fs*timmes_signal(:,1));
            dn=diff(In);
            pwrr.single.signal(Itilt)=sum(y(1,In(1):In(2)).^2)/dn;  %This has been prefiltered...
            pwrr.filtered.signal(Itilt)=sum(yout(In(1):In(2),1).^2)/dn;
            
            y1=yout(In(1):In(2),1);
            y2=yout(In(1):In(2),2);
            pwrr.filtered.xcor(Itilt)=y1'*y2./(norm(y1).*norm(y2));
            
            
            fprintf('Xcor for bottom sensor, tilt %6.2f: %6.4f, power: %6.2f \n',  ...
                tilt_offset(Itilt),abs(pwrr.filtered.xcor(Itilt)),10*log10(pwrr.single.signal(Itilt)));
            
            if Itilt==length(tilt_offset)
                figure(45)
                plot(tilt_offset,10*log10(pwrr.filtered.signal));grid on
                xlabel('Tilt Range Estimate');
                ylabel('dB re 1uPa');
                title('Measured power of filtered mode 1 signal');
                keyboard
            end
            
        end
        
    end  %plot_mode_filtering..

%Inner function for plotting xcorr
    function plot_xcorr
        pause(0.5);
        MMM(Itilt)=getframe(gcf);
        figure(3)
        % for Im=1:3
        subplot(3,1,Im)
        plot(lags/Fs,xc);xlim(10*[-0.05 0.05]);
        grid on
        title(sprintf('Mode %i, leads mode 1 by %10.6f ms',Im,1000*talign{Imodel}(Itilt,Ir,Im)));
        xlabel('Relative delay (s)');
        ylabel('Cross-correlation');
        %end
    end
end  %function modal_filtering
