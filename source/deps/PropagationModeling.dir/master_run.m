%%%%%master_run.m%%%%%%%%%%
%  Code to compute propagation from a horizontal planar source array
%  like that of a seismic airgun
% This version cycles through frequency, then nests azimuth, in order to
% save computer memory....

%Mar 7, made Tmin and Tmax vectors so I can plot multiple ranges easily...
%path('../../matlab_IO.dir',path)

%function master_run

clear
close all
%%case_output specifies type of acoustic output.  Categories include
%       TL_C:  coherent transmission loss
%       TL_I:  incoherent transmission loss
%       TL_S:  semicoherent transmission loss (Lloyd mirror source pattern)-BELLHOP
%       time_series:  time series
%       separate_modes: times series with separate modes broken out
%       group_velocity:  group velocity curves
%       ray:  ray tracing.

%%Settings for Arctic vertical array
% flag_suppress_io=0;
% Icase_all=-1;  %selects a scenario without menu input.  If -1, establish menu input
% case_output='TL';  %time_series, TL, or group_velocity..
% run_options.incoherent_addition=0;  %If one, add modes incoherently, and establish propagation law...
% nfreq_TL=25;  %Default number of frequency computions for TL plots

%%Settings for San Ignacio
%flag_suppress_io=1;
%Icase_all=1:4;  %selects a scenario without menu input.  If -1, establish menu input
%case_output='TL_C';  %time_series, TL, or group_velocity..
%run_options.incoherent_addition=0;  %If one, add modes incoherently.

%%Settings for Southern California glider scenario
%flag_suppress_io=0;
%Icase_all=-1;  %selects a scenario without menu input.  If -1, establish menu input
%case_output='TL_I';  %time_series, TL, or group_velocity..
%case_output='group_velocity';  %ray, time_series, TL, or group_velocity..
%run_options.max_kr=79;
%run_options.incoherent_addition=1;  %If one, add modes incoherently.

%%Settings for generic shallow water modal scenario
flag_suppress_io=0;
Icase_all=-1;  %selects a scenario without menu input.  If -1, establish menu input
%case_output='TL_I';  %time_series, TL, or group_velocity..
case_output='group_velocity';  %ray, time_series, TL, or group_velocity..
run_options.max_kr=8;
%run_options.incoherent_addition=1;  %If one, add modes incoherently.

%Set environmental variables...
setup_environment;

%%%THESE NEVER CHANGE-ARE HISTORICAL VARIABLES
case_array='single_element'; %single_element, Fairfield_1680, Kondor_3930, Ewing_20
azi=[0 ];  %Azimuth of zero means toward ship bow, 90 is directly from the side...
case_spectrum='impulse';  %impulse, Hanning, sine_wave, Bolt_airgun_80

for Icase=Icase_all
    close all;
    [case_output,savename,code_chc,case_bath, ...
        case_ssp,case_bottom,~,ro,sd,rd,D,fl,fu,fs,Tmin,Tmax,tilt]=TOC_scenario(Icase,case_output,flag_suppress_io);
    
    separate_modes=strcmp(case_output,'separate_modes');
    %To ensure that number of time series points are a power of 2,
    %  adjust the freqeuncy spacing, and thus the upper time limit...
    if strcmp(case_output,'time_series')||separate_modes||strcmp(case_output,'group_velocity')
        df = 1./max(Tmax-Tmin);
        Nsamp =   floor(fs/df);
        Nsamp=2^ceil(log10(Nsamp)/log10(2));
        df=fs/Nsamp;
        Tmax=Tmin+Nsamp/fs;
        fprintf('window length: %8.4f sec, df: %6.2f Hz\n',Nsamp/fs,df);
        F=linspace(0,fs,Nsamp+1);
        flo=round(fl*Nsamp/fs);
        fhi=round(fu*Nsamp/fs);
        freq=F(flo:fhi);
        nfreq = length(freq);
    else
        if fl==fu
            nfreq=1;
            freq=fl;
        else
            prompt={'Min (Hz)','Max (Hz)','Nfreq','log?'};
            name='Enter frequency range for TL calculations (override previous entries):';
            numlines=1;
            defaultanswer={'25','300','200','no'};
            %answer=inputdlg(prompt,name,numlines,defaultanswer);
            answer=defaultanswer;
            
            fl=eval(answer{1});fu=eval(answer{2});
            nfreq_TL=eval(answer{3});  %Default number of frequency computions for TL plots
            nfreq=nfreq_TL;
            %
            if strcmp(answer{4},'no')
                freq=linspace((fl),(fu),nfreq);
            else
                freq=logspace(log10(fl),log10(fu),nfreq);
            end
        end
    end
    
    disp(sprintf('Source depth: %6.2f, Receiver depth: %6.2f',sd,rd'));
    bath=TOC_bath(case_bath,D);
    [lthick,svp,sspopt,cmedtop,cmedbot,rho,pwavatten,plottitle,case_ssp,case_bottom,rms_roughness]=TOC_env(case_ssp,case_bottom,D,bath);
    [x,y,Sf_pos,X,Y]=TOC_array(case_array);
    
    
    
    %%%There are two ways of storing environmental sound speed profiles.
    %%% The svp variable can store a complex sound speed profile (e.g. CTD)
    %%% Otherwise, for isovelocity or linear gradient profiles cmedtop and cmedbot
    %%%  are used (see TOC_env.m)
    if ~isempty(svp)
        %%%Clean it up
        [cc,Iunq]=unique(svp(:,1));
        c_offset=interp1(svp(Iunq,1),svp(Iunq,2),sd);
    else
        c_offset=interp1([0 D]',[cmedtop(1) cmedbot(1)]',sd);
    end
    %c_offset is a single number
    
    %%if code_chc is RAM, set grid factors
    dr_factor=1;  %For RAM, how small to make the step size? dr=dr_factor*1500/freq;
    dz_factor=0.1;  % dz=dz_factor*dr;
    
    if strcmp(code_chc,'RAM')&&flag_suppress_io==0
        prompt={'dr_factor','dz_factor'};
        def={'0.2','0.05'};
        tittle='Parameters for RAM grid construction';
        answer=inputdlg(prompt,tittle,1,def);
        dr_factor=str2num(answer{1});
        dz_factor=str2num(answer{2});
    end
    
    for Id=1:size(bath,1)
        kr{Id}=zeros(nfreq,run_options.max_kr);
    end
    
    run_options.nofield=strcmpi('group_velocity',case_output);
    run_options.incoherent_addition=strcmpi('tl_i',case_output);
    run_options.case_output=case_output;
    h = waitbar(0,'Please wait...');
    
    first_time=true;
    
    for If=1:nfreq
        
        if rem(If, 50)==0,fprintf('Frequency: %6.2f Hz \n', freq(If));waitbar(If/nfreq,h);end
        
        [p,rplot,zplot,sd,Umode{If}]=compute_replicas(code_chc,case_bottom,freq(If),sspopt,lthick,cmedtop,cmedbot, rho, pwavatten, ...
            svp,bath,ro,rd,sd,dr_factor,dz_factor,fu,rms_roughness,run_options);
        
        if strcmp(code_chc,'KRAKEN')
            if isempty(Umode{If}{1}.phi)
                disp(sprintf('No modes at %6.2f Hz',freq(If)));
                continue
            end
            for Id=1:size(bath,1)
                kr{Id}(If,1:length(Umode{If}{Id}.kr))=Umode{If}{Id}.kr;
            end
        end
        nr=length(rplot);
        nz=length(zplot);
        ns=length(sd);
        if ~isempty(p)&&first_time
            if separate_modes
                TL=zeros(nz,nr,run_options.max_kr,nfreq); %TL dimensions are [rd ro sd nfreq]
            else
                TL=zeros(nz,nr,ns,nfreq); %TL dimensions are [rd ro sd nfreq] 
       
            end
            first_time=false;
        end
        if ~isempty(p)
            TL(:,:,1:size(p,3),If)=p;  %TL dimensions are [rd ro sd nfreq], or [rd ro nkr nfreq], if run_options.separate_modes
        end
        %if rem(If,500)==0,save(savename);end
    end %If frequency
    close(h)
    
    try
        switch case_output
            case 'group_velocity'
                
                compute_and_plot_group_velocity;
                clear extrct Xfft Xsource Xfft B x0
                save([savename '.mat']);
                
            case {'TL_C','TL_I','TL_S','TL_CG','TL_IG','TL_SG'}
                if size(TL,4)>1
                    disp('TL_C assumes one frequency, squeezing out source depth');
                end
                for Is=1:length(sd)
                    TLdB0=(-20*log10(abs(squeeze(TL(:,:,Is,:)))+eps));  %Assumes single frequency
                    plotTL;
                    disp('Warning!  TL_C under multiple source depths has not been checked...');
                    keyboard
                end
                
                %save([savename '_TL.mat'],'TLdB0','zplot','ro','freq');
            case {'time_series','separate_modes'}
                
                %%Create FM sweep if desired
                for Isource=2:2
                    if Isource==1  %FM Sweep
                        def1={sprintf('[%i %i]',fl,fu),'1','2','20'};
                        
                        if flag_suppress_io==0
                            prompt1={'[min_freq max_freq](Hz)','duration (s)','start time (s)', 'dB SNR'};
                            dlgTitle1='Parameters for scenario selection';
                            answer=inputdlg(prompt1,dlgTitle1,1,def1);
                            
                        else
                            close all
                            answer=def1;
                        end
                        tmp=str2num(answer{1});
                        fstart=tmp(1);fend=tmp(2);
                        tduration=str2num(answer{2});
                        tstartt=str2num(answer{3});
                        SNRR=str2num(answer{4});
                        
                        ysource=simulate_FMsweep_with_noise(fs,fstart,fend,Tmax-Tmin,tduration,tstartt,SNRR,'PSD');
                        figure
                        spectrogram(ysource,hanning(2048),round(0.7*2048),2048,fs,'yaxis')
                        title('source signature')
                        x_sweep=plot_time_series(TL,sd,rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,0,ysource,tilt);
                        
                    else  %No source spectrum
                        ysource=[];
                        Tmin0=Tmin;
                        
                        if separate_modes
                            nloops=run_options.max_kr;
                        else
                            nloops=length(sd);
                        end
                            
                          
                        for Is=1:nloops
                            Tmin=Tmin0;
                            if Is==1
                                tlimm=plot_time_series(squeeze(TL(:,:,Is,:)),sd,rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,0,ysource,tilt);
                            end
                            [x, Tmin]=plot_time_series_all(squeeze(TL(:,:,Is,:)),rplot,zplot,azi,Tmin,Nsamp,fs,freq,nfreq,flo,fhi,ysource,tilt,tlimm);
                            Tmax=Tmin+diff(tlimm);
                            xall(:,:,Is,:)=x;
                            %Quick check
                            Nfft=2048;
                            Icount=1;
                            figure
                            for I=1:size(x,1)
                                for J=1:size(x,2)
                                    try
                                    subplot(size(x,1),size(x,2),Icount);
                                    catch
                                        keyboard
                                    end
                                    Icount=Icount+1;
                                    x0=squeeze(x(I,J,:));
                                    [~,FF,TT,B] = spectrogram(x0,hanning(Nfft/4),round(.9*Nfft/4),Nfft,fs);
                                    %figure('name','Pulse spectrogram');
                                    B=10*log10(B);
                                    imagesc(TT,FF/1000,B-max(max(B)));%
                                    set(gca,'fontweight','bold')
                                    grid on
                                    axis('xy');
                                    %spectrogram(x(Ir,:),hanning(2048),round(0.9*2048),2048,fs,'yaxis')
                                    ylim([min(freq) max(freq)]/1000);
                                    
                                    colormap(jet);caxis([-20 0]);
                                    %xlabel('Time (sec)');ylabel('kHz');
                                    %title(sprintf('zs: %6.2f zr: %6.2f m rs: %6.2f km',sd(Is), zplot(I),rplot(J)/1000));
                                end
                            end
                            orient landscape
                            print('-djpeg','-r300',sprintf('%s_%i.jpg',savename,Is));
                            
                        end %sd
                        x=xall;clear xall
                    end %if Isource==1
                end  %for Isource
                
                %Save data to file
                clear extrct Xfft Xsource Xfft B x0 
                %rplot and zplot describe dimensions of x
                %rd=sd;
                save([savename '.mat']);
                
                
        end  %case_output, Icase
    catch
        disp('crash')
    end %try
    
end




%

