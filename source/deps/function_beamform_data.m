function handles=function_beamform_data(app,handles)
tdate_start=handles.tdate_start;
tlen=handles.tlen;
mydir=pwd;
Ichan='all';  %Hardwire first channel
%set(handles.togglebutton_ChannelBeam,'String','Channel');
[x,t,Fs,~,~,head]=load_data(handles.filetype, tdate_start,tlen,Ichan,handles,app);

%If not multichannel data, return
if isempty(x)||~head.multichannel
    uiwait(msgbox('CSDM requires multichannel data or non-zero x value'));
    return
end
%Will need to have columns be channels, rows time
if floor(tlen*Fs)~=size(x,1)
    x=x.';
end
%%For MDAT file, channel 1 is already shallowest channel, so remove data below...

Fs=round(Fs);
%disp(sprintf('Fs=%i',Fs));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

contents=get(handles.popupmenu_ovlap,'String');
ovlap=str2double(contents{get(handles.popupmenu_ovlap,'Value')})/100;
ovlap=min([1-1/Nfft ovlap]);

fmin=str2double(get(handles.edit_fmin,'String'));
fmax=str2double(get(handles.edit_fmax,'String'));
%chann=input('Enter channel indicies (remember that channels have been flipped, 1 is now shallowest) (1:Nchan):');
%if isempty(chann)
chann=1:head.Nchan; %We now assume that any bad channels are marked in the LOG files.
%end

Ichc=menu('Enter type of frequency processing:','Contour', ...
    'All (averaged)','Individual Frames (and ray tracing)', ...
    'Time delay and sum','Cross correlation between two phones');

delay_and_sum_flag=false;  %Boolean variable is true if delay and sum option selected
correlation_flag=false;
if Ichc==1
    Nf=input('Enter number of harmonics:');
    if isempty(Nf)
        Nf=1;
    end

    ftmp=ginput(Nf);
    frange=ftmp(:,2)*1000;

    Nbad=input('Enter number of bad frequencies: ');
    if ~isempty(Nbad)
        fbad=ginput(Nbad);
        fbad=fbad(:,2)*1000;
    else
        fbad=[];
    end
    fr=input('Enter frequency band to search (10):');
    if isempty(fr)
        fr=10;
    end

    [Ksout,Ns]=extractKsbest_contour(app, x,ovlap,Nfft,chann, frange,fr,fbad,Fs,Nfft,'keep_zeros');
    %    Ksout: structure array containing
    %     Kstot CSDM size(Nel,Nel,Nfreq)
    %     freq: frequencies corresponding to Ks in third dimension
    %     fcontour: frequencies detected, (Is,If)--can use to make contours..
    %     fcount: Number of times each frequency has been averaged..
    %     fpower: Power in each bin

    %[Kstot,Ks_eig,Nsnap,freq_all,pwr_est,SNRest]=extractKsexact_MDAT(x,ovlap,Nfft,chann,1000*[fmin fmax],Fs,-1,Nfft);

    axes(handles.axes1)
    hold on
    plot(Ksout.tcontour,Ksout.fcontour/1000,'o');
    hold off
elseif Ichc==2  %extractKsexact
    %%process all frequencies over a given range
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);

    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);


    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Threshold (dB)', };
    def1={'-Inf'};
    answer=inputdlg(prompt1,'CSDM threshold?',1,def1);

    threshold=eval(answer{1});

    [Ksout.Kstot,Ksout.freq,Ksout.VV,Ksout.EE]=extractKsexact(x,ovlap,Nfft,chann,frange,Fs,-1,Nfft,0,threshold);
    if ~isempty(Ksout.EE)
        %Ksout.SNR=Ksout.EE(1,:)./Ksout.EE(2,:);
        Ksout.SNR=Ksout.EE(1,:)./sum(Ksout.EE(2:end,:));
    else
        Ksout.SNR=Inf*ones(size(Ksout.freq));
    end

    figure
    plot(Ksout.freq,10*log10(Ksout.SNR),'-x');ylabel('dB SNR');xlabel('Freq (Hz)')
    grid on;title('estimated SNR of CSDM frequency components: first eignvalue/sum(rest)')

    %%%%%%%%%Extract individual frames and ray trace%%%%%%%%%
elseif Ichc==3 %%Extract frames and ray trace
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);

    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);


    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Threshold (dB)', 'optional FFT for longer FFT of small segments'};
    def1={'-Inf',''};
    answer=inputdlg(prompt1,'CSDM threshold?',1,def1);

    threshold=eval(answer{1});
    if isempty(answer{2})
        Nfft2=Nfft;
    else
        Nfft2=eval(answer{2});
    end

    [Ksout.Kstot,Ksout.freq,Ksout.t]=extractKsframes(x,ovlap,Nfft2,chann,frange,Fs,Nfft,0,threshold);

elseif Ichc==4 %% simple time and delay on filtered signal...

    delay_and_sum_flag=true;
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);

    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);


    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Min Freq (Hz)', 'Max Freq (Hz)'};
    def1={num2str(frange(1)),num2str(frange(2))};
    answer=inputdlg(prompt1,'Check frequency selection',1,def1);

    frange=[eval(answer{1}) eval(answer{2})];

    bpFilt = designfilt('bandpassiir', 'FilterOrder', 20, ...
        'HalfPowerFrequency1', frange(1), 'HalfPowerFrequency2', frange(2),...
        'SampleRate', Fs);

    %fvtool(bpFilt);
    x=filter(bpFilt,x);

elseif Ichc==5 %% simple cross correlation on filtered signal...

    correlation_flag=true;
    disp('Select bounding box for time and range:');
    ftmp=ginput(2);

    %Trim x...
    Igood=( t>=min(ftmp(:,1))&t<=max(ftmp(:,1)));
    x=x(Igood,:);
    t=t(Igood);

    frange=sort(ftmp(:,2)*1000);
    fprintf('frange: %6.2f to %6.2f Hz\n',frange);
    prompt1={'Min Freq (Hz)', 'Max Freq (Hz)','Max delay (msec)', ...
        'scale option [biased, unbiased, coeff]'};
    def1={num2str(frange(1)),num2str(frange(2)),num2str(1000*2*head.geom.spacing(2)/1472),'unbiased'};
    answer=inputdlg(prompt1,'Check Xcorr selections',1,def1);

    if isempty(answer)
        return
    end
    frange=[eval(answer{1}) eval(answer{2})];
    maxlag=ceil(eval(answer{3})*Fs/1000);
    scaleopt=answer{4};
    %maxopt=answer{5};

    bpFilt = designfilt('bandpassiir', 'FilterOrder', 20, ...
        'HalfPowerFrequency1', frange(1), 'HalfPowerFrequency2', frange(2),...
        'SampleRate', Fs);

    %fvtool(bpFilt);
    x=filtfilt(bpFilt,x);

end


if ~delay_and_sum_flag&~correlation_flag
    yes=menu('Beamform?','No','Conventional','MV','Both','Reflection Coefficient Estimation','MFP');
else
    yes=2;
end
yes_eigen=0;
while yes>1

    if yes<6&~correlation_flag

        prompt1={'Vector of angles (deg) [(-20:0.2:20)]','Vector of sin angles (deg) [(-0.2:0.005:0.2)]', ...
            'hydrophone indicies [all]','sound speed (m/sec)', 'matched filter? (1=yes, 0=no)' ...
            };
        def1={'-10:0.2:10', '-0.3:0.005:0.3',sprintf('[1:%i]',length(chann)),'1480','1'};

        answer=inputdlg(prompt1,'Beamforming parameters',1,def1);
        try
            if isempty(answer{2})
                angles=eval(answer{1});
            else
                angles=180*asin(eval(answer{2}))/pi;

            end
            fprintf('Angular range is %i points between %6.2f to %6.2f degrees\n',length(angles),min(angles),max(angles));
            Igood_el=eval(answer{3});
            cc=eval(answer{4});

            %matched filter--Ludovic
            if eval(answer{5})==1
                prompt1={'Signal type [Chirp]:','Frequency range (Hz):', ...
                    'Duration (sec):'};
                def1={'chirp','[3000 9000]','1'};

                answer=inputdlg(prompt1,'Beamforming parameters',1,def1);
                ff=eval(answer{2});
                tau=eval(answer{3});
                matched_filter=chirp(0:1/Fs:tau,ff(1),tau,ff(2))';
            end

        catch
            errdlg('Could not understand your beamforming parameters');
            return

        end
    end

    %%Assign beam_str value based on earlier choices
    switch yes
        case 2
            beam_str='CV';

            if delay_and_sum_flag
                beamchc='delaynsum';
            elseif correlation_flag
                beamchc='timedomaincorrelation';
            else
                Ndim=ndims(Ksout.Kstot);

                if Ndim==3
                    beamchc='freqvspwr';
                elseif Ndim==4
                    beamchc='angvstime';  %migration
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%  Main beamforming decision loop  %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if strcmp(beamchc,'timedomaincorrelation')
                %%Simple time domain cross-correlation
                [cc,lags]=xcov(x(:,1),x(:,2),maxlag,scaleopt);
                figure('Name',sprintf('%s Xcov result',scaleopt));
                for I=1:2
                    subplot(3,1,I)
                    plot(t,x(:,I));
                    grid on;title(sprintf('Signal %i',I));
                end
                subplot(3,1,3);
                plot(1000*lags/Fs,cc,'k',1000*lags/Fs,abs(hilbert(cc)),'r');
                [maxx,Ishifth]=max(abs(hilbert(cc)));

                [maxx,Ishift]=max(cc);


                grid on;xlabel('Delay (msec)');ylabel('Amplitude');
                title(sprintf('Max lag is %8.4f msec, %8.4f msec hilbert',1000*lags(Ishift)/Fs,1000*lags(Ishifth)/Fs));


            elseif strcmp(beamchc,'freqvspwr')
                B=conventional_beamforming(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);

                figure(20);
                imagesc(Ksout.freq,[],10*log10(abs(Ksout.EE)));xlabel('frequency (Hz)');ylabel('Eigenvalue');colorbar

                Ieig=input('Pick eigenvector [return yields none]:');
                if ~isempty(Ieig)
                    V1=squeeze(Ksout.VV(:,Ieig,:));
                    for If=1:length(Ksout.freq)
                        Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';
                    end
                    B_eig=conventional_beamforming(Ksout.Kstot_eig(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
                    yes_eigen=1;
                else
                    close(20);
                end

                ray_angles=plot_beamforming_results(true); % peakpicking

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%Individual snapshots and ray tracing%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            elseif strcmp(beamchc,'angvstime') %Individual snapshots
                Bsum=zeros(length(angles),length(Ksout.t));
                df=Ksout.freq(2)-Ksout.freq(1);
                %plot(Bsum,angles,'k');

                for Isnap=1:length(Ksout.t)
                    if rem(Isnap,10)==0,disp(Isnap);end
                    B=conventional_beamforming(squeeze(Ksout.Kstot(Igood_el,Igood_el,:,Isnap)),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
                    %Bsum(:,Isnap)=10*log10(sum(abs(B)))/length(Ksout.freq);
                    Bsum(:,Isnap)=(sum(10*log10(abs(B))))/length(Ksout.freq);
                end

                %plot vs sin angle

                figure
                imagesc((Ksout.t-min(Ksout.t))*1000,sin(angles*pi/180),Bsum)
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('sine of Elevation angle');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                set(gca,'xtick',0:5:xlimm(2));

                %plot migration angle...
                figure(21)
                imagesc((Ksout.t-min(Ksout.t))*1000,angles,Bsum)
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('Elevation angle (deg)');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                set(gca,'xtick',0:5:xlimm(2));
                printstr=sprintf('AngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
                saveas(gcf,[printstr '.fig'],'fig')




                %%Permit ray tracing if desired
                prep_ray_trace;


            elseif strcmp(beamchc,'delaynsum') %Individual snapshots


                Bsum=zeros(length(angles),size(x,1));
                space=head.geom.rd(Igood_el);
                space=space(2)-space(1);
                for Isnap=1:length(angles)
                    if rem(Isnap,10)==0,disp(Isnap);end
                    xtot=delaynsum(x,-angles(Isnap),space,Fs,Igood_el,cc);
                    % Bsum(Isnap,:)=abs(hilbert(xtot'));
                    Bsum(Isnap,:)=xtot';
                end


                %Optional matched filter Ludovic...
                if exist('matched_filter','var')
                    %x_filt=zeros(N,Nel);
                    %h=waitbar(0,'Filtering data ...');
                    %Bsum_org=filter(flipud(matched_filter),1,Bsum')';
                    for K=1:size(Bsum,1)
                        Bsum(K,:)=filter(flipud(matched_filter),1,Bsum(K,:)');
                    end


                end
                Bsum_org=Bsum;
                Bsum=abs(hilbert(Bsum.'))';
                Bsum=Bsum/max(max(Bsum));
                Bsum=20*log10(Bsum);

                %plot vs sin angle
                tt=(1:length(xtot))/Fs;
                figure(20);clf;
                pcolor(tt*1000,sin(angles*pi/180),Bsum);shading flat
                caxis([-20 0]);colorbar
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('sine of Elevation angle');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                if diff(xlimm)<100
                    set(gca,'xtick',0:5:xlimm(2));
                else
                    set(gca,'xtick',0:100:xlimm(2));

                end
                printstr=sprintf('SinAngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
                saveas(gcf,[printstr '.fig'],'fig')

                %plot migration angle...
                figure(21);clf;
                pcolor(tt*1000,angles,Bsum);shading flat;axis('ij')
                caxis([-20 0]);colorbar
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('Elevation angle (deg)');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                if diff(xlimm)<100
                    set(gca,'xtick',0:5:xlimm(2));
                else
                    set(gca,'xtick',0:100:xlimm(2));

                end
                printstr=sprintf('AngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);
                saveas(gcf,[printstr '.fig'],'fig')

                %Repeat using cross-correlation result (to remove
                %   complexities in time waveform)
                figure(21)
                uiwait(msgbox('Please select a box to cross-correlate with other beams (matched filter), hit return to skip:'));
                tmp=ginput(2);
                % tmp=[];
                if isempty(tmp)
                    prep_ray_trace;
                    return
                end

                Bsum_corr=zeros(length(angles),size(Bsum,2));
                angle_want=mean(tmp(:,2));
                indd=round(tmp(1,1)*Fs):round(Fs*tmp(2,1));
                [~,Iang]=min(abs(angle_want-angles));
                x_match=fliplr((Bsum_org(Iang,:)));Nx=length(x_match);

                for Iang=1:length(angles)
                    Bsum_corr(Iang,:)= conv(x_match,Bsum_org(Iang,:),'same');

                end
                Bsum_corr=20*log10(abs(hilbert(Bsum_corr.')))';
                Bsum_corr=Bsum_corr-max(max(Bsum_corr));
                %Bsum_conv=Bsum_conv./max(max(Bsum_conv));
                %plot vs sin angle
                figure(22)
                imagesc(1000*(1:size(Bsum_corr,2))/Fs,(angles),Bsum_corr);
                caxis([-20 0]);colorbar
                set(gca,'fontweight','bold','fontsize',14)
                xlabel('Time (msec)');ylabel('sine of Elevation angle');
                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                titstr=sprintf('%s: Nfft: %i, %6.2f to %6.2f kHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),Nfft,min(frange)/1000,max(frange)/1000);
                title(titstr);grid on;orient landscape
                xlimm=xlim;
                if diff(xlimm)<100
                    set(gca,'xtick',0:5:xlimm(2));
                else
                    set(gca,'xtick',0:100:xlimm(2));

                end
                printstr=sprintf('MatchSinAngleVsTime_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                % print(gcf,'-djpeg',[printstr '.jpg']);
                % saveas(gcf,[printstr '.fig'],'fig')
                %%Permit ray tracing if desired
                prep_ray_trace;

                return
            end

        case 3
            beam_str='MV';

            B=MV_beamforming(app, Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
        case 4
            beam_str='CVnMV';

            B=conventional_beamforming(Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);

            B2=MV_beamforming(app, Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
        case 5
            beam_str='RC';

            R=derive_reflection_coefficient2(app, Ksout.Kstot(Igood_el,Igood_el,:),angles,Ksout.freq,head.geom.rd(Igood_el),cc);
        case 6
            beam_str='MFP';

            prompt1={'Model file..','tilt offset between top and bottom phone (m)','ranges (m)', 'depths (m):','plot intermediate images?', ...
                'Nfft for plotting:','Frequency SNR_cutoff (dB), or array with specific frequencies to process..','Primary Eigenvector only?','receiver depths (m):'};

            dlgTitle1='Parameters for matched field processing...';

            [MFP_replica_file,range_str,default_tilt,default_SNR]=load_MFP_scenario;

            if length(MFP_replica_file)>1
                error('Cannot select multiple environmental models during MFP');
            end

            %head.geom.rd=get_VA_2010_depths_and_offset(2010,handles.myfile,MFP_replica_file{1});
            head.geom.rd=head.geom.rd;

            def1={MFP_replica_file{1}, default_tilt{1},range_str{1},'1:1:55','0','2048',['[' num2str(default_SNR{1}) ']'],'yes',num2str(head.geom.rd)};
            answer=inputdlg(prompt1,dlgTitle1,1,def1);
            model_name=answer{1};
            tilt_offset=str2num(answer{2});
            ranges=str2num(answer{3});
            depths=str2num(answer{4});
            plot_chc=str2num(answer{5});
            Nfft2=str2num(answer{6});
            SNRmin=str2num(answer{7});  %May also be frequency
            head.geom.rd=str2num(answer{9});
            Ksout.Nfft=Nfft;
            Ksout.ovlap=ovlap;

            if ~isempty(Ksout.VV)
                V1=squeeze(Ksout.VV(:,1,:));
                for If=1:length(Ksout.freq)
                    Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';

                end
            else
                Ksout.SNR=Inf*ones(size(Ksout.freq));
                Ksout.Kstot_eig=[];
            end

            srcname=sprintf('CSDM_Nfft%i_%s',Nfft,datestr(tdate_start,30));

            if strcmp(answer{8},'yes')&&~isempty(Ksout.VV)
                save_result=matched_field_processor(app, model_name,tilt_offset,ranges,depths,Ksout.Kstot_eig,Ksout.freq,10*log10(Ksout.SNR),head.geom.rd,SNRmin,srcname);
            else
                save_result=matched_field_processor(app, model_name,tilt_offset,ranges,depths,Ksout.Kstot,Ksout.freq,10*log10(Ksout.SNR),head.geom.rd,SNRmin,srcname);
            end


            if ~isempty(save_result)
                save(srcname,'Ksout','head','tdate_start','tlen');
                write_covmat(app, Ksout.freq,Ksout.Kstot,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
                if ~isempty(Ksout.VV)
                    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
                    write_covmat(app, Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
                end
            end
            return
    end

    if ~correlation_flag
        plot_beamforming_results(true); %ppeak picking
        yes=menu('Beamform?','No','Conventional','MV','Both','Reflection Coefficient Estimation','MFP');
    else
        yes=0;
    end
end  %while yes

if correlation_flag
    return
end

yes=menu('Save CSDM?','Yes','No');
if yes==1

    Ksout.Nfft=Nfft;
    Ksout.ovlap=ovlap;
    V1=squeeze(Ksout.VV(:,1,:));
    for If=1:length(Ksout.freq)
        Ksout.Kstot_eig(:,:,If)=V1(:,If)*V1(:,If)';

    end


    %     srcname=sprintf('CSDM_Nfft%i_%s_flipped',Nfft,datestr(tdate_start,30));
    %     write_covmat(Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));
    %
    figure
    imagesc(Ksout.freq,[],10*log10(abs(Ksout.EE)));
    xlabel('frequency (Hz)')
    ylabel('eigenvalue:')
    title(sprintf('eigenvectors of %s',srcname))
    orient landscape
    for I=1:3
        figure(I+1)
        imagesc(abs(squeeze(Ksout.VV(:,I,:))));

    end
    srcname=sprintf('CSDM_Nfft%i_%s',Nfft,datestr(tdate_start,30));
    save(srcname,'Ksout','head','tdate_start','tlen');
    write_covmat(app, Ksout.freq,Ksout.Kstot,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));

    srcname=sprintf('CSDM_Nfft%i_%s_eigenvector',Nfft,datestr(tdate_start,30));
    write_covmat(app, Ksout.freq,Ksout.Kstot_eig,head.geom.rd,srcname,sprintf('Nfft: %i ovlap %6.2f chann: %s',Nfft,ovlap,mat2str(chann)));

end

%%Inner function for prepping a ray trace and/or invariant trace
    function prep_ray_trace(ray_angles2,dt_data)

        repeat=1;
        while repeat==1
            ButtonName = questdlg('What to do?', 'Migration Analysis','Ray Tracing', 'Invariant Tracing', 'Nada', 'Nada');
            switch ButtonName,
                case 'Nada'
                    return
                case 'Ray Tracing'
                    scenarios=raytrace_MATLAB;
                    [Selection,OK]=listdlg('liststring',scenarios,'SelectionMode','single');
                    if ~OK
                        return
                    end

                    %figure(21);pause(0.5);
                    %%Interpret results depending on plot
                    if strcmp(beamchc,'angvstime')||strcmp(beamchc,'delaynsum')

                        prompt1={'Number of paths to pick?'};
                        def1={'0'};
                        %Other defaults: MURI_Feb2014_20140218T140039,
                        answer=inputdlg(prompt1,'Ray picks',1,def1);
                        Nrays=str2num(answer{1});

                        tmp=ginput(Nrays);


                        dt=max(tmp(:,1))-tmp(:,1);  %Rays that arrive first will have positive numbers
                        [dt_meas,Isort]=sort(dt);  %Arrange by first-arriving rays
                        ang_meas=tmp(Isort,2);
                    elseif strcmp(beamchc,'freqvspwr')

                        ang_meas=(ray_angles2);
                        dt_meas=[dt_data];  %No relative arrival information

                        dt_meas=dt_meas-min(dt_meas);
                        [dt_meas,Isort]=sort(dt_meas);
                        ang_meas=ang_meas(Isort);
                    end

                    [x,dt,dz]=raytrace_MATLAB(ang_meas,dt_meas,scenarios{Selection});
                case 'Invariant Tracing'

                    prompt1={'Ranges (km)','Beta','sound speed at receiver (m/sec)'};

                    dlgTitle1='Parameters for invariant tracing...';

                    def1={'[40:10:60 90:10:120]','-9','1481'};
                    answer=inputdlg(prompt1,dlgTitle1,1,def1);

                    try
                        figure(20);
                        pause(0.5);
                        invariant_range_estimate(eval(answer{1}),eval(answer{2}),eval(answer{3}));
                    catch
                        uiwait(errordlg('Invariant Processing failed!'));
                    end


            end % switch

            repeat=menu('Another analysis?','Yes','No');
        end


    end

%%Inner function for plot_beamforming_results
    function ray_angles=plot_beamforming_results(peak_picking)
        ray_angles=[];
        %Plot beamforming output
        try
            close(1);
        end
        figure(1);clf
        if yes==4||yes_eigen*yes==2
            subplot(2,1,1)
        end

        %Image of beamform output vs look angle and frequency.
        imagesc(Ksout.freq/1000,angles,10*log10(B'));
        colorbar
        cmap=colormap;
        caxis([60 100])
        caxis('auto');
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Frequency (kHz)');ylabel('Angle from horizontal (deg)');grid on;
        title(sprintf('%s, %i FFT, %i elements',datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
        set(gcf,'colormap',cmap(1:4:64,:));

        figure(2);clf
        if yes==4||yes_eigen*yes==2
            subplot(2,1,1)
        end

        %%%Plot summed beampattern
        df=Ksout.freq(2)-Ksout.freq(1);
        Bsum=sum(10*log10((B)))/length(Ksout.freq);
        plot(Bsum,angles,'k');
        set(gca,'fontweight','bold','fontsize',14);axis('ij');
        xlabel('Mean dB Beampower ');ylabel('Angle from horizontal (deg)');grid on;
        title(sprintf('%s, %i FFT, %i elements',datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));

        if ~peak_picking
            return
        end
        %%%%%%%%%%%
        %%%Look for local maxima:
        not_happy=true;
        hplot=[];
        while not_happy
            if ~isempty(hplot)
                set(hplot,'vis','off');
            end
            prompt1={'Automatic peak pick?','Half-beamwidth (deg)', 'threshold SNR'};
            def1={'no','3', '3'};

            answer=inputdlg(prompt1,'Peakpicking parameters',1,def1);
            if isempty(answer)
                return
            end
            try
                peak_pick_chc=strcmpi(answer{1},'yes');
                halfbeamwidth=eval(answer{2});
                thresholddB=eval(answer{3});

            catch
                errdlg('Could not understand your peakpicking parameters');
                return

            end

            if peak_pick_chc
                peaks=peak_picker_Thode(Bsum,angles,halfbeamwidth,[min(angles) max(angles)],thresholddB);
                figure(2)
                hold on
                hplot=plot(peaks{1}.adp.PdB,peaks{1}.adp.F,'go');
                for I=1:length(peaks{1}.adp.PdB)
                    fprintf('Adp Path %i: %6.2f deg beampower %6.2f dB\n',I,peaks{1}.adp.F(I),peaks{1}.adp.PdB(I));
                    %try
                    %fprintf('ISI Path %i: %6.2f deg beampower %6.2f dB\n',I,peaks{1}.isi.F(I),peaks{1}.isi.PdB(I));
                    %end
                end
                [PdB,Isort]=sort(peaks{1}.adp.PdB,'descend');
                ray_angles=peaks{1}.adp.F(Isort);
            else
                prompt1={'Number of picks'};
                def1={'3'};
                try
                    answer=inputdlg(prompt1,'Peakpicking parameters',1,def1);

                    tmp=ginput(eval(answer{1}));
                catch
                    return
                end
                PdB=tmp(:,1);
                ray_angles=tmp(:,2);

                %Find local maxima near each selection
                fprintf('Original angles selected: %s\n',mat2str(ray_angles,3));
                sigma_ang=2;
                for Iray=1:length(ray_angles)
                    Iseg=find(angles>ray_angles(Iray)-sigma_ang&angles<ray_angles(Iray)+sigma_ang);
                    [~,Imax]=max(Bsum(Iseg));
                    ray_angles(Iray)=angles(Iseg(Imax));

                end
                fprintf('Revised angles: %s\n',mat2str(ray_angles,3));

                %[PdB,Isort]=sort(PdB,'descend');
                %ray_angles=ray_angles(Isort);
            end
            not_happy = questdlg('Redo Peak Pick?', ...
                ' Question', ...
                'Yes', 'No', 'No');
            not_happy=strcmp(not_happy,'Yes');
        end %while not_happy
        %%save ray angles, sorted by power

        if length(ray_angles)>1
            [Bshift,dt_data]=conventional_beam_timedelay(Ksout.Kstot(Igood_el,Igood_el,:),ray_angles,Ksout.freq,head.geom.rd(Igood_el),cc,Nfft,Fs);
            prep_ray_trace(ray_angles,dt_data);
        end


        if yes==4
            figure(1)
            subplot(2,1,2);
            imagesc(Ksout.freq,angles,10*log10(abs(B2)'));
            colorbar
            %cmap=colormap;
            caxis([60 100])
            caxis('auto');
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');ylabel('Angle from horizontal (deg)');grid on;
            title('MV processor');
            %set(gcf,'colormap',cmap(1:4:64,:));

            figure(2);
            subplot(2,1,2)
            plot(angles,sum(10*log10((B2)))/length(Ksout.freq),'k');
            set(gca,'fontweight','bold','fontsize',14);
            ylabel('Mean Beampower (dB)');xlabel('Angle from horizontal (deg)');grid on;
            title(sprintf(' %s, %i FFT, %i elements',datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));

        elseif yes_eigen*yes==2
            figure(1)
            subplot(2,1,2)
            imagesc(Ksout.freq,angles,10*log10(B_eig'));
            colorbar
            %cmap=colormap;
            caxis([60 100])
            caxis('auto');
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Frequency (Hz)');ylabel('Angle from horizontal (deg)');grid on;
            title(sprintf('Eigenvector %i only: %s, %i FFT, %i elements',Ieig,datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));
            %set(gcf,'colormap',cmap(1:4:64,:));

            figure(2);
            subplot(2,1,2)
            plot(angles,sum(10*log10((B_eig)))/length(Ksout.freq),'k');
            set(gca,'fontweight','bold','fontsize',14);
            ylabel('Mean dB Beampower ');xlabel('Angle from horizontal (deg)');grid on;
            title(sprintf('Eigenvector %i only: %s, %i FFT, %i elements',Ieig,datestr(tdate_start,'dd-mmm-yyyy HH:MM:SS.FFF'),Nfft,length(head.geom.rd)));

        end

        yes_print=menu('Print?','Yes','No');
        if yes_print==1
            for III=1:3
                figure(III)
                orient tall

                ttt=tdate_start+datenum(0,0,0,0,0,min(ftmp(:,1)));
                printstr=sprintf('Beamforming_%s_%ito%ikHz',datestr(ttt,'yyyymmddTHHMMSS.FFF'),floor(min(frange)/1000),floor(max(frange)/1000));
                print(gcf,'-djpeg','-r300',[printstr '.jpg']);

                print(gcf,'-djpeg','-r300',sprintf('Beamforming%s_%s_%ito%ideg_%4.2fres_%i.jpg', ...
                    beam_str,datestr(ttt,'yyyymmddTHHMMSS.FFF'),min(angles),max(angles),angles(2)-angles(1),III));
            end
        end
    end %function plot_beamforming_results

end