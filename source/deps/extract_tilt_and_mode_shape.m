
%depth_shift=-2;  %How much to shift my estimated depth by
%freq_want=[84 104];
%freq_want=[120:5:145];
prompt={'Frequencies at which to extract modes (Hz):'};
def={'[85:100]'};
dlgTitle	=	sprintf('Frequencies to extract modes');
lineNo		=	ones(size(prompt));
answer		=	inputdlg(prompt,dlgTitle,lineNo,def);

if	isempty(answer)
    return;
end

freq_want=eval(answer{1});
Nfft_filt=2^nextpow2(MM(1));
Fwarp=(0:Nfft_filt-1)*Fs/Nfft_filt;
Igood=find(Fwarp>=(minfreq)&(Fwarp<=(maxfreq)));

Uwarp=fft(mode_stack,Nfft_filt);  %Spectrum of each mode arrival
Uwarp=Uwarp(Igood,:,:);

%Sum energy in each modal arrival. Assign sign later.
Umode_broadband=squeeze(sqrt((sum(abs(mode_stack).^2,1))))';

if MM(2)==1
    Umode_broadband=Umode_broadband';
end

Nc=length(Nchan);
Nf=length(Igood);
%Ifreq_plot=[1 Nf/4 Nf/2 3*Nf/4 Nf];

%Option: pick frequency with most power.
%Fpower=squeeze(sum(sum(shiftdim(abs(Uwarp).^2,1))));
%Ifreq_plot=[1 Nf/2  Nf];
%Ifreq_plot=round(([ 84 98 105]-minfreq)*Nfft_filt/Fs);
Ifreq_plot=round((freq_want-minfreq)*Nfft_filt/Fs);

freq_plot=Fwarp(Igood(Ifreq_plot));
Nplot=length(freq_plot);

%Cyle through each mode and extract mode shape from various frequencies,
% as well as broadband time series
legendCell=cellstr(num2str(freq_plot', '%2.1f Hz'))';

for Im=1:MM(2)
    
    %normalize both individual-frequency based and broadband modes.
    Umode=squeeze(Uwarp(:,Im,:)).';
    Umode_broadband(:,Im)=Umode_broadband(:,Im)/norm(Umode_broadband(:,Im));

    for II=1:Nf
        Umode(:,II)=Umode(:,II)/norm(Umode(:,II));
    end
    
    %Extract phase along array
    Uang_uncorrected=angle(Umode./(ones(Nc,1)*Umode(1,:)));
    
    if Im==1
        Umode1=Umode;  %We assume that mode 1 phase arises from tilt
    end
    
    tmp=(Umode.*conj(Umode1));  %Remove tilt or phase offsets
    Uang=angle(tmp./(ones(Nc,1)*tmp(1,:)));  %Reference phase to first element
    
    %Determine sign of mode: observe when corrected phase
    %   (relative to bottom mode) exceeds pi/2
    sgn_chnge=-(abs(Uang(:,Ifreq_plot))>pi/2);
    sgn_chnge=2*sgn_chnge+1;  %Convert 0 and -1 to 1 and -1;
    if Nplot>1
        sgn_chnge=median(sgn_chnge')';
    end
    
    Umode_broadband(:,Im)=sgn_chnge.*Umode_broadband(:,Im);
    Umode_shape=(sgn_chnge*ones(1,Nplot)).*abs(Umode(:,Ifreq_plot));
    
    %end
    %tmp=sum(Umode.^2);
    %[junk,Imax]=max(tmp);
    
    %save('vars','Umode','hdr','Ifreq_plot')
    
    figure(20);
    subplot(MM(2),3,3*Im-2);
    plot(abs(Umode(:,Ifreq_plot)),hdr.geom.rd,'o-');grid;axis('ij')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('%s: abs value of mode %i',datestr(tdate_start),Im));

    if Im==1
        ylabel('depth (m)');
        legend(legendCell,'location','westoutside');   
    end
    
    subplot(MM(2),3,3*Im-1);
    plot((Uang_uncorrected(:,Ifreq_plot)),hdr.geom.rd,'o-');grid;axis('ij')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('uncorrected phase of mode %i',Im));
    if Im==1
        ylabel('depth (m)');
    end
    
    subplot(MM(2),3,3*Im);
    plot((Uang(:,Ifreq_plot)),hdr.geom.rd,'o-');grid;axis('ij')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('tilt-corrected phase of mode %i',Im));
    xlabel('Phase (rad)');
    %     plot(Umode_fin(:,Ifreq_plot),hdr.geom.rd,'o-');grid;axis('ij')
    %     title(sprintf('Recovered mode %i',Im));
    if Im==1
        ylabel('depth (m)');
    end   

    %% To delete
     try
        load C:/Users/cedri_000/Desktop/Kraken/Sims/Separated/mode_t_sim mode_t_sim
        freqs1=fieldnames(mode_t_sim);
        tmp1=zeros(55,Nplot);
        for ff=1:Nplot
            [~,minf]=min(abs(str2double(strrep(strrep(freqs1,'f',''),'_','.'))-freq_plot(ff)));
            tmp1(:,ff)=mode_t_sim.(char(freqs1(minf)))(:,Im);
        end
        sim=real(tmp1);
     catch 
        warning('Simulation does not exist!: %s','C:/Users/cedri_000/Desktop/Kraken/Sims/Separated/mode_t_sim ' );
        sim=0.1;
     end
    
    col='rgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbkyrgbky';
    figure(19)
    subplot(1,MM(2),Im)
    ylim([0 55])
    hold on,grid on
    for ff=1:Nplot
        plot(Umode_shape(:,ff),hdr.geom.rd,'Color',col(ff)),
        try
            plot(sim(:,ff)*max(Umode_shape(:,ff))/max(sim(:,ff)),1:55,'--','Color',col(ff))
        catch
            warning('Simulation does not exist!: %s','C:/Users/cedri_000/Desktop/Kraken/Sims/Separated/mode_t_sim ' );
       
        end
    end
%     plot(Umode_shape(:,:),hdr.geom.rd),
%     plot(sim,1:55,'--')
    plot(Umode_broadband(:,Im),hdr.geom.rd,'ko')
    set(gca,'fontweight','bold','fontsize',14);
    title(sprintf('Normalized mode %i',Im));
    ylabel('depth(m)');
    title(sprintf('Inverted mode %i',Im));
    axis('ij')
    if Im==1; legend([legendCell,{'broadband','simulation'}],'location','westoutside'); end
end

%%Plot recovered tilt
U1_phase=angle(Umode1./(ones(Nc,1)*Umode1(1,:)));
U1_tilt=unwrap(U1_phase)./(ones(Nc,1)*2*pi*Fwarp(Igood)/c1);

%Simple broadband calculations
xt1=squeeze(mode_stack(:,1,:));
RR=40;
maxlag=round(50*RR*Fs/c1);

for I=1:(Nc-1)
    y1=resample(xt1(:,I),RR,1);
    y2=resample(xt1(:,I+1),RR,1);
    [cc,lags]=xcov(y1,y2,maxlag);
    [junk,Imax]=max(abs(cc));
    tilt_broadband(I+1)=c1*lags(Imax)/(RR*Fs);
    
end
tilt_broadband=cumsum(tilt_broadband);
tilt_broadband=tilt_broadband-tilt_broadband(1);  %Make first phone the origin

figure(18)
subplot(1,2,1)
plot(U1_phase(:,Ifreq_plot),hdr.geom.rd);axis('ij');grid
xlabel('rad');ylabel('depth (m)');
title('Array tilt in terms of mode 1 phase');
if Im==1; legend(legendCell,'location','westoutside');end

subplot(1,2,2);
plot(U1_tilt(:,Ifreq_plot),hdr.geom.rd);axis('ij');grid
hold on
plot(tilt_broadband,hdr.geom.rd,'ko')
title('Array tilt in terms of mode 1 offset');
xlabel('m offset');ylabel('depth (m)');
grid on
legend([legendCell,{'broadband'}],'location','best')

for ff=1:Nplot
    freq=freq_plot(ff);
    mode_extract.(strrep(num2str(freq,'f%1.1f'),'.','_'))=squeeze(Umode_shape(:,ff,:));
end

%save mode_extract

for II=18:20
    figure(II)
    saveas(gcf,sprintf('ModeTiltWarpInversion_%i.jpg',II))
end