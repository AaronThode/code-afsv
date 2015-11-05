%function compute_and_plot_group_velocity

for Id=1:size(bath,1)
    c{Id}.phase=real(2*pi*freq'*ones(1,run_options.max_kr)./kr{Id});
    c{Id}.group=(2*pi*(freq(2)-freq(1))./diff(real(kr{Id})));
    c{Id}.group=[c{Id}.group; ones(1,run_options.max_kr)];
    c{Id}.s12=1./c{Id}.group(:,2)-1./c{Id}.group(:,1);
    c{Id}.s13=1./c{Id}.group(:,3)-1./c{Id}.group(:,1);
    c{Id}.s23=1./c{Id}.group(:,3)-1./c{Id}.group(:,2);
end

if size(bath,1)>1
    dr=diff(bath);dr=dr(:,1);
    dr=dr./(bath(end,1)-bath(1,1));
    cgroup=zeros(size(c{1}.group));
    s12=zeros(size(c{1}.s12));
    s13=zeros(size(c{1}.s13));
    s23=zeros(size(c{1}.s23));
    for Id=1:(size(bath,1)-1)
        cgroup=cgroup+dr(Id)./c{Id}.group;
        s12=s12+dr(Id)*c{Id}.s12;
        s13=s13+dr(Id)*c{Id}.s13;
        s23=s23+dr(Id)*c{Id}.s23;
    end
    cgroup=1./cgroup;
else
    cgroup=c{1}.group;
    s12=c{1}.s12;
    s13=c{1}.s13;
    s23=c{1}.s23;
end

Maxmodes=30;
U=zeros(length(rd),length(freq),Maxmodes);
for If=1:length(freq)
    for Im=1:Maxmodes
        try
            U(:,If,Im)=Umode{If}{1}.phi(:,Im);
            %U(:,If,2)=Umode{If}{1}.phi(:,2);
            %U(:,If,3)=Umode{If}{1}.phi(:,3);
        catch
            %disp(sprintf('Freq %6.2f no mode %i',freq(If), Im));
        end
    end
end

figure('name','Group Velocities')
subplot(2,1,1)
plot(freq,cgroup,'o');ylim([1350 1500]);xlim([0 500]);grid on
set(gca,'fontweight','bold','fontsize',14);
xlabel('Frequency (Hz)');ylabel('group velocity (m/s)')
title(sprintf('D: %i, cmedtop: %s svp profile: %s',D,mat2str(cmedtop),case_ssp))
xlimm=xlim;
set(gca,'yminorgrid','on');set(gca,'xminorgrid','on')

subplot(2,3,4);
imagesc(freq,rd,squeeze(U(:,:,1)));grid on;colormap(jet)
set(gca,'fontweight','bold','fontsize',14);
xlabel('Frequency (Hz)');ylabel('Depth (m)');title('Mode 1 at receiver');
set(gca,'yminorgrid','on');%set(gca,'xminorgrid','on');


subplot(2,3,5);
imagesc(freq,rd,squeeze(U(:,:,2)));grid on;colormap(jet)
set(gca,'fontweight','bold','fontsize',14);
xlabel('Frequency (Hz)');ylabel('Depth (m)');title('Mode 2 at receiver');
set(gca,'yminorgrid','on');%set(gca,'xminorgrid','on');

subplot(2,3,6);
imagesc(freq,rd,squeeze(U(:,:,3)));grid on;colormap(jet)
set(gca,'fontweight','bold','fontsize',14);
xlabel('Frequency (Hz)');ylabel('Depth (m)');title('Mode 3 at receiver');
set(gca,'yminorgrid','on');%set(gca,'xminorgrid','on');


orient tall
print('-djpeg',[savename '_GroupVelocity.jpg']);

%%%%%%%%%%%%%%%%%%%%
%  Plot difference in group velocities 1 and 2, 2 and 3; useful for estimating range from spectrogram
figure('name','Differences in group slowness');
plot(freq,1./(1000*s12),'k',freq,1./(1000*s23),'g',freq,1./(1000*s13),'r');
set(gca,'fontweight','bold','fontsize',14); ylim([0 500])
xlabel('Frequency (Hz)');ylabel ('cXY (km/sec)');grid on
title(sprintf('R=[c12,c13]*dt: %i, cmedtop: %s svp profile: %s',D,mat2str(cmedtop),case_ssp))
xlim(xlimm);
set(gca,'yminorgrid','on');set(gca,'xminorgrid','on');
legend('Mode 1-2','Mode 2-3','mode 1-3')

orient landscape
print('-djpeg',[savename '_GroupSlowness.jpg']);


%             figure
%             [junk,Iwant]=min(abs(freq-85));
%             plot(rd,squeeze(U(:,Iwant,:)),'-o')
%             [junk,Iwant]=min(abs(freq-105));hold on;
%             plot(rd,squeeze(U(:,Iwant,:)),'--o');grid on

