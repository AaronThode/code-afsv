function R = derive_reflection_coefficient2(app, Kstot, angles, freq, rd, c)
%function [B,wout]=conventional_beamforming(Ks,angles,freq,Lz,c,yesnorm)

tilt=input('Enter tilt estimate in deg (0):');
if isempty(tilt)
    tilt=0;
end

% angles=angles0+tilt;
%
% Igood=find(angles<=90&angles>=-90);
% angles=angles(Igood);
%
% Nup=Iang-1;
% Ndown=length(Iang:length(angles));
% if Ndown>=Nup
%    angles=angles(1:(Iang+Nup));
% else
%     angles=angles((Nup-Ndown+2):end);
%end
Lz=rd-rd(1);

B=conventional_beamforming(Kstot,angles,freq,Lz,c,'yesnorm');
%B=MV_beamforming(Kstot,angles,freq,Lz,c,'yesnorm');
B=real(B).';  % B now [angles freq]
figure
imagesc(freq,angles,10*log10(abs(B)))
caxis([-20 0]);
axis('xy')
set(gca,'fontweight','bold','fontsize',14);
xlabel('Frequency (Hz)');
ylabel('Angle (positive is towards surface');
grid on
%title(sprintf('Theoretical conventional beamforming response, even spacing, source angle %6.2f',angle))
colorbar('fontweight','bold','fontsize',14)


for It=1:length(tilt)
    [junk,Iang]=min(abs(angles-tilt(It)));

    %R is bottom/surface, or negative over positive
    %upp=flipud(B(1:(Iang-1),:));  %Negative angles point towards bottom
    %downn=B((Iang+1):end,:);
    %Nmin=min([size(downn,1) size(upp,1)]);
    %R=downn(1:Nmin,:)./upp(1:Nmin,:);

    downn=flipud(B(1:(Iang-1),:));  %Negative angles point towards bottom
    upp=B((Iang+1):end,:); %Positive angles point towards surface
    Nmin=min([size(downn,1) size(upp,1)]);
    R=downn(1:Nmin,:)./upp(1:Nmin,:);

    %%Note that R is the power reflection coefficient (Harrison and Simmons, Eq. (2)), so 10*logR please..
    R=-10*log10(abs(R))';
    angle_graze=angles(Iang+(1:Nmin))-tilt(It);
    figure

    subplot(3,1,1)
    imagesc(10*log10(abs(upp')))
    subplot(3,1,2)
    imagesc(10*log10(abs(downn')))
    subplot(3,1,3)
    imagesc(angle_graze,freq,R);
    caxis([0 15]);
    set(gca,'fontweight','bold','fontsize',14);
    grid on;colorbar;axis('xy')
    xlabel('Grazing angle (deg)');
    ylabel('Frequency (Hz)')
    title(sprintf('reflection loss (dB), tilt %6.2f',tilt(It)))

end

yes=input('Save?');
if ~isempty(yes)
    save(sprintf('Refl_%s.mat',num2str(mean(tilt))),'angle_graze','freq','R','tilt','rd','c','Kstot','angles','Lz');
end
% disp('Select a frequency:');
% tmp=ginput(1);
% [junk,Iwant]=min(abs(tmp(2)-freq));
% figure
% plot(angle_graze,R(Iwant,:));
% keyboard
end
