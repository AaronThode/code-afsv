%%%%%%%%%yout=display_automated_crosslink_data(local,Icall,param,filename,yes_overlay)
% If yes_overlay==1, assume that spectrograms already up, overplot
%   local:  The 'locations' variable in the automated processing results
%   Icall: index of a particular call in 'local'
%   param:  structure with the following fields:  Nfft, ovlap, Fs
%   filename: cell array of strings, each string containing absolute path and filename of raw data    
function [yout,tabs]=display_automated_crosslink_data(local,Icall,param,filename,yes_overlay)
%close all;
yout=[];
%strr='abcdefghijklmnop';
strr=[];
for II=1:length(filename)
    Islash=max(findstr(filename{II},'/'));
    %disp(goodFile{1}{II}(Islash+6-1));
    strr=[strr filename{II}(Islash+6-1)];
end

extra_time=2;
Ianchor=find(local{Icall}.dt==0);
%Inot=find(isnan(local{Icall}.dt));
Nstations=length(local{Icall}.dt);
if ~exist('yes_overlay')
    yes_overlay=zeros(1,Nstations);
end

Igood=find(local{Icall}.Nsegments>0);
dt=local{Icall}.dt(Igood);

duration=local{Icall}.Totalduration(Igood);
flo=[local{Icall}.Totalfmin(Igood)];
fhi=[local{Icall}.Totalfmax(Igood)];
SNR=local{Icall}.SNR(Igood);
ctime=local{Icall}.ctime_min(Igood);
bearing=local{Icall}.bearing(Igood);
tabs=datenum(1970,1,1,0,0,ctime);
cplot_start=local{Icall}.ctime_debug(Igood);
%tlen=local{Icall}.duration_debug(Igood)+param.energy.bufferTime;


disp(sprintf('Automated localized call %i at mean time of %s',Icall,datestr(mean(tabs))));
disp(sprintf('%i stations used: %s', length(Igood),strr(Igood)));

disp(sprintf('Relative arrival times at %s is %s',strr(Igood),mat2str(dt,3)));
for I=1:length(Igood),
    disp(sprintf('Station %s: Time: %s bearing: %6.2f Total Flo = %6.2f, Total Fhi = %6.2f, duration=%6.2f, SNR=%6.2f \n',strr(Igood(I)),datestr(tabs(I)),bearing(I),flo(I),fhi(I),duration(I),SNR(I)));
end

linewidth=2;
Nrows=4;
Nfig=1;
Iplot=0;
if all(yes_overlay==0)
    figure
    set(gcf,'units','normalized','pos', [0.05 0.05 0.45 0.9]);
end
%for I=1:Nstations
for I=Nstations:-1:1
    Iplot=Iplot+1;
    if Iplot==2*Nrows+1,
        Iplot=1;Nfig=Nfig+1;
        figure;
        set(gcf,'units','normalized','pos', [0.05+(Nfig-1)*0.45 0.05 0.45 0.9]);
    end
    
    subplot(Nrows,2,Iplot);
    Iwant=(find(I==Igood));
    
    try
        if all(yes_overlay==0)
            if ~isempty(Iwant),
                y=extract_signal_from_localization(local,I, Icall,filename{I},'debug',1,extra_time);
            else
                y=readfile(filename{I},mean(cplot_start)-10,20,1,'ctime','calibrate');
                
            end
            yout{I}=y;
            
            [SS,FF,TT,PP]=spectrogram(y,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,param.Fs,'yaxis');
            imagesc(TT,FF,10*log10(abs(PP)));axis('xy');
            set(gca,'fontweight','bold','fontsize',14);xlabel('Time (sec)');ylabel('Hz');
            if ~isempty(Iwant)
                title(sprintf('Station %s index: %i Date: %s',upper(strr(I)),local{Icall}.station_indicies(Igood(Iwant)),datestr(tabs(Iwant))));
            else
                title(sprintf('Station %s',upper(strr(I))));
            end
        end
        
        if ~isempty(Iwant)
            hold on
            cstartt=ctime(Iwant)-cplot_start(Iwant)+extra_time;
            if yes_overlay(I)>0
                offset=cplot_start(Iwant)-yes_overlay(I);
            else
                offset=0;
            end
            hh=line(cstartt + offset+[0 duration(Iwant)],[flo(Iwant) flo(Iwant)]);set(hh,'Color',[0 0 0],'linewidth',linewidth);
            hh=line(cstartt + offset+[0 duration(Iwant)],[fhi(Iwant) fhi(Iwant)]);set(hh,'Color',[0 0 0],'linewidth',linewidth);
            hh=line(cstartt + offset+[0 0],[flo(Iwant) fhi(Iwant)]);set(hh,'Color',[0 0 0],'linewidth',1);
            hh=line(cstartt + offset+duration(Iwant)*[1 1],[flo(Iwant) fhi(Iwant)]);set(hh,'Color',[0 0 0],'linewidth',linewidth);
            %title(sprintf('Call %i, station %s, start time %s,   duration %6.2f, SNR: %6.2f, ',Icall,strr(Igood(Iwant)), ...
            %    ctime2str(cplot_start(Iwant)),duration(Iwant), SNR(Iwant)));
        else
            %title(sprintf('Call %i, station %s, start time %s, ',Icall,strr(I),ctime2str(mean(cplot_start))));
            
        end
        %caxis([-20 60]);
        %caxis([40 110]);
        caxis([40 120])
        
        colorbar;
    catch
        disp(sprintf('Could''nt find filename'));
    end
    
end

end