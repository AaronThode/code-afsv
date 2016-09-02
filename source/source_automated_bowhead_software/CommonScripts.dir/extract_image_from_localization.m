%function [y,start_ctime,t,head]=extract_image_from_locations(locations,Istation,Icall,station,filename,params,nchan)
%% Give a station and index number, extract templated image
%% for signal processing
function Image=extract_image_from_localization(locations,Istation,Icall,station,filename,params,nchan,Idebug);  %Assumes no buffertime
%(locations,K,J,

if ~exist('nchan')
    nchan=1;
end
start_ctime=locations{Icall}.ctime_min(Istation);

if start_ctime==0
    Image=[];
    return
end
%tlen=mean(locations{Icall}.Totalduration(Istation));
%params=station(Istation).param;
Iindex=locations{Icall}.station_indicies(Istation);
[y,start_ctime]=extract_signal_from_station(station(Istation),Iindex,filename,'debug',nchan);
for I=1:length(nchan),
    [SS,FF,TT,PP]=spectrogram(y(:,nchan(I)),hanning(params.Nfft),round(params.ovlap*params.Nfft),params.Nfft,params.Fs,'yaxis');

    Igood=find(TT>=params.morph.eq_time);
    PP=10*log10(abs(PP(:,Igood)));
    SS=SS(:,Igood);
    TT=TT(Igood);

    tmp=station(Istation).Image{Iindex};
    %     if islogical(tmp)
    %         keyboard;
    %     end
    Image1=sign(double(tmp));
    mincol=min([size(Image1,2) size(SS,2)]);

    Image{I}=Image1(:,1:mincol).*flipud(SS(:,1:mincol));

    if Idebug>0
        figure;
        subplot(2,1,1)
        imagesc(TT,FF,PP);colorbar
        subplot(2,1,2)
        imagesc(10*log10(abs(Image{I})));colorbar
        title(['Channel ' int2str(I)]);
        pause
        close
    end
    %Image=[];

end

end
