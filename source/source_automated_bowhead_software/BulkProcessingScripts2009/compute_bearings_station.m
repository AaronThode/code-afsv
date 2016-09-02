%  function station=compute_bearings_station(station,goodFile,param,algchc,Idebug,Nsamples)
%  add bearing, bearing_rel, kappa, and bearing_sd (standard deviation) to station object.
%  For example, bearing of 100th detection on station 5 is station(5).bearing(100);

function station=compute_bearings_station(station,goodFile,param,algchc,Idebug,Nsamples)

if isempty(station.ctime_min)
    station.bearing=[];
    station.kappa=[];
    station.bearing_sd=[];
    
    station.bearing_rel=[];
    return
end
if strcmp(algchc,'sel_ratio')
    buffertime=param.morph.eq_time;
elseif strcmp(algchc,'sel_ratio_FFT')
    buffertime=param.morph.eq_time;
else
    buffertime=0;
end


if ~isempty(findstr('T2007',goodFile))
    cal07flag=1;
    brefa_table=calibrate_bearing_Shell2007(param.calibration_dir,goodFile);
else
    cal07flag=0;
end

for J=1:length(station.ctime_min),
    if rem(J,50)==0,disp(sprintf('%5.2f percent done',100*J/length(station.ctime_min)));end
    station.bearing(J)=0;
    station.kappa(J)=0;
    station.bearing_sd(J)=0;
    station.bearing_rel(J)=0;
    
    ctime=station.ctime_min(J)-buffertime;
    if ctime<=0,
        station.bearing(J)=NaN;
        station.bearing_rel(J)=NaN;
        continue
    end
    
    tlen=station.Totalduration(J)+2*buffertime;
    [y,t,head]=readfile(goodFile,ctime,tlen,1:3,'ctime','calibrate');
    
    
    
    if isempty(y)
        station.bearing(J)=NaN;
        station.bearing_rel(J)=NaN;
        continue
    end
    
    [thet0,kappa,sd]=extract_bearings(y,buffertime,param.Nfft,param.Fs,station.Totalfmin(J),station.Totalfmax(J),algchc,Nsamples);
    
    
    
    %Convert measured bearing into true bearing...
    if cal07flag==0
        thet=bnorm(thet0+head.brefa);
    else
        [junk,Icol]=calibrate_bearing_Shell2007(param.calibration_dir,goodFile,1);
        thet= bnorm(interp1(0:360,brefa_table(:,Icol),bnorm(thet0)));
    end
    if Idebug>0,
        info_str=sprintf('%s, DASAR %i from alg_chc %s relative bearing: %6.2f true bearing: %6.2f +/- %6.2f', ...
            ctime2str(station.ctime_min(J)), J,algchc,thet0,thet,sd);
        disp(info_str);
        % plot_spectrogram(full(y{1}))
        % disp(' ');
        
    end
    
    station.bearing(J)=thet;
    station.bearing_rel(J)=thet0;
    if ~isempty(kappa)
        station.kappa(J)=kappa;
        station.bearing_sd(J)=sd;
    end
end %J

end









