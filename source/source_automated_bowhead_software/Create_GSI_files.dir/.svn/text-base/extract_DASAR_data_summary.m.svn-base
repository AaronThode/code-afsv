function [did,damean,azone,utmzone]=extract_DASAR_data_summary(cdpath)

[numeric,txt,raw] = xlsread(cdpath,'Main','a1:x8','basic');
[idmax ncols] = size(numeric);
mdatebase = datenum('30-dec-1899');
% get DASAR labels
for id = 1:idmax
    did(id).lbl = raw(id+1,1);
    unit(id) = numeric(id,1);
    did(id).dser = raw(id+1,3);
    dv1 = datevec(numeric(id,3)+mdatebase);
    seconds =86400*numeric(id,4);
    dv1(4)= floor(seconds/3600);
    seconds = seconds -3600*dv1(4);
    dv1(5) = floor(seconds/60);
    dv1(6) = seconds-60*dv1(5);
    did(id).ctstart = mat2c_tm(dv1);      %A/D start time
    did(id).ctinst=mat2c_tm(datevec(numeric(id,5)+mdatebase)); %Installed
    did(id).ctrecov=mat2c_tm(datevec(numeric(id,6)+mdatebase));    %recovered time
        % get data start and end times if present
    s= numeric(id,7);
    if s>0
        did(id).ctdstart=mat2c_tm(datevec(s+mdatebase));              %data start time
        did(id).ctdend=mat2c_tm(datevec(numeric(id,8)+mdatebase));  %data end time
    else
        did(id).ctdstart=did(id).ctinst;
        did(id).ctdend=did(id).ctinst;
    end
    did(id).sectors=numeric(id,9);   %sectors written
    did(id).magbrg=numeric(id,10);   %mag bearing
    did(id).utmx=numeric(id,11);    %UTM position
    utmx(id)=numeric(id,11);
    did(id).utmy=numeric(id,12);
    utmy(id)=numeric(id,12);
    did(id).ddepth=numeric(id,13);   %Dasar depth, meters
    did(id).chps = numeric(id,14);  %channels per sample for DASAR
    did(id).use = raw(id+1,16);     % use only if not 'F'
    did(id).brefa=numeric(id,16);    %Dasar reference bearing 
    did(id).xgf=numeric(id,17);      %sinch gain correction facyor
    s = numeric(id,18);
    if s > 0
        did(id).ctref=mat2c_tm(datevec(s + mdatebase)); %ref for time correction
    else
        did(id).ctref=ctstart(id);
    end
    did(id).tint = numeric(id,19);    %zero intercept
    did(id).tdrift=numeric(id,20);  %clock drift sec/day
    did(id).brefb=numeric(id,21);   %ref bearing after change
    did(id).stda=numeric(id,22);    % std dev of bearing cal before change
    did(id).stdb=numeric(id,23);    % std after change
end

damean=[mean(utmx(1:idmax)) mean(utmy(1:idmax))];
% recover UTM zone from column label
azone = char(raw(1,12));
utmzone = str2num(azone(length(azone)));
disp([num2str(idmax),' Dasars to be Processed'])