
function obuf=ReadRAWDataWriteGsi(tba,nsample,chpsin,chpsout,datapath,yearsite,DasarId,fstr)
%tba-byte offset desired from start of raw data files
%nsample: number of samples to read in

obuf = zeros(chpsin,buffer_samples,'uint16'); %Total output buffer size in adcrs
%

iraw=floor(tba/512000000);
baddr=tba-512000000*iraw;
fstr=int2str(1000+iraw);
fn = char(strcat(datapath,'Raw',yearsite,'\',DasarID,'\', DasarID,fstr(3:4),'.raw'));
%    disp(['baddr=' sprintf('%12.0d',baddr)])
% number of samples to read
%K=nsample*chpis;
%K = nsample;  %samples
%read time series data into array
[fid,message] = fopen(fn,'r','ieee-be') ;
status = fseek(fid,baddr,-1);
dataok(id)=status~=-1;
if(status~=0);
    ferror(fid)
end
%if file opened and seek was ok read data.
if status ~= -1
    % data is big-endian (Motorola)
    [x,count] = fread(fid,[chpsin,nsample],'uint16'); % 8 = 4 ch x 2 bytes/sample
    %count
    fclose(fid);
    
    obuf(:,1:size(x,2))=x;
    
    %if record is short, read rest from start of next raw file
    if count < chpsin*nsample
        iraw=iraw+1;
        fstr=int2str(1000+iraw);
        fn = char(strcat(datapath,'RAW',yearsite,'\',DasarID,'\',DasarID,fstr(3:4),'.raw'));
        [fid,message] = fopen(fn,'r','ieee-be');
        status=fseek(fid,0,-1);
        [obuf(:,size(x,2)+(1:nsample)),count2] = fread(fid,[chpsin,nsample-size(x,2)],'uint16');
        fclose(fid);
        
    end
    
    
end
end
