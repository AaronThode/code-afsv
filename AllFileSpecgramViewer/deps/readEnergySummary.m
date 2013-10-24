% readEnergySummary.m
% function [data,head,Igood]=readEnergySummary(fn, index, interval_suppress);
% matlab script to read *.detsum files, or Energy detection summary files
%
% Inputs:
% fn is filename string,
% index is a vector of indicies to read in data...
% interval_suppress is a two-element vector containing a time range in
%   seconds.  If the time lapse between the previous and current detection
%   lies within this interval, do not include in final result...
%Outputs:
%   NOTE: subband outpts have been disables
%  data: structure of detections that contain the following fields:
%        .nstart:  The start sample of when the energy first exceeded
%        threshold at any detector
%
%        .npt:  The number of samples during which the energy exceeds the
%        threshold for any subdetector
%
%        .magnitude: Peak SEL reached by the detection.  Organized as above.
%
%        .ctime:  C-time of point where energy exceeded threshold, using
%        corrected sample rate (i.e. clock drift incorporated).
%
%        .interval:  time in samples between this detection and previous
%        detection.
%
%        .names: cell array of strings describing features stored in
%               "features"
%        {'min_freq','min_duration','max_freq','max_duration','peak_freq',
%               'peak_duration','peak','ctime_peak','total_duration'};

%
%        .features: [Nfeatures x length(index)] matrix of detection
%           features.  Note that 'peak' is in dB units.  All other units
%           in Hz or sample length.  Note that 'total_duration' is the same
%           as data.npt.
%   head:  Header file of detection function input parameters, including Nfft and
%           Fs.
%   Igood: indicies of "index" vector trimmed to those units that pass the interval
%   check.
%

function [data,head,Igood]=readEnergySummary(fn, index, interval_suppress)

Igood=[];
head=readEnergy_header(fn);
fmid=0.5*(head.flow+head.fhigh);
fid = fopen(char(fn),'r','ieee-be');
fseek(fid,0,1);
flen=ftell(fid);
fseek(fid,512,-1);

% In a summary block there are three longs, three ints (fmax, etc)
% three longs (nstarts), one double (amp_max), one long (peak_index), and one double ('0')
summary_bytes=3*8+3*4+3*8+8+8+8;
%summary_bytes=8*3+head.Ndetectors*8*3+8;
nsamples=ceil(flen/summary_bytes);
fprintf('There are %i samples in %s\n',nsamples,fn);
%In a summary block there is two longs, one double, one int and one long
%per detector, 

%data.nstart=zeros(1+head.Ndetectors,length(index));
%data.npt=data.nstart;
%data.magnitude=zeros(head.Ndetectors,length(index));
if ~isinf(index)
    data.ctime=zeros(1,length(index));
else
    data.ctime=zeros(1,nsamples);
end
data.npt=data.ctime;
data.nstart=data.ctime;
data.names={'min_freq','min_duration','max_freq','max_duration','peak_freq','peak_duration','peak','ctime_peak','total_duration'};
data.features=zeros(length(data.names),length(index));

% Check if current block is desired
I=1;
Imark=1;


while (I<=max(index))
    if rem(I,10000)==0
        fprintf('%6.2f percent done\n',100*ftell(fid)/flen);
    end
    if (isinf(index)|I==index(Imark))
        try
            error_flag_type=0;
            vec=fread(fid,2,'int64');
            data.nstart(1,Imark)=vec(1,1); %nstart_total
            data.npt(1,Imark)=vec(2,1);
            data.ctime(Imark)=fread(fid,1,'double');
           % keyboard;
            %Read suband data
            error_flag_type=1;
            for J=1:3,
                data.features(2*J-1,Imark)=fmid(1+fread(fid,1,'int32'));
                data.features(2*J,Imark)=fread(fid,1,'int64');
            end
            error_flag_type=2;
            data.features(2*J+1,Imark)=10*log10(fread(fid,1,'double'));
            data.features(2*J+2,Imark)=fread(fid,1,'double');
            data.features(2*J+3,Imark)=data.npt(1,Imark)/head.Fs;
            %for J=1:head.Ndetectors,
            %    data.nstart(1+J,Imark)=fread(fid,1,'int64');
            %    data.npt(1+J,Imark)=fread(fid,1,'int64');
            %    data.magnitude(J,Imark)=fread(fid,1,'double');
            %end
            fseek(fid,8,0);
            Imark=Imark+1;
        catch
            disp(sprintf('End of file reached after detection %i, error flag type %i',I,error_flag_type));
            break;

        end
    else
        fseek(fid,summary_bytes,0);
    end
    I=I+1;

end

%%Trim data

data.nstart=trim(data.nstart,Imark);
data.npt=trim(data.npt,Imark);
data.features=data.features(:,1:(Imark-1));
data.ctime=trim(data.ctime,Imark);
if length(data.ctime)>0,
    data.interval=[0 diff(data.nstart(1,:))];
else
    data.interval=[];
end
fclose(fid);

if length(data.nstart)<1,
    Igood=[];
    disp('No calls detected in file');
    return;
end

%%Remove repetitive signals if desired...
if nargin>2,

    ngood=round(head.Fs*interval_suppress);
    Igood=find(data.interval<ngood(1)|data.interval>ngood(2));

    %Now get really clever, and pick up repetitive calls
    %  that were interrupted by an original call...
    Ifirst=Igood(1);
    Igood=Igood(2:end);

    two_call_interval=data.interval(Igood)+data.interval(Igood-1);
    Igood2=find(two_call_interval<ngood(1)|two_call_interval>ngood(2));
    Igood=[Ifirst Igood(Igood2) ];

    data.nstart=data.nstart(Igood);
    data.npt=data.npt(Igood);
    %data.magnitude=data.magnitude(:,Igood);
    data.ctime= data.ctime(Igood);
    data.interval=data.interval(Igood);
    if ~isinf(index),
        Igood=index(Igood);
    end
    disp(sprintf('Removing repetitive signals, %i left',length(Igood)));

else
    Igood=1:length(data.ctime);
end
end

function x=trim(x,I)
x=x(:,1:(I-1));

%N = 128;
%spectrogram(omi,hanning(N),N/2,N,head.Fs,'yaxis') % gram
%end


end
