% readEnergySnips.m
%function [x,nstarts,npts,feq,eq,Ireturn,head]=readEnergySnips(fn,index,wordlen,stack_chc,keep_file_open)
%
% matlab script to read *.snip files, or Energy detection time series files
%
% Input:
%   fn is filename string,
%   index is a vector of indicies to read in data.  The values do not have to be contiguous.
%   wordlen is the primitive data type stored in the file,  must either be 'short' or 'double'.
%   stack_chc is either 'cat' or 'cell'.  If the former, the time series are concatenated into
%       a single vector, otherwise arranged in a cell matrix.
%   keep_file_open: if exists, keep file open for a persistant read (saves
%       time if a large file and you plan to return to collect more snips)
% Output:
%    x: output data
%    nstarts: sample at which each time series starts in fn
%    npts: number of samples in each snippet
%    feq:  Frequencies associated with eq in Hz--the center frequency of band
%       used to compute mean PSD
%    eq: mean PSD across detector bandwidth used to flag this snippet.
%       Combined with frequency information in the summary file, can be used
%       to estimate equalized PSD that was used to normalize the data.  For
%       example, the SNR of this snippet can be compared by dividing eq
%       into  |fft(x,Nfft)|^2/(Fs*Nfft);
%    Ireturn: current index in the file..
%    head: header file information--the input parameters for the detector..
% 

function [x,nstarts,npts,feq,eq,Ireturn,head]=readEnergySnips(fn,index,wordlen,stack_chc,keep_file_open)

persistent fid Icurrent

if isempty(index)
    disp('Error: your indicie request is empty');
    nstarts=[];npts=[];Ireturn=[];x=[];
end
if nargin<4,
    stack_chc='cat';
end
if strcmp(wordlen,'short')||strcmp(wordlen,'uint16')
    nbyte=2;
    prec='uint16';
elseif strcmp(wordlen,'double')
    nbyte=8;
    prec='double';
else
    error('readEnergySnips: wordlen not recognized');
end
if strcmp(stack_chc,'cat')
    x=[];eq=[];
elseif strcmp(stack_chc,'cell')
    x{1}=[];eq{1}=[];

end
if ~isinf(index)
    nstarts=zeros(1,length(index));
    npts=nstarts;
end

head=readEnergy_header(fn);
feq=0.5*(head.flow+head.fhigh).';
if ~isempty(fid),
    try
        ftell(fid);
    catch
        disp('File pointer no longer valid, resetting...');
        fclose('all');
        fid=[];
        Icurrent=[];
    end
end

if isempty(fid)|isempty(Icurrent)
    fclose('all');
    %keyboard;
    fid = fopen(char(fn),'r','ieee-be');
    fseek(fid,512,-1);
    % Check if current block is desired
    Icurrent=1;
    disp('Opening file');
else
    %disp(sprintf('Current byte position upon file entry: %i number of calls into file: %i',ftell(fid),Icurrent));
end
Imark=1;

while (Icurrent<=max(index))
    %if I==2145,
    %   keyboard;
    %end
    try
        ctime=fread(fid,1,'double');

        %Read equalization
        Npsd=fread(fid,1,'int16');
        eq0=fread(fid,Npsd,'double');
        checkk=fread(fid,1,'int16');
        vec=fread(fid,2,'int64');
        nstarts(Imark)=vec(1);
        npts(Imark)=vec(2);
        %if npt<0
        %    keyboard;
        %end
    catch
        disp(sprintf('End of file reached at sample %i',Icurrent));
        Ireturn=Icurrent;
        fclose(fid);fid=[];
        Icurrent=[];  %Important, need to reset file

        return
    end %try
    %keyboard;
    if (isinf(index)|Icurrent==index(Imark))
        try
            x0=fread(fid,[head.nchan npts(Imark)],wordlen);
            if strcmp(stack_chc,'cat');
                x=[x x0];
                eq=[eq eq0];
            else
                x{Imark}=x0;
                eq{Imark}=eq0;
            end
            %disp(sprintf('Call succesfully read: %i th in file',Imark));

            Imark=Imark+1;

        catch
            %disp(sprintf('requested index %i too high',index(Imark)));
            disp('Could not read time segment...');
            fclose(fid);fid=[];
            Ireturn=Icurrent;
            return
        end

    else
        fseek(fid,nbyte*head.nchan*npts(Imark),0);
    end %if I=index(Imark);
    %fseek(fid,8,0); %Skip '0' marker
    zerock=fread(fid,1,'double');
    if (zerock~=0),
        disp('zero check fails');
        keyboard;
    end
    Icurrent=Icurrent+1;

end


Ireturn=Icurrent;

if ~exist('keep_file_open')
    disp('Closing file...');
    fclose(fid);
    fid=[];Icurrent=[];
end

if isempty(x{1}),
    % keyboard;
end
%N = 128;
%spectrogram(omi,hanning(N),N/2,N,head.Fs,'yaxis') % gram
%end

if 1==0,
    Nfft=128;ovlap=0.75;Fs=1000;
    figure
    for I=1:length(data),
        subplot(length(x),1,I);
        [S,F,T,P]=spectrogram(x{I},hanning(Nfft),round(ovlap*Nfft),Nfft,Fs,'yaxis');
        imagesc(T,F,10*log10(abs(P)));axis('xy')
        title(nstarts(I)/Fs);
    end

end

end
