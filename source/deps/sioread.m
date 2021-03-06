function [x, Header] = sioread(filename,p1,npi,channels)
%-----------------------------------------------------------------------
% sioread.m
%
% This program runs under windows, unix, and macs.  
%
% function x=sioread(filename,p1,npi,channels);
%
% Inputs:
% 	filename: Name of sio file to read
% 	p1:	Point to start reading ( 0 < p1 < np)
% 	npi: 	Number of points to read in (if 0, read in all points)
% 	channels: Single number or vector containing the channels to read
% 		(example-to read channels 4 thru 10 enter 4:10)
%       First channel is channel 1.
%
% Example:  xdata = sioread('../data.dir/test1.sio',10,100,[1 2 4:6]);
% xdata will be a matlab array of 100 points and 5 channels.  The
% points start at the siofile point #10 and the channels extracted are
% channels 1,2,4,5,6
%
%
% originally written by Aaron Thode.  Modified by Geoff Edelmann and 
% James Murray. (Final Version: July, 2001)
%
% read bug fixed by Dave Ensberg (July 2003)
%
% Modified 2014-02, by Jit Sarkar
%	Added output for header info, and option to output just header
%
%-----------------------------------------------------------------------
 
% pathname=[pwd '/' filename];

% since the majority of sio files will be created as 'SUN' sio
% files, we will assume they are 'big endian' and check for
% compliance.
%
endian='b';
fid=fopen(filename,'r',endian);
fseek(fid,28,-1);
bs=fread(fid,1,'long');  % should be 32677
if bs ~= 32677
  fclose(fid);
  endian='l';
  fid=fopen(filename,'r',endian);
  fseek(fid,28,-1);
  bs=fread(fid,1,'long'); % should be 32677
  if bs ~= 32677
    error(['Problem with byte swap constant:' num2str(bs)])
  end
end

% I can just use fseek(fid,0,-1) to position to beginning of file
% but I think that closing and reopening is cleaner.
fclose(fid);
fid=fopen(filename,'r',endian);
 

 id=fread(fid,1,'long');  % ID Number
 nr=fread(fid,1,'long');  % Number of Records in File
 rl=fread(fid,1,'long');  % Record Length in Bytes
 nc=fread(fid,1,'long');  % Number of Channels in File
 sl=fread(fid,1,'long');  % Bytes/ Point
 type = 'float';
 if sl==2,
        type='short';
 end
 f0=fread(fid,1,'long');  % 0 = integer, 1 = real
 np=fread(fid,1,'long');  % Number of Points per Channel
 bs=fread(fid,1,'long');  % should be 32677
 fn=fread(fid,24,'char'); % Read in the file name
 com= fread(fid,72,'char'); % Read in the Comment String
 year=fread(fid,1,'long');
 day=fread(fid,1,'long');
 hour=fread(fid,1,'long');
 min=fread(fid,1,'long');
 temp=fread(fid,10,'float');
 Fs=temp(4);
 AD=temp(6);
 
 rechan=ceil(nr/nc);     %records per channel
 ptrec=rl/sl;            %points/record
 
 %	Header object, for output - Jit (2014-02-17)
 Header.ID			=	id;
 Header.N_Rec		=	nr;
 Header.BperRec		=	rl;
 Header.N_Chan		=	nc;
 Header.BperPoint	=	sl;
 Header.IsReal		=	f0;
 Header.PperChan	=	np;
 Header.RperChan	=	rechan;
 Header.PperRec		=	ptrec;
 Header.fname		=	char(fn(:).');
 Header.comment		=	char(com(:).');
 Header.date=datenum(year,0,day,hour,min,0);
 Header.Fs=Fs;
 
%  so for a random MFNA data file:    RAVA6.14049032400.000.sio  
% 
% usr(  1)= 0.282222E-41      usi(  1)=         2014       year
% usr(  2)= 0.686636E-43      usi(  2)=           49         Jday
% usr(  3)= 0.420390E-44      usi(  3)=            3          Hour
% usr(  4)= 0.336312E-43      usi(  4)=           24         Min
% usr(  8)=  25000.0          usi(  8)=   1187205120       sampling frequency
% usr(  9)= 0.224208E-43      usi(  9)=           16         A/D number of bits
% usr( 10)=  2.50000          usi( 10)=   1075838976     A/D peak voltage range ( 2x this is full scale )

 %	If no channels or points requested, just return header
 if nargin==1|| isempty(npi) || isempty(channels) || (npi==0) || (max(channels)==0)
	x	=	[];
	return;
 end
	 
 
 

% check to make sure the channel list is valid
for ii = 1:length(channels)
  if (channels(ii) <= 0) 
    error(['Channel must be positive: ' channels(ii)])
  end
  if (channels(ii) > nc) 
    error(['Channel does not exist:' channels(ii)])
  end
end

%-- First calculate what time period is desired

 r1=floor((p1-1)/ptrec)+1;		%record where we wish to start recording
 if npi <= 0			%if npi is 0, then we want all points		
   npi = np-p1+1;
 end
 p2=p1+npi-1;			%ending position in points for each channel.
 r2=ceil(p2/ptrec);		%record where we end recording

 if p2>np,			%-- Error checking to prevent reading past EOF
       disp('Warning:  p1 + npi > np');
       p2=np;r2=rechan;
 end
 totalrec=r2-r1+1;		%number of records desired read.

% we're going to read in entire records into the variable tmpx, then
% we'll worry about the start point.

 pp1 = mod(p1-1,ptrec)+1;		% p1's position in record

 x = zeros(p2-p1+1,length(channels));	% "allocate" a matrix for final result
 tmpx = zeros(totalrec*ptrec,1);	% make a temporary matrix
% Loop over the desired channels
 for J=1:length(channels)
    count = 1;				% Start 
    trec = (r1-1)*nc + channels(J);	% First Record to read for this channel
    for R = 1:totalrec
      status = fseek(fid,rl*trec,-1);        % position to the desired record
      if status == -1
        error(ferror(fid))
      end
      tmpx(count:count+ptrec-1) = fread(fid,ptrec,type);	% Read in a record's worth of points
      count = count+ptrec;	% adjust for the next set of points
      trec = trec + nc;		% Next record for this channel is nc records away
    end
    x(1:p2-p1+1,J) = tmpx(pp1:pp1+p2-p1);
 end


 fclose(fid);

%-- end of program
