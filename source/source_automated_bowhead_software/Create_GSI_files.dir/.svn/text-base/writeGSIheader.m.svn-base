function writeGSIheader(ofid,did,ctbc,ctec,utmzone,tsamp,ndc)
%use,tdrift,utmx,utmy,ddepth,brefa,
%ndc: number of data channels.
%
% There is a 1-sector (512 bytes) header containing the following descriptive data:
%
% Bytes	Format		Contents
% 10	char		DASAR designation in form SSyydr
% 1	byte		Number of data channels (3)
% 1	char		Ch1 contents (P)
% 1 	char		Ch2 contents (X)
% 1	char		Ch3 contents(Y)
% 1	char		Ch4 contents if present(Z)
% 1	char		Proc Status of DASAR
% 			=T use for all
% 			=F do not use
% 			=D use for detection; not location
% 48	byte		Reserved
% 8	Double		ctbc, c-time of first sample(absolute)
% 8	Double		ctec, c-time of sample after last(absolute)
% 8	Double		ctdrift, clock drift in sec/day
% 8	Double		samprn, Nominal sample rate (1,000 samples/sec)
% 8	Double		UTM X,m of DASAR location
% 8	Double		UTM Y,m
% 8	Double		Depth, m
% 8	Double		UTM Zone
% 8	Double		Brefa, bearing in degrees of DASAR Cosine Axis
% 376	byte		Reserved

%count1 = fwrite(ofid,char(strcat(did.lbl, tseries )));
ncl = length(char(did.lbl));
blstr = '          ';
%ndc = 3;
label = [char(did.lbl) blstr(1:10-ncl) char(ndc) 'PXY ' char(did.use)];
count1 = fwrite(ofid,label);    %bytes
fill(1:48)=0;
prec = char('double');
count2 = fwrite(ofid,fill(1:48));   %bytes
count3 = fwrite(ofid,[ctbc ctec did.tdrift tsamp did.utmx did.utmy did.ddepth utmzone did.brefa ], 'double');  %doubles=8bytes
fill(1:376)=0;
count4 = fwrite(ofid,fill(1:376));  %bytes
end
