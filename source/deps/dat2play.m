%%	Read in boat noise sample and convert
%	Input location
dir_name	=	'.';%'/Users/jit/Jonah/Data/SEASWAP/SEASWAP_2012/PLAYBACK_SAMPLES/EngineHaulingSounds_ForAaron';
file_name	=	'11181941.ADI';%'07061144.ADI';
[pathstr, fname, ext]	=	fileparts(file_name);

[x,tfs,tfe,fs]	=	read_adi_file(dir_name,file_name,[],0,3600);



%%	Covert file to uint16 limits
Nmax	=	(2^16)-1;
fs_new	=	100e3;%48e3;		%Actual playback rate

if fs ~= fs_new
    [q,p]	=	rat(fs/fs_new);
    x		=	resample(x,p,q);
end

%x	=	x - min(x);		x	=	x * Nmax /max(x);	%	uint binary output
x	=	x - mean(x);	x	=	x / max(abs(x));	%	+/- 1 float output


%%	Write out new file
% fid		=	fopen(['PLAY.' file_name '.bin'],'w','ieee-be');
% COUNT	=	fwrite(fid,x(round(end*0.9):end),'uint16');
% fclose(fid);
% disp(['Bytes written = ' num2str(COUNT)]);

audiowrite([fname '.wav'], x, fs_new);