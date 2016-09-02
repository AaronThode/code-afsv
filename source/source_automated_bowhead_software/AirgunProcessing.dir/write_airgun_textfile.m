%%%%%%%%%write_airgun_textfile.m%%%%%%%%%%%
% function write_airgun_textfile(fname,data)
% Definition of columns:
%	   A -1 for any value means the program was unable to compute a value...
%       (1)  Date Time: date and time in local AK time
%       (2)  interval (s): estimated interval in s between airgun shots.
%       (3)  peak pressure (dB re 1uPa).  The peak pressure attained by the pulse (noise not subtracted).
%       (4)  bandlimited signal rms without noise subtraction(dB re 1uPa). Evaluated between 10 and 450 Hz.
%       (5)  bandlimited signal SEL without noise subtraction (dB re 1uPa^2-s).Evaluated between 10 and 450 Hz.
%       (6)  bandlimited signal rms with noise subtraction(dB re 1uPa).Evaluated between 10 and 450 Hz.
%       (7)  bandlimited signal SEL with noise subtraction(dB re 1uPa^2-s). Evaluated between 10 and 450 Hz.
%       (8)  signal duration (s): time between 5 and 95% cumulative energy of pulse, bandlimited between 10 and 450 Hz,
%               when noise estimate is subtracted.
%       (9)  noise rms (dB re 1uPa):  Estimated rms level of noise sample before pulse onset.
%       (10) noise duration (s): duration of noise used to estimate (9).  NOT used to estimate (11)
%               Very small values indicate noise sample before pulse was highly nonstationary.
%       (11) noise SEL (dB re 1uPa^2-s): Quantity equal to column(8)*[column(9)].^2, when columns are converted
%               to linear units.  The intention of this column is to compare with column (5) or (7) to get SNR.
%       (12) peak frequency (Hz): Frequency between 10 and 450 Hz with greatest power spectral density of
%               pulse.
%       (13) number clips (formerly "clipped?"):  The number of RAW samples that lie within 3 samples of 65536 or within 3 samples of
%               zero.  The RAW samples are data samples before equalization, filtering, and scaling.  I'd recommend removing associated values from consideration if values are greater than 5.
%       (14) bearing: The azimuth of the pulse in degrees from true north.
%
% Some relationships when converted to linear units:
%       (5)=(8)*(4).^2;  %Signal SEL without noise equals (signal duration) * square (signal rms without noise)
%       (7)=(8)*(6).^2;  %Signal SEL without noise equals (signal duration) * square (rms without noise)
%       (4)=(7)+(8)*(9).^2; %Relationship between SEL without noise subtracted and SEL with noise subtracted.
%       (4)=(7) + (11) ;  % How the "noise SEL relates" to the other SEL quantities.
%   This expression is true in dB units.
%     SNR (dB) = (7)-(11)=(6)-(9); % %       
% fname: string with name of output file
% data: structure of data produced by airgun_processing.  In particular
%   should have the following fields:
%   .ctime
%   .ICI
%   .bearing
%   .clipped
% Either
%   .level.rms_Malme
%   .level.SEL_Malme
% or
%   .level.rms_FFT_band
%   .level.SEL_FFT_band
%
%   .level.peak
%   .level.t_Malme  %duration of signal computed by Malme method...
%
%   .noise.rms /broadband
%   .noise.SEL_FFT /broadband
%
%   .bearing: azimuthal arrival of airgun energy
%

function write_airgun_textfile(fname,data)

bandlimited_flag=1;  %If 1, output the bandlimited SEL and rms, otherwise output broadband time domain (SEL_Malme).

bearing_flag=0;


% head_str='Date Time, interval(s), peak pressure (dB re 1uPa), ';
% head_str=[head_str 'bandlimited signal rms without noise subtraction(dB re 1uPa), bandlimited signal SEL without noise subtraction(dB re 1uPa^2-s),'];
% head_str=[head_str 'bandlimited signal rms with noise subtraction(dB re 1uPa), bandlimited signal SEL with noise subtraction(dB re 1uPa^2-s),'];
% head_str=[head_str 'signal duration (s), noise rms (dB re 1uPa), noise duration (s), noise SEL (dB re 1uPa^2-s), peak frequency (Hz), number clips'];

try
   
    %Nsamples=length(data.ctime);
    fid=fopen([fname '.txt'],'w');
   
    if isempty(data)
        fclose(fid);
        return
    end
    
    if isempty(data.ICI)
        fclose(fid);
        return
    end
    
    %Convert to dB:
    level.peak=20*log10(data.level.peak);
        
    if bandlimited_flag==0
        level.rms_Sonly=20*log10(data.level.rms_Malme);
        level.SEL_Sonly=10*log10(data.level.SEL_Malme);
        
        level.SEL_SplusN=data.level.SEL_Malme+(data.noise.rms.^2).*data.level.t_Malme;
        level.rms_SplusN=10*log10(level.SEL_SplusN./data.level.t_Malme);  %20log(sqrt)=10log10
        level.SEL_SplusN=10*log10(level.SEL_SplusN);
        
        noise.rms_Nonly=20*log10(data.noise.rms);
        noise.SEL=10*log10(data.noise.SEL);
        
        head_str='Date Time, airgun interval(s),   peak pressure (dB re 1uPa), ';
        head_str=[head_str '  signal rms S+N(dB re 1uPa),   signal SEL S+N(dB re 1uPa^2-s),'];
        head_str=[head_str '  signal rms S only(dB re 1uPa),   signal SEL S only(dB re 1uPa^2-s),'];
        head_str=[head_str '  signal duration (s),   noise rms (dB re 1uPa),   noise duration (s),   noise SEL (dB re 1uPa^2-s), peak frequency (Hz), number clips'];
        
    else
        level.rms_Sonly=20*log10(data.level.rms_FFT_band(1,:));
        level.SEL_Sonly=10*log10(data.level.SEL_FFT_band(1,:));
        
        noise.rms=20*log10(data.noise.rms_FFT_band(1,:));
        noise.SEL=10*log10(data.noise.SEL_FFT_band(1,:));
        
        level.SEL_SplusN=data.level.SEL_FFT_band(1,:)+(data.noise.rms_FFT_band(1,:).^2).*data.level.t_Malme;
        level.rms_SplusN=10*log10(level.SEL_SplusN./data.level.t_Malme);  %20log(sqrt)=10log10
        level.SEL_SplusN=10*log10(level.SEL_SplusN);
        %airgun_shots.level.SEL_FFT_band(1:length(features.SEL_FFT_band),Icall)
        
        head_str=sprintf('Bandlimited data:%6.2f-%6.2f Hz: Date Time, airgun interval(s),  peak pressure (dB re 1uPa), ',data.freq_bandwidth(1),data.freq_bandwidth(2));
        head_str=[head_str ' signal rms S+N(dB re 1uPa),  signal SEL S+N(dB re 1uPa^2-s),'];
        head_str=[head_str ' signal rms S only(dB re 1uPa),  signal SEL S only(dB re 1uPa^2-s),'];
        head_str=[head_str 'signal duration (s), noise rms (dB re 1uPa), noise duration (s), noise SEL (dB re 1uPa^2-s), peak frequency (Hz), number clips'];
        
    end
    
    if isfield(data,'bearing')
        head_str=[head_str ' ,bearing(deg)'];
        bearing_flag=1;
    end
    fprintf(fid,[head_str '\n']);
    

    %%Remove level.SEL less than zero
    Igood=find(imag(level.SEL_SplusN)==0);
    %save temp_useme
    for I=Igood
        try
            fprintf(fid,'%s,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%6.3f,%10.6f,%6.3f,%10.6f,%10.6f,%i', ...
                ctime2str(data.ctime(I)),data.ICI(I), level.peak(I), ...
                level.rms_SplusN(I),level.SEL_SplusN(I), ...
                level.rms_Sonly(I), level.SEL_Sonly(I), ...
                data.level.t_Malme(I),noise.rms(I),data.noise.duration(I),noise.SEL(I), data.level.peakF(I),data.clipped(I));
            if bearing_flag==1
                fprintf(fid,', %10.6f \n',data.bearing(I));
            else
                fprintf(fid,'\n');
            end
        catch
            disp(sprintf('Sample %i failed to write.\n',I));
        end
    end
    
    
    fclose(fid);
    
    
    %%%Delphine, you can comment out the following section if octave levels
    %%%     are not neeeded...
    fid=fopen([fname '_OctaveLevels.txt'],'w');
    fprintf(fid,'1/3 octave rms levels (dB re 1uPa) at center frequencies specified by second line.\n');
    fprintf(fid,'%i,',data.f_center');
    fprintf(fid,'\n');
    
    
    Nsamples=length(data.ctime);
    
    %Convert to dB:
    
    rms=20*log10(data.level.rms_FFT_third_octave);
    for I=1:Nsamples,
        try
            fprintf(fid,'%s,',ctime2str(data.ctime(I)));
            fprintf(fid,'%10.6f,', rms(:,I));
            fprintf(fid,'\n');
        catch
            disp(sprintf('Sample %i of octave write failed.\n',I));
        end
    end
    fclose(fid);
    
    %%End of comment out section
catch
    
    disp(sprintf('Crash! Failure to write %s',fname));
end

