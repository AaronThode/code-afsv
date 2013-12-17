
% function write_Java_script(pchc,param)
% Melania Guerra and Aaron Thode
%
% Function will write a Cshell script for calculating SEL and PSD using the Java
% PreProcessor
% Oct 2007-first version
% July 24, 2009, extensively revised for multielement files.
% Input:
%       pchc: 'SEL' or 'PSD'
%       param: structure that contains the following fields.
%           The number of PSD or SEL is determined by length of sec_avg and f_low, respectively
%           file_dir,  string containing directory holding files of
%               interest.
%           dir_out-desired output directory
%           exten:  string containing file extension of form '*.MDAT' or
%               '*.gsi' or '*0828*.gsi' etc d
%           Nfft,
%           ovlap: number of samples overlapping between snapshots
%           Fs: sampling rate of file in Hz.  May be overwritten based on
%               particular file header.
%           channel:   If 1, use first channel, etc. Note this is MATLAB
%                       indexing convention--this script substracts one
%                       from this value so that channel 1 becomes 0, etc.
%                       If channel<0, then import beamforming parameters.
%
%                       If channel is a vector, add a for loop in script to generate
%                       an output file for every channel in the vector...
%           dumpsize: number of samples per channel to read into RAM memory
%               at a time.
%           nstart:  Samples per channel to skip before processing.  If 0 process from
%               file start.
%           nsamples: Samples per channel to process.  I 0 process entire
%               file.
%           f_high: for SEL, upper frequency to sum over.  Also used by PSD averaging
%           f_low   for SEL, lower frequency to sum over
%               (f_low and f_figh and sev_avg may be vectors, defining the number of SEL  or PSD objects
%                to be created).
%           sec_avg:  for PSD, seconds over which to average
%           filter_type: for PSD, 'default' or 'autocorrelation';
%           maxT:  maximum frequency or output time for PSD, depending on 'filter_type' chosen
%           (beamforming)-optional for beamforming.  If nstart<0
%                   beamforming.angles: vector of angles to beamform.  A
%                       value of zero is broadside, +-90 degrees is endfire.
%                   beamforming.badel:  channels not to use for
%                       beamforming.  Uses MATLAB convention (1=first
%                       element, etc.)
%                   beamforming.Lspace:  For MDAT file, spacing in m
%                       between elements.
%
%    
%
%
% Output:
%       file masterSEL.scr contains a command line with the following...
%        java -jar <executable_path>/PreProcessV1.jar $data_filename $dirout $Npsd $Nsel $Ndet  ...
%           $channel $nstart $npt $dumpsize $Nfft $Fs $ovlap',param.main_dir);


function write_Java_script(pchc,param)
Nrun=[];
if strcmp(lower(pchc),'sel')&&isfield(param,'f_low')
    Nrun=length(param.f_low);
    if length(param.f_low)~=length(param.f_high);
        error('f_low and f_high must be same lengths');
    end
elseif strcmp(lower(pchc),'psd')&&isfield(param,'sec_avg')
    Nrun=length(param.sec_avg);
    
end
if isempty(Nrun)
    error('fields in param do not match your choice of function');
end

if ~isfield(param,'program_dir')
    param.program_dir=load_program_location;
end


try
    f_sel = fopen(sprintf('master%s.scr',pchc),'wt');
    fprintf(f_sel, '#!/bin/csh \n\n\n');
    
    fprintf(f_sel,'set dirout = "%s"    \n',param.dir_out);
    if strcmp(lower(pchc),'sel')
        fprintf(f_sel,'set Nsel = %i        \n',Nrun);
        command_str=sprintf('java -jar %s $i $dirout 0 $Nsel 0 ',param.program_dir);
        
    elseif strcmp(lower(pchc),'psd')
        fprintf(f_sel,'set Npsd = %i        \n',Nrun);
        command_str=sprintf('java -jar %s $i $dirout $Npsd 0 0 ',param.program_dir);
        
    end
    
    %%%%%%Beamforming additions
    if (isfield(param,'beamforming'))
        nang=length(param.beamforming.angles);
        nbad=length(param.beamforming.badel);
        if isfield(param,'f_high')
            maxfreq=max(param.f_high);
        else
            maxfreq=param.Fs/2;
        end
        if isfield(param,'f_low')
            minfreq=min(param.f_low);
        else
            minfreq=0;
        end
        Ntotal=2+nang+1+nbad+1;
        if isfield(param.beamforming,'Lspace')
            Lspace=param.beamforming.Lspace;
            Ntotal=Ntotal+1;
        else
            Lspace=[];
        end
        
        
        fprintf(f_sel,'set array_params = %i \n',-Ntotal);
        command_str=sprintf('%s $array_params %i ',command_str,nang);
        for I=1:nang
            command_str=sprintf('%s %i ',command_str,param.beamforming.angles(I));
        end
        command_str=sprintf('%s %i ',command_str,nbad);
        for I=1:nbad
            command_str=sprintf('%s %i ',command_str,param.beamforming.badel(I)-1);
        end
        command_str=sprintf('%s %6.2f %6.2f %6.2f ',command_str,minfreq,maxfreq,param.beamforming.Lspace);
        
        
        
    elseif ~isfield(param,'channel')
        fprintf(f_sel,'set channel = 0\n');
        command_str=sprintf('%s $channel ',command_str);
        
    elseif length(param.channel)==1&&param.channel>=0  %just access a single channel
        fprintf(f_sel,'set channel = %i     \n',param.channel-1);
        command_str=sprintf('%s $channel ',command_str);
    elseif length(param.channel)>1
        %% don''t define anything here.
        command_str=sprintf('%s $channel ',command_str);
    end
    
    
    
    if ~isfield(param,'nstart')
        fprintf(f_sel,'set nstart = 0\n');
    else
        fprintf(f_sel,'set nstart = %11.0f     \n',param.nstart);
    end
    
    if ~isfield(param,'nsamples')
        fprintf(f_sel,'set npt = 0\n');
    else
        fprintf(f_sel,'set npt = %11.0f       \n',param.nsamples);
    end
    
    
    if ~isfield(param,'dumpsize')
        fprintf(f_sel,'set dumpsize = 1000000\n');
    else
        fprintf(f_sel,'set dumpsize = %11.0f    \n',param.dumpsize);
    end
    
    
    fprintf(f_sel,'set Nfft = %i        \n',param.Nfft);
    
    if ~isfield(param,'Fs')
        fprintf(f_sel,'set Fs = 2048\n');
    else
        fprintf(f_sel,'set Fs = %i    \n',param.Fs);
    end
    
    fprintf(f_sel,'set ovlap = %i        \n',param.ovlap);
    
    
    command_str = sprintf('%s $nstart $npt $dumpsize $Nfft $Fs $ovlap',command_str);
    
    if ~isfield(param,'filter_type')
        param.filter_type='default';
    end
    
    
    if Nrun>0
        for I = 1:Nrun
            
            if strcmp(lower(pchc),'sel')
                command_str = sprintf('%s $flo%i $fhi%i',command_str,I,I);
                fprintf(f_sel,'set flo%i =  %i          \n',I,param.f_low(I));
                fprintf(f_sel,'set fhi%i = %i          \n\n',I,param.f_high(I));
            elseif strcmp(lower(pchc),'psd')
                
                if isfield(param.filter_type,'autocorrelation')&&I==1
                    if ~isfield(param,'f_low')|~isfield(param,'f_high')
                        error('Need to define f_low and/or f_high for autocorrelation option');
                    end
                    if ~isfield(param,'maxT')
                        param.maxT=Nfft/2;
                    end
                    fprintf(f_sel,'set sec_avg%i =  %i          \n',I,-param.sec_avg(I));  %Put a minus sign to indicate further PSD parameters...
                    command_str = sprintf('%s $sec_avg%i',command_str,I);
                    command_str = sprintf('%s %10.6f %s %i %i',command_str,param.maxT,'autocorrelation',param.f_low(1),param.f_high(1));
                end
                
                
                %%%%%What if regular PSD wants to define an upper frequency limit?
                if isfield(param,'f_high')&&~isfield(param.filter_type,'autocorrelation')
                    fprintf(f_sel,'set sec_avg%i =  %i          \n',I,param.sec_avg(I));
                    fprintf(f_sel,'set fhi%i =  %i          \n',I,param.f_high(I));
                    command_str = sprintf('%s $sec_avg%i $fhi%i',command_str,I,I);
                    
                end
                
                
            end  %if PSD
        end
        
        
        
    end  %if Nrun>0


    
    
    
    if isfield(param,'channel')&&length(param.channel)>1&&~isfield(param,'beamforming')
        temp=mat2str(param.channel);
        temp=temp(2:(end-1));
        fprintf(f_sel,'foreach channel (%s)\n',temp);
        fprintf(f_sel,'foreach i (%s/%s)       \n',param.file_dir,param.exten);
        fprintf(f_sel,'echo " Starting processing of " $i, channel $channel\n');
    else
        fprintf(f_sel,'foreach i (%s/%s)       \n',param.file_dir,param.exten);
        fprintf(f_sel,'echo " Starting processing of " $i\n');
        
    end
    
    command_str = [command_str ' >> out.txt\n'];
    fprintf(f_sel,command_str);
    fprintf(f_sel,'echo " Finished processing of " $i\n');
    fprintf(f_sel,'end\n');
    
    if isfield(param,'channel')&&length(param.channel)>1&&~isfield(param,'beamforming')
        fprintf(f_sel,'end\n');
    end
    fclose(f_sel);
    if strcmp(lower(pchc),'sel')
        !chmod +x mastersel.scr
    else
        !chmod +x masterpsd.scr
    end
    
catch
    disp('Could not write file!');
end
    
    
    
