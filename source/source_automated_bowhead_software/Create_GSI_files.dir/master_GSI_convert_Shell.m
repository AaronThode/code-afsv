% rdsumry_c.m
% script to prepare environment for gensiof Processing
% called from compiled version of supergensiof.
% 1.reads Operator ID
% 2.Reads Site desired to process
% 3.Opens and reads applicable dsumry file to get all Dasar-related information.
%
%=======================================
clear all
% Request operator ID
%OperatorID=input('Enter Operator ID:','s');
%tlastfn=[OperatorID 'TL.txt'];
%
% Have operator select site for processing

chpsin = 3;  %Channels per sample, read in
chpsout = 3;   %channels per sample, read out
buffer_samples = 2*3600*1000*chpsin;  %samples to read into memory at once
%buffer_samples = 400*1000*chpsin;  %samples to read into memory at once


for Isite=1:5
    try
        paths.siteyear = sprintf('S%i10',Isite);
        paths.raw0 = sprintf('/Volumes/Shell_2010/%s/Raw%s',paths.siteyear,paths.siteyear);
        paths.output='/Volumes/Shell2010_GSI_Data_copy';
        paths.dsumry='/Volumes/Shell_2010/dsumryFiles2010';
        
        cdpath = get_summary_data_path(paths.dsumry,paths.siteyear);
        [did,damean,azone,utmzone]=extract_DASAR_data_summary(cdpath);
        
        for id = 1: length(did)
            if ~strcmp(did(id).use,'F')
                paths.DasarID = char(did(id).lbl);
                disp(sprintf('Processing DASAR %s',paths.DasarID));
                paths.raw=sprintf('%s/%s',paths.raw0,paths.DasarID);
                paths.sitefolder=sprintf('%s/%sgsif',paths.output,paths.siteyear);
                if ~exist(paths.sitefolder,'dir')
                    status1 = mkdir(paths.sitefolder);
                end
                paths.dasarfolder=sprintf('%s/%s',paths.sitefolder,paths.DasarID);
                
                if ~exist(paths.dasarfolder,'dir')
                    status2 = mkdir(paths.dasarfolder);
                end
                
                mtvstart = c2mat_tm(did(id).ctdstart);
                mtvend   = c2mat_tm(did(id).ctdend);
                
                
                totime=0;
                dnstart = datenum(mtvstart);
                dnend   = datenum(mtvend);
                
                for dn = floor(dnstart):floor(dnend)  %For each day...
                    tic
                    dngo = max(dn,dnstart);
                    mtbc = datevec(dngo);
                    disp(sprintf('Processing day %s of DASAR %s',datestr(mtbc),paths.DasarID));
                    mtec = datevec(dn +1);
                    try
                        if dn ~= dnend
                            tba=gen_gsi(did(id),utmzone,paths,mtbc,mtec,chpsin,chpsout,buffer_samples);
                        end
                    catch
                        disp(sprintf('Crash! Processing day %s of DASAR %s failed',datestr(mtbc),paths.DasarID));
                    end
                    totime=totime+toc;
                    toc
                end
                % clear DasarID
                %beep
                %donesound
                
            end  %if
        end %for
        
    catch
        disp(sprintf('Crash: Site %i failed',Isite));
    end %try
end %Isite
