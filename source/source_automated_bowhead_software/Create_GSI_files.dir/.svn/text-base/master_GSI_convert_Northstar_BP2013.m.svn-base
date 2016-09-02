%master_GSI_convert_Northstar.m
%% Top script to run to convert raw DASAR to gsi
chpsin = 3;  %Channels per sample, read in
chpsout = 3;   %channels per sample, read out
buffer_samples = 2*3600*1000*chpsin;  %samples to read into memory at once
%buffer_samples = 400*1000*chpsin;  %samples to read into memory at once

for Isite=1:1
try
%         paths.siteyear = sprintf('NA08');
%         paths.raw0 = sprintf('/Volumes/BP_GSI_Data/%s/Raw%s',paths.siteyear,paths.siteyear);
%         paths.output='/Volumes/BP_GSI_Data';
%         paths.dsumry='/Volumes/BP_GSI_Data/dsumryFiles2008';
        paths.siteyear = sprintf('NA13');
        paths.raw0 = sprintf('/Volumes/BP_2013/%s/Raw%s',paths.siteyear,paths.siteyear);
        paths.output='/Volumes/BP_2013';
        paths.dsumry='/Volumes/BP_2013/NA13/StatusNA13';
        
        cdpath = get_summary_data_path(paths.dsumry,paths.siteyear);
        [did,damean,azone,utmzone]=extract_DASAR_data_summary(cdpath);
        
        for id = 1: length(did)
%             if ~isempty(findstr(did(id).lbl{1},'J'))
            if ~isempty(findstr(did(id).lbl{1},'C'))
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
                        %elseif dn==dnend % Conditional added by khkim, 10 Jan 2011
                       %     mtec = mtvend;
                       %     tba=gen_gsi(did(id),utmzone,paths,mtbc,mtec,chpsin,chpsout,buffer_samples);
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
                
            end
        end
        
    catch
        disp(sprintf('Crash: Site %i failed',Isite));
    end
end
