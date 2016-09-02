function [goodDASAR,goodFile,goodNames,DASAR_str_out]=find_DASAR_dates_Processed(date_str,Isite,DASAR_str,outputdir,Icase)
%function [goodDASAR,goodFile,goodNames,DASAR_str]=find_DASAR_dates_Processed(date_str,Isite,DASAR_str,rawdatadir,Icase)
% Given a site and date, produce strings of DASAR names and filenames that
%   exist at that site and during that date.
%   Input:
%       date_str: string of form '20070902'
%       Isite: integer of Site number
%       DASAR_str: string of letters of desired DASARS within a site, e.g.:
%           'abc'.  If '*' all available DASARS are desired.  DASARS must end
%           in a letter 'a-z'.
%       rawdatadir: String giving complete path to raw data top level
%           (folder that contains 'Site_01', etc.)
%       Icase: descriptive scenario string.  Must contain a substring
%           'ShellYY' or BPYY', where YY is a two-digit year.
%   Output:
%       goodDASAR: cell array of available DASAR names
%       goodFile: full path name and filename of each file that fits criteria
%       goodName: short name of file, minus the pathname...
%       DASAR_str: string of letters representing DASARS actually present
%           in data...

%%What DASARS are available at the site?
%keyboard;

goodDASAR=[];
goodFile=[];
goodNames=[];
DASAR_str_out=[];
%if ~exist('detectiondir')
%    detectiondir=[];
%end


if strcmp(outputdir(end),'/')
    outputdir=outputdir(1:(end-1));
end

Idot=min(findstr(Icase,'.'))-1;

for KK=1:8
    year_str=sprintf('%02i',6+KK); %Search starting in 2007
    if findstr(['Shell' year_str],Icase)
        if findstr(DASAR_str,'*')>0,
            DASAR_str=upper('abcdefghijklmnopqrsuvwxyz');  %No 'T', because of 'T' in filename
        end
        DASAR_str=upper(DASAR_str);
        localdir=sprintf('%s/Site_%02i/%s/morph/',outputdir,Isite,Icase(1:Idot));
        Icount=1;
        for J=1:length(DASAR_str)
            available_DASARS=dir(sprintf('%s/*S%i*%s*%s*morph*',localdir,Isite,DASAR_str(J),date_str));
            %available_DASARS=dir(sprintf('%s/*S%i*%s*%s*',localdir,Isite,DASAR_str(J),date_str));
            
            if length(available_DASARS)>0
                
                DASAR_str_out=[DASAR_str_out DASAR_str(J)];
                goodDASAR{Icount}=available_DASARS(1).name(1:6);
                goodFile{Icount}=[localdir  available_DASARS(1).name];
                Idot=findstr(available_DASARS(1).name,'.')-1;
                goodNames{Icount}=available_DASARS(1).name(1:Idot(1));  %Remove the *.sio label
                Icount=Icount+1;
            else
                
            end
        end
    end
end
end









