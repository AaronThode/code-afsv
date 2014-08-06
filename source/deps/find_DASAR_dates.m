function [goodDASAR,goodFile,goodNames,DASAR_str_out]=find_DASAR_dates(date_str,Isite,DASAR_str,rawdatadir,Icase)
%function [goodDASAR,goodFile,goodNames,DASAR_str]=find_DASAR_dates(date_str,Isite,DASAR_str,rawdatadir,Icase)
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
%       goodNames: short name of file, minus the pathname...
%       DASAR_str_out: string of letters representing DASARS actually present
%           in data...
%     For AURALS:
%       goodDASAR: cell array of files that were recorded on this single data
%       goodFile: full path name and filename of each file that fits criteria
%       goodNames: short name of file, minus the pathname and extension...
%       DASAR_str_out: string of letters representing DASARS actually present
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


if strcmp(rawdatadir(end),'/')
    rawdatadir=rawdatadir(1:(end-1));
end

%%Shell...

for KK=1:20 
    year_str=sprintf('%02i',6+KK); %Search starting in 2007
    if findstr(['Shell' year_str],Icase)
        if findstr(DASAR_str,'*')>0,
            DASAR_str=upper('abcdefghijklmnopqrstuvwxyz');
        end
        DASAR_str=upper(DASAR_str);
        disp(sprintf('Searching for %s/S%i%sgsif/S%i*',rawdatadir,Isite,year_str,Isite));
        available_DASARS=dir(sprintf('%s/S%i%sgsif/S%i*',rawdatadir,Isite,year_str,Isite));
        Icount=1;
        DASAR_str_out=[];
        for I=1:length(available_DASARS)
            DASARname=available_DASARS(I).name;
            if findstr(DASARname(end-1),DASAR_str)>0
                localdir=sprintf('%s/S%i%sgsif/%s',rawdatadir,Isite,year_str,DASARname);
                fnames=dir(sprintf('%s/*%s*',localdir,date_str));
                if length(fnames)==1
                    DASAR_str_out=[DASAR_str_out DASARname(end-1)];
                    goodDASAR{Icount}=DASARname;
                    goodFile{Icount}=[localdir '/' fnames(1).name];
                    Idot=findstr(fnames(1).name,'.')-1;
                    goodNames{Icount}=fnames(1).name(1:Idot(1));  %Remove the *.sio label
                    Icount=Icount+1;
                end
            end
       end
   end
end

%AURALs...

for KK=1:20 
    year_str=sprintf('%02i',6+KK); %Search starting in 2007
    if findstr(['AURAL_OW' year_str],rawdatadir)
        %if findstr(DASAR_str,'*')>0,
        %    DASAR_str=upper('abcdefghijklmnopqrstuvwxyz');
        %end
        %DASAR_str=upper(DASAR_str);
        %         available_DASARS=dir(sprintf('%s/A%i/S%i*',rawdatadir,Isite,year_str,Isite));
        Icount=1;
        %         DASAR_str_out=[];
        %         for I=1:length(available_DASARS),
        %             DASARname=available_DASARS(I).name;
        %             if findstr(DASARname(end-1),DASAR_str)>0,
        localdir=sprintf('%s/A%i/',rawdatadir,Isite);
        fnames=dir(sprintf('%s/A%i/AURAL*%s*.wav',rawdatadir,Isite,date_str));
        for fi=1:length(fnames)
            %DASAR_str_out=[DASAR_str_out DASARname(end-1)];
            goodDASAR{Icount}=fnames(fi).name;
            goodFile{Icount}=[localdir '/' fnames(fi).name];            %was 1 instead of fi but changed to fi so that fname (goodFile) changed along with shortname (goodNames)
            Idot=findstr(fnames(fi).name,'.')-1;
            goodNames{Icount}=fnames(fi).name(1:Idot(1));  %Remove the *.wav label
            Icount=Icount+1;
        end
    end
end

%end

%%%%%%%%%%%BP%%%%%%%%%%%%%%%%%%
for KK=1:8 
    year_str=sprintf('%02i',6+KK); %Two-character Search starting in 2007
    if findstr(['BP' year_str],Icase)
        if findstr(DASAR_str,'*')>0,
            DASAR_str=upper('abcdefghijklmnopqrstuvwxyz');
        end
        DASAR_str=upper(DASAR_str);
        
        available_DASARS=dir(sprintf('%s/NA%sgsif/NA%s*',rawdatadir,year_str,year_str));
        
        
        Icount=1;
        DASAR_str_out=[];
        for I=1:length(available_DASARS),
            DASARname=available_DASARS(I).name;
            if findstr(DASARname(end-1),DASAR_str)>0,
                localdir=sprintf('%s/NA%sgsif/%s',rawdatadir,year_str,DASARname);
                fnames=dir(sprintf('%s/*%s*',localdir,date_str));
                if length(fnames)==1,
                    DASAR_str_out=[DASAR_str_out DASARname(end-1)];
                    goodDASAR{Icount}=DASARname;
                    goodFile{Icount}=[localdir '/' fnames(1).name];
                    Idot=findstr(fnames(1).name,'.')-1;
                    goodNames{Icount}=fnames(1).name(1:Idot(1));  %Remove the *.sio label
                    Icount=Icount+1;
                    
                end
            end
        end
        
        
    end
end




if findstr('NorthStar08',Icase)
    disp('Using NorthStar08 conventions');

    if findstr(DASAR_str,'*')>0,
        DASAR_str=upper('abcdefghijklmnopq');
    end
    DASAR_str=upper(DASAR_str);

    available_DASARS=dir(sprintf('%s/NA08*',rawdatadir));


    Icount=1;
    DASAR_str_out=[];
    for I=1:length(available_DASARS),
        DASARname=available_DASARS(I).name;
        if findstr(DASARname(end-1),DASAR_str)>0,
            localdir=sprintf('%s/%s',rawdatadir,DASARname);

            fnames=dir(sprintf('%s/*%s*',localdir,date_str));
            if length(fnames)==1,
                DASAR_str_out=[DASAR_str_out DASARname(end-1)];
                goodDASAR{Icount}=DASARname;
                goodFile{Icount}=[localdir '/' fnames(1).name];
                Idot=findstr(fnames(1).name,'.')-1;
                goodNames{Icount}=fnames(1).name(1:Idot(1));  %Remove the *.sio label
                Icount=Icount+1;

            end
        end
    end
end






end





