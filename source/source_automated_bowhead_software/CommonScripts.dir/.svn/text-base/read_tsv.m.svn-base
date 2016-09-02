%%%read_tsv.m%%%%%
%function [ind,localized]=read_tsv(fname,ctmin,ctmax,DASAR_list);
%  Based on documentation in WCoutFormat03.doc by Bill McLennan
%Input:
%   fname:  full filename (including path) of *.tsv file.  Must have tsv
%       extension in string.  e.g. {'D07s3d'}
%   ctmin: minimum c-time desired
%   ctmax: maximum c-time desired
%   DASAR_list: cell matrix of DASAR names to extract
%Output:
%   ind: row is call number, column is corresponding station in DASAR_list
%       if particular call not found in DASAR, ind.wctype will be 0
% %   localized:  localization information about call across all DASARS, see below.
% localized.ctev(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.atev(Icall,1)=datenum(tabfield(tline,J)).';J=J+1;
% localized.utmx(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.utmy(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.wctype(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.area(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.axmajor(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.axminor(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% 
% localized.Baxis(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.Ang(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
% localized.Nused(Icall,1)=str2num(tabfield(tline,J)).';J=J+5;
% localized.Outcome{Icall,1}=(tabfield(tline,J)).';J=J+1;
% localized.comment{Icall,1}=(tabfield(tline,J)).';J=J+1;
% localized.OperatorID{Icall,1}=(tabfield(tline,J)).';J=J+1;
% localized.DateTimeProc{Icall,1}=(tabfield(tline,J)).';J=J+1;
%  %count(bin)=count(bin)+1;
%  for Idd=1:length(DASAR_list),
%      Igood=findstr(DASAR_list{Idd},tline);
%      if ~isempty(Igood),
%         ind.ctime(Icall,Idd)=str2num(tabfield(tline(Igood:end),5));
%         ind.wctype(Icall,Idd)=localized.wctype(Icall);
%         ind.sigdb(Icall,Idd)=str2num(tabfield(tline(Igood:end),6));
%         ind.stndb(Icall,Idd)=str2num(tabfield(tline(Igood:end),7));
%         ind.flo(Icall,Idd)=str2num(tabfield(tline(Igood:end),8));
%         ind.fhi(Icall,Idd)=str2num(tabfield(tline(Igood:end),9));
%         ind.duration(Icall,Idd)=str2num(tabfield(tline(Igood:end),10));
%         %keyboard;
%     end
% end
function [ind,localized]=read_tsv(fname,ctmin,ctmax,DASAR_list)
fid=fopen(fname,'r');
badbins=0;
nline=0;
Icall=0;
ind=[];
localized=[];
if fid ~= -1    % did file open?
    while feof(fid)==0
        tline=fgets(fid);
        nline=nline+1;
        if tline ==-1
            break
        end
        % read computed time of call
        ctime=str2num(tabfield(tline,1));
        % wctype=str2num(tabfield(tline,5));
        %compute bin number
        %bin=ceil((ctime -ctstart)/binsize);
        %%Aaron addition


        if (ctime>=ctmin)&(ctime<=ctmax) %& wctype <=8
            Icall=Icall+1;
            if rem(Icall,500)==0, disp(Icall);end
            J=1;
            localized.ctev(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.atev(Icall,1)=datenum(tabfield(tline,J)).';J=J+1;
            localized.utmx(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.utmy(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.wctype(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.nseq(Icall,1)=0;
            if localized.wctype(Icall,1)-floor(localized.wctype(Icall,1))>0
                tmp=tabfield(tline,J-1);
                Idot=findstr(tmp,'.')+1;
                localized.nseq(Icall,1)=str2num(tmp(Idot:end));
            end
            localized.area(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.axmajor(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.axminor(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;

            localized.Baxis(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.Ang(Icall,1)=str2num(tabfield(tline,J)).';J=J+1;
            localized.Nused(Icall,1)=str2num(tabfield(tline,J)).';J=J+5;
            localized.Outcome{Icall,1}=(tabfield(tline,J)).';J=J+1;
            localized.comment{Icall,1}=(tabfield(tline,J)).';J=J+1;
            localized.OperatorID{Icall,1}=(tabfield(tline,J)).';J=J+1;
            localized.DateTimeProc{Icall,1}=(tabfield(tline,J)).';J=J+1;

            %count(bin)=count(bin)+1;
            for Idd=1:length(DASAR_list),
                Igood=findstr(DASAR_list{Idd},tline);
                if ~isempty(Igood),
                    ind.wgt(Icall,Idd)=str2num(tabfield(tline(Igood:end),2));
                    ind.bearing(Icall,Idd)=str2num(tabfield(tline(Igood:end),4));
                    ind.bref(Icall,Idd)=str2num(tabfield(tline(Igood:end),3));
                    ind.ctime(Icall,Idd)=str2num(tabfield(tline(Igood:end),5));
                    ind.wctype(Icall,Idd)=localized.wctype(Icall);
                    ind.sigdb(Icall,Idd)=str2num(tabfield(tline(Igood:end),6));
                    ind.stndb(Icall,Idd)=str2num(tabfield(tline(Igood:end),7));
                    ind.flo(Icall,Idd)=str2num(tabfield(tline(Igood:end),8));
                    ind.fhi(Icall,Idd)=str2num(tabfield(tline(Igood:end),9));
                    ind.duration(Icall,Idd)=str2num(tabfield(tline(Igood:end),10));
                    ind.sdm(Icall,Idd)=str2num(tabfield(tline(Igood:end),11));
                    ind.kappa(Icall,Idd)=str2num(tabfield(tline(Igood:end),12));
                    
                    %keyboard;
                end
            end

        else
            % display(tline);
            badbins=badbins+1;
            % display(['bin,wctype,badbins='  num2str(bin) blanks(2) num2str(wctype_all(Icall)) blanks(2) num2str(badbins)])
        end
    end   %of file read loop(while)
    fclose(fid);
    
    %%If no detections, return
    if isempty(ind)
        return
    end
    %Just in case the last DASAR in DASAR_list had no detections, make sure
    %number of columns in ind same as DASAR_list
    ncol_needed=length(DASAR_list)-size(ind.ctime,2);
    
    if ncol_needed>0,
        ind.wgt=[ind.wgt zeros(Icall,ncol_needed)];
        ind.bearing=[ind.bearing zeros(Icall,ncol_needed)];
        ind.bref=[ind.bref zeros(Icall,ncol_needed)];
        ind.ctime=[ind.ctime zeros(Icall,ncol_needed)];
        ind.wctype=[ind.wctype zeros(Icall,ncol_needed)];
        ind.sigdb=[ind.sigdb zeros(Icall,ncol_needed)];   
        ind.stndb=[ind.stndb zeros(Icall,ncol_needed)];
        ind.flo=[ind.flo zeros(Icall,ncol_needed)];
        ind.fhi=[ind.fhi zeros(Icall,ncol_needed)];
        ind.duration=[ind.duration zeros(Icall,ncol_needed)];
        ind.sdm=[ind.sdm zeros(Icall,ncol_needed)];
        ind.kappa=[ind.kappa zeros(Icall,ncol_needed)];

    end
        
end

function test
fname='TSV_files_Shell07/wc090107s3.tsv';
[ind,localized]=read_tsv(fname,0,Inf,{'D07s3a','D07s3b','D07s3c','D07s3d','D07s3e','D07s3f','D07s3g'});

%function fs = TabField(s,n)
%
%  Where s is a string containing Tab-delimited fields,
%  this function will return the nth field, with trailing
%  and leading spaces removed.  If n is an array, then this
%  function returns a cell array, in which each cell is
%  the field indicated by the corresponding element in n.

function fs = tabfield(s,n)
%mbcharvector(s);
%mbint(n);
if isempty(n)
    fs = '';
elseif length(n) == 1
    xx = find((s == char(9)));

    if n > (length(xx) + 1)
        fs = [];
    else
        if isempty(xx)
            fs = s;
        else
            xx = [0;xx(:);length(s) + 1];
            fs = s((xx(n) + 1):(xx(n + 1) - 1));
        end
    end
    fs = strtrim(fs);
else
    fs = cell(size(n));
    for ii = 1:prod(size(n))
        fs{ii} = TabField(s,n(ii));
    end
end
return