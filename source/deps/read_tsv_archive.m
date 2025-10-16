%%%read_tsv_archive.m%%%%%
%function [ind,localized]=read_tsv_archive(fname,ctmin,ctmax,DASAR_list);
%  Based on documentation in WCoutFormat03.doc by Bill McLennan
%  Use this to read files labeled 'manual_archive' because additional
%   fields beyond Bill's document exist.
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
% localized.ctev(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;  %ctime date
% localized.atev(Icall,1)=datenum(tabfield(tline,J)).';J=J+1; %Date string
% localized.ctev_UTC(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;  %ctime date UTC
% localized.atev_UTC(Icall,1)=datenum(tabfield(tline,J)).';J=J+1; %Date string UTC

% localized.utmx(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.utmy(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.utm_zone(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;

% localized.lat(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.long(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.range(Icall,1)=str2double(tabfield(tline,J)).';J=J+1; %Range to nearest DASAR  

% localized.wctype(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.area(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.axmajor(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.axminor(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
%
% localized.Baxis(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.Ang(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
% localized.Nused(Icall,1)=str2double(tabfield(tline,J)).';J=J+5;
% localized.Outcome{Icall,1}=(tabfield(tline,J)).';J=J+1;
% localized.comment{Icall,1}=(tabfield(tline,J)).';J=J+1;
% localized.OperatorID{Icall,1}=(tabfield(tline,J)).';J=J+1;
% localized.DateTimeProc{Icall,1}=(tabfield(tline,J)).';J=J+1;
%
% ind.wgt(Icall,Idd)=str2double(tabfield(tline(Igood:end),2));
% ind.bearing(Icall,Idd)=str2double(tabfield(tline(Igood:end),4));
% ind.bref(Icall,Idd)=str2double(tabfield(tline(Igood:end),3));
% ind.ctime(Icall,Idd)=str2double(tabfield(tline(Igood:end),5));
% ind.sigdb(Icall,Idd)=str2double(tabfield(tline(Igood:end),6));
% ind.stndb(Icall,Idd)=str2double(tabfield(tline(Igood:end),7));
% ind.flo(Icall,Idd)=str2double(tabfield(tline(Igood:end),8));
% ind.fhi(Icall,Idd)=str2double(tabfield(tline(Igood:end),9));
% ind.duration(Icall,Idd)=str2double(tabfield(tline(Igood:end),10));
% ind.sdm(Icall,Idd)=str2double(tabfield(tline(Igood:end),11));
% ind.kappa(Icall,Idd)=str2double(tabfield(tline(Igood:end),12));

function [ind,localized]=read_tsv_archive(fname,ctmin,ctmax,DASAR_list)
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
            fprintf('Bad line: %s\n',tline);
            break
        end
        % read computed time of call
        ctime=str2double(tabfield(tline,1));
        % wctype=str2double(tabfield(tline,5));
        %compute bin number
        %bin=ceil((ctime -ctstart)/binsize);
        %%Aaron addition

        try
            if (ctime>=ctmin)&&(ctime<=ctmax) || isnan(ctime)
                Icall=Icall+1;
                if rem(Icall,1000)==0, disp(Icall);end
                J=1;
                localized.ctev(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                try
                    localized.atev(Icall,1)=datenum(tabfield(tline,J));J=J+1;
                catch
                    localized.atev(Icall,1)=NaN;J=J+1;
                end
                %  J=J+2;
                localized.ctev_UTC(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                try
                    localized.atev_UTC(Icall,1)=datenum(tabfield(tline,J));J=J+1;
                catch
                    localized.atev_UTC(Icall,1)=NaN;J=J+1;
                end


                localized.utmx(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.utmy(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.utm_zone{Icall,1}=(tabfield(tline,J));J=J+1;

                localized.lat(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.long(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.range(Icall,1)=str2double(tabfield(tline,J)).';J=J+1; %Range to nearest DASAR

                localized.wctype(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.nseq(Icall,1)=0;
                if localized.wctype(Icall,1)-floor(localized.wctype(Icall,1))>0  %If there is a decimal in the call type
                    tmp=tabfield(tline,J-1);
                    Idot=strfind(tmp,'.')+1;
                    localized.nseq(Icall,1)=str2double(tmp(Idot:end));
                end
                localized.area(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.axmajor(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.axminor(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;

                localized.Baxis(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.Ang(Icall,1)=str2double(tabfield(tline,J)).';J=J+1;
                localized.Nused(Icall,1)=str2double(tabfield(tline,J,1)).';J=J+5;
                localized.Outcome{Icall,1}=(tabfield(tline,J));J=J+1;
                localized.comment{Icall,1}=(tabfield(tline,J));J=J+1;
                localized.OperatorID{Icall,1}=(tabfield(tline,J));J=J+1;
                localized.DateTimeProc{Icall,1}=(tabfield(tline,J));J=J+1;

                %count(bin)=count(bin)+1;
                for Idd=1:length(DASAR_list)
                    Igood=strfind(tline,DASAR_list{Idd});
                    if isempty(Igood)
                        DASAR_temp=DASAR_list{Idd};
                        DASAR_temp(end)='1';
                        Igood=strfind(tline,DASAR_temp);

                    end
                    if isempty(Igood)
                        continue
                    end

                    %else
                    ind.wgt(Icall,Idd)=str2double(tabfield(tline(Igood:end),2));
                    ind.bref(Icall,Idd)=str2double(tabfield(tline(Igood:end),3));

                    ind.bearing(Icall,Idd)=str2double(tabfield(tline(Igood:end),4));
                    ind.ctime(Icall,Idd)=str2double(tabfield(tline(Igood:end),5));
                    %ind.wctype(Icall,Idd)=localized.wctype(Icall);
                    ind.sigdb(Icall,Idd)=str2double(tabfield(tline(Igood:end),6));
                    ind.stndb(Icall,Idd)=str2double(tabfield(tline(Igood:end),7));
                    ind.flo(Icall,Idd)=str2double(tabfield(tline(Igood:end),8));
                    ind.fhi(Icall,Idd)=str2double(tabfield(tline(Igood:end),9));
                    ind.duration(Icall,Idd)=str2double(tabfield(tline(Igood:end),10));
                    ind.sdm(Icall,Idd)=str2double(tabfield(tline(Igood:end),11));
                    ind.kappa(Icall,Idd)=str2double(tabfield(tline(Igood:end),12));

                    %keyboard;
                    %end
                end


            else  %If not in time span
                display(tline);
                badbins=badbins+1;
                % display(['bin,wctype,badbins='  num2str(bin) blanks(2) num2str(wctype_all(Icall)) blanks(2) num2str(badbins)])
            end
        catch
            disp('try problem');
            disp(tline);
            Icall=Icall-1;
            keyboard

        end %try
    end   %of file read loop(while)
    fclose(fid);

    %%If no detections, return
    if isempty(ind)
        return
    end
    %Just in case the last DASAR in DASAR_list had no detections, make sure
    %number of columns in ind same as DASAR_list
    ncol_needed=length(DASAR_list)-size(ind.ctime,2);

    if ncol_needed>0
        disp('Fixing column length because last DASAR had no detections');
        Nrow=size(ind.bearing,1);
        ind.wgt=[ind.wgt zeros(Nrow,ncol_needed)];
        ind.bearing=[ind.bearing zeros(Nrow,ncol_needed)];
        ind.bref=[ind.bref zeros(Nrow,ncol_needed)];
        ind.ctime=[ind.ctime zeros(Nrow,ncol_needed)];
        %ind.wctype=[ind.wctype zeros(Nrow,ncol_needed)];
        ind.sigdb=[ind.sigdb zeros(Nrow,ncol_needed)];
        ind.stndb=[ind.stndb zeros(Nrow,ncol_needed)];
        ind.flo=[ind.flo zeros(Nrow,ncol_needed)];
        ind.fhi=[ind.fhi zeros(Nrow,ncol_needed)];
        ind.duration=[ind.duration zeros(Nrow,ncol_needed)];
        ind.sdm=[ind.sdm zeros(Nrow,ncol_needed)];
        ind.kappa=[ind.kappa zeros(Nrow,ncol_needed)];

    end

    %%%Replaces zeros with NaN when no data present to avoid confusion
    fieldnamess=fieldnames(ind);
    Ibad=find(ind.ctime==0);
    for JJ=1:length(fieldnamess)
        ind.(fieldnamess{JJ})(Ibad)=NaN;
    end
end

function test
fname='TSV_files_Shell07/wc090107s3.tsv';
[ind,localized]=read_tsv(fname,0,Inf,{'D07s3a','D07s3b','D07s3c','D07s3d','D07s3e','D07s3f','D07s3g'});

%function fs = TabField(s,n,nc)
%
%  Where s is a string containing Tab-delimited fields,
%  this function will return the nth field, with trailing
%  and leading spaces removed.  If n is an array, then this
%  function returns a cell array, in which each cell is
%  the field indicated by the corresponding element in n.

function fs = tabfield(s,n,nc)
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
    if exist('nc','var')
        fs=fs(nc);
    end
else
    fs = cell(size(n));
    for ii = 1:prod(size(n))
        fs{ii} = TabField(s,n(ii));
    end
end
return