function [brefa_table,Icol]=calibrate_bearing_Shell2007(cal_dir,fname,no_load_table)
%function [brefa_table,Icol]=calibrate_bearing_Shell2007(cal_dir,fname,no_load_table)
% Return nonlinear calibration data for 2007 Dasars
% Input:
%          cal_dir: base directory of calibration information
%          fname: string similar in structure to 'S108H0T20080921T000000'
%          no_load_table:  if exists, don't load table.
% Note:  the calibration is assigned based on name of file, not name of enclosing directory..
brefa_table=[];
Islash=1+max(strfind(fname,'/'));
if ~isempty(Islash)
    fname=fname(Islash:end);
end
site=str2double(fname(2));
dasar=lower(fname(5));
Icol=double(dasar)-96;


%%Date checking...
%% In 2007 if a DASAR was redeployed it was reassigned a number...
% if site==3&&strcmp(dasar,'g')
%     datte=datenum(fname(8:end),'yyyymmddTHHMMSS');
%     if datte>=datenum(2007,9,18,3,15,0)
%         Icol=8;
%     end
% elseif site==5&&strcmp(dasar,'d')
%     datte=datenum(fname(8:end),'yyyymmddTHHMMSS');
%     if datte>=datenum(2007,9,18,19,0,0)
%         Icol=8;
%     end
% elseif site==5&&strcmp(dasar,'e')
%     datte=datenum(fname(8:end),'yyyymmddTHHMMSS');
%     if datte>=datenum(2007,9,19,0,0,0)
%         Icol=9;
%     end
% end


if nargin==2
    cal_dir=sprintf('%s/S%i07gsif/Caloutput07s%i/Bgrid07s%i.mat',cal_dir,site,site,site);
    data=load(cal_dir);
    brefa_table=data.bgridt;
    
    %     if site==1  %DASAR 1B failed, not reflected in bgrid table...
    %        brefa_table=[brefa_table(:,1) zeros(size(bgridt,1),1) brefa_table(:,2:end)];
    %     end
end

end %calibrate_bearing_Shell2007
