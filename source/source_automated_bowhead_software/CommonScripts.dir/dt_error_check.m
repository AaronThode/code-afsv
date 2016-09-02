function [dt_err,dt_est]=dt_error_check(VM,DASAR_coords,dt_true)
cwater=1490;
%%Note, do not include DASAR_coords and actual_dtthat you do not want used..
Nunits=size(DASAR_coords,1);
dt_err=NaN*ones(Nunits,1);
dt_est=dt_err;

ranges=sqrt(sum((ones(Nunits,1)*VM-DASAR_coords).^2,2));
%Ianchor=find(locations{I}.dt==0);
Ianchor=min(find(dt_true==0));
if isempty(Ianchor)
    dt_err=Inf*ones(Nunits,1);
    return
end
dt_est=(ranges(Ianchor)-ranges)./cwater;
dt_err=abs(dt_true'-dt_est);

end
