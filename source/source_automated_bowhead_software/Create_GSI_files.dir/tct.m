function ds = tct(ct)
%tct.m function to translate c-time to matlab time vector
mt = c2mat_tm(ct);
ds = datestr(mt);
