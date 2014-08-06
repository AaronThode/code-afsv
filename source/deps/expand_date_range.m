%%expand_date_range.m%%
%% Given two dates, fill in all days in between as year/month/day strings.
%% e.g. '20070829'
%function date_range=expand_date_range(date_str,Site_vector),
%  Input:
%   date_str: two element cell array containing two date strings
%  Output:
%   date_range, matrix of date strings (row is date).
%   date_range_cell: date strings as a cell matrix (useful for using menu
%       command).
%Sept. 7 2008


function [date_range,date_range_cell]=expand_date_range(date_str)

date_range=[];
Idate=1;
for K=1:length(date_str{1})
    for I=1:2
        dateval(I)=datenum(str2num(date_str{I}{K}(1:4)),str2num(date_str{I}{K}(5:6)),str2num(date_str{I}{K}(7:8)));
    end

    date_vec=dateval(1):datenum(0,0,1,0,0,0):dateval(2);
    date_range=[date_range ; datestr(date_vec,'yyyymmdd')];
    for I=1:length(date_vec)
        date_range_cell{Idate}=datestr(date_vec(I),'yyyymmdd');
        Idate=Idate+1;
    end

end