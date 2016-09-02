%function DasarID=chooseDASAR(did)
% chooseDASAR.m
% In-line script to put available DASARs in a menu
% and have operator choose one
%
% Generate menu for operator to choose a Dasar

function [DasarID,id]=chooseDASAR(did)

itemlist = did(1).lbl;
for id = 2:length(did)
    itemlist = [itemlist  did(id).lbl];
end
id = menu('Choose a DASAR:',itemlist)
DasarID=char(did(id).lbl)

