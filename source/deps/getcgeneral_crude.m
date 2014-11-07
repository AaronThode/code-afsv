
%%%%%%%%%%%%%%%%%%%%getcgeneral_crude.m%%%%%%%%%%%%%%%%
% Create a linear sound speed profile
function [c,dcdz]=getcgeneral_crude(z)
persistent Data

dz=0.1;

if isempty(Data)
    %Data=load('../MetaData.dir/SSP_smooth.mat');
    Data=load('SSP_MURI_Feb2014_smooth');
end

if z<0
    z=0.05;
end

if z>=(max(Data.Z)-dz)
    z=max(Data.Z)-2*dz;
    
end

c=interp1(Data.Z,Data.C,z);
cp=interp1(Data.Z,Data.C,z+dz);
dcdz=(cp-c)/dz;