function H = getSensitivity(f,sensor)
% applies the tabulated sensitivity for the chosen sensor
% INPUTS:
%   f:      desired frequencies in Hz
%   sensor: desired sensor ('GTI-M35-300-omni', 'GTI-M35-300-directional',
%               or 'HTI-92WB'
% OUTPUTS:
%   H:      sensivity values for in dB at frequencies f

% Alison B. Laferriere

Vmax=5;  %Assumed peak-to peak range of the ADC

switch sensor
    case 'GTI-M35-300-directional'
        ftab =  [100 4e3 7e3 10e3 12e3 20e3];
        htab = [-192 -160 -156 -155 -155.5 -160];
    case 'GTI-M35-300-omni'
%         ftab =  [100  3e3  3.9e3 5e3   6e3  9e3  10e3 12e3 20e3 200e3];
%         htab = [-163 -164 -163   -165.5 -165 -164 -165 -163 -164 -164];
        ftab =  [100 2e3 3e3 20e3];
        htab = [-163 -163 -164 -164];
    case 'HTI-92WB'
        ftab = [100 50e3];
        htab = [-145 -145];
end
if any(f==0)
    fidx = f>=100;
    H(fidx) = interp1(log10(ftab),htab,log10(f(fidx)),'linear','extrap');
    H(~fidx) = htab(1);
else
    H = interp1(log10(ftab),htab,log10(f));
end

 H =  (10.^(-H/20));  %%%Units of uPa/V
 
 H=  (Vmax./(2^32))*H;  %%Units of uPa/(count);
 