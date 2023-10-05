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
    case {'GTI-M35-300-directional','M-35-NS','M-35-EW'}
        ftab =  [100 4e3 7e3 10e3 12e3 20e3];
        htab = [-192 -160 -156 -155 -155.5 -160];
        phitab = -90*ones(size(htab));
    case {'GTI-M35-300-omni','M-35-omni'}
%         ftab =  [100  3e3  3.9e3 5e3   6e3  9e3  10e3 12e3 20e3 200e3];
%         htab = [-163 -164 -163   -165.5 -165 -164 -165 -163 -164 -164];
        ftab =  [100 2e3 3e3 20e3];
        htab = [-163 -163 -164 -164];
        phitab = 0*ones(size(htab));
    case 'HTI-92WB'
        ftab = [100 50e3];
        htab = [-145 -145];
        phitab = 0*ones(size(htab));
    case 'VS-209-omni'
        ftab =  [1e3 4e3 7e3 8e3 10e3 20e3];
        htab = [-162 -163 -164 -165 -165 -165];
        phitab = 0*ones(size(htab));
    case {'VS-209-directional','VS-209-X','VS-209-Y','VS-209-Z'}
        ftab =  [1e3 4.5e3 10e3 20e3];
        htab = [-185 -171 -160 -160];
        %phitab = [-90 -90 -160 -160];
        phitab = -90*ones(size(htab));
    case 'VS-301-omni'
        % based on nominal - measured sensitivity TBD
        ftab =  [3 1e3 2e3 20e3];
        htab = [-162 -162 -162 -162];
        phitab = 0*ones(size(htab));
    case {'VS-301-directional','VS-301-X','VS-301-Y','VS-301-Z'}
        % based on nominal - measured sensitivity TBD
        ftab = 3:10:20e3;
        omega = 2*pi*ftab;

        a_sens_g = 10; % V/g

        rho = 1000;  % kg/m^3
        c = 1500; % m/s

        % the conversion from V/g to V/uPa is to DIVIDE by the factor:
        %  g * (9.8 m/s^2 /g) * (rho kg/m^3) * (c m/s) / (i * omega 1/s) * (1 x 10^6 uPa / Pa)
        g_to_uPa = -1i*9.8*rho*c*1e6./omega;
        
        h_uPa = a_sens_g./g_to_uPa; % units of V/uPa
        htab = 20*log10(h_uPa); % units of dB re V/uPa
        phitab = -90*ones(size(htab));
end

if any(f==0)
    fidx = f>=3;
    HdB(fidx) = interp1(log10(ftab),htab,log10(f(fidx)),'linear','extrap');
    HdB(~fidx) = htab(1);
    
    phi_deg(fidx) = interp1(log10(ftab),phitab,log10(f(fidx)),'linear','extrap');
    phi_deg(~fidx) = phi_deg(1);
else
    HdB = interp1(log10(ftab),htab,log10(f));
    phi_deg = interp1(log10(ftab),phitab,log10(f));
end

 H =  (10.^(-HdB/20)).*exp(1i*phi_deg*pi/180);  %%%Units of uPa/V
 
 H=  (Vmax./(2^32))*H;  %%Units of uPa/(count);
 