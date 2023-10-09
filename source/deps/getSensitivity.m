function H = getSensitivity(FF,sensor,varargin)
% function H = getSensitivity(f,sensor,nbits,Vmax)
% applies the tabulated sensitivity for the chosen sensor pressure (not
% power) spectrum
% INPUTS:
%   f:      desired frequencies in Hz, a row vector
%   sensor: desired sensor ('GTI-M35-300-omni', 'GTI-M35-300-directional',
%               or 'HTI-92WB'
%   additional name-value pair inputs: 'nbits', 'vmax','units' can be specified as:
%       getSensitivity(f,'GTI-M35-300-omni','nbits',16,'Vmax',2,'uPa/V') etc.
%
% OUTPUTS:
%   H:      sensivity values for dB (amplitude) at frequencies f; a row vector

% Alison B. Laferriere

p = inputParser;

p.PartialMatching = true; % this makes it case insensitive

addParameter(p,'nbits',24);
addParameter(p,'vmax',5);
addParameter(p,'units','uPa/count') % or 'uPa/V'
parse(p,varargin{:});

switch sensor
    case 'DASAR-omni'  %%%Remember energy is arriving as normal modes so a slight vertical offset.

        DASAR_chc='Aaron';  %Aaron or Alison
        H=ones(1,length(FF));

        switch DASAR_chc
            case 'Aaron'
                %%%% slopee measured at 00:06:30 local time on 17 August, 2014 5G
                %                 %%%%23:57:28  1 seconds 10/1/2014 5G
                %                 %%% Also checked 04-Oct-2010 02:45:47.4 DASAR 5G, within 10° below
                %                 %%%     375 Hz.  This signal is lower received level.
                %                 %slopee=1.1*4.1e-04;  %%%radians phase per radians frequency, for
                %                 slopee=1.2*4.1e-04;  %%%radians phase per radians frequency, for 1 sec averages
                %
                slopee=1.2*4.1e-04;  %%%radians phase per radians frequency, for 1 sec averages


                %Phasee=10*pi/180+slopee*2*pi*(FF-75);  %Phase is 10 degrees at 75 Hz, linear
                Phasee=slopee*2*pi*FF;  %Nearly identical to above (y-intercept is
                %           -1 degrees in equation above.


            case 'Alison'
                fe = [0 93 148 200 250 360 500]';
                pex = [0 -3 -7.8 -18.3 -25 -60 -90]';
                phaseex = interp1(fe,pex,FF);
                phaseex(isnan(phaseex)) = 0;
                Phasee=exp(1i*(pi/180)*phaseex);
        end
        H=H.*exp(-1i*Phasee.');

        return
    case {'DASAR-directional','DASAR-X','DASAR-Y'}
        %%%%Thought normal modes might explain this, but better using a
        %%%%frequency-independent factor.
        H=ones(1,length(FF));

        D=18; c=1500;
        k=2*pi*FF/c;
        cos_elevation=real(sqrt(k.^2-(pi/D).^2)./k);
        f_trans=150;
        cos_elevation(FF>=f_trans)=cos_elevation(FF>=f_trans).*0.8;

        phasse=-7*ones(size(FF));
        phasse(FF>=350)=4;
        %Gains(:,2:3)=Gains(:,2:3).*cos_elevation.*exp(1i*2*pi*(phasse)/180); %10-Apr-2020 00:02:03.000 287 (110 and -6 null) deg Arctic5G_2014 5G 100-150 Hz 17 m water depth
        
        %H=H.*cosd(25).*exp(1i*2*pi*(phasse.')/180); %10-Apr-2020 00:02:03.000 287 (110 and -6 null) deg Arctic5G_2014 5G 100-150 Hz 17 m water depth

        return

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

    otherwise
        H=ones(1,length(FF));
        return
end
if any(FF==0)
    fidx = FF>=3;
    HdB(fidx) = interp1(log10(ftab),htab,log10(FF(fidx)),'linear','extrap');
    HdB(~fidx) = htab(1);

    phi_deg(fidx) = interp1(log10(ftab),phitab,log10(FF(fidx)),'linear','extrap');
    phi_deg(~fidx) = phi_deg(1);
else
    HdB = interp1(log10(ftab),htab,log10(FF));
    phi_deg = interp1(log10(ftab),phitab,log10(FF));
end

H =  (10.^(-HdB/20)).*exp(1i*phi_deg*pi/180);  %%%Units of uPa/V

switch p.Results.units
    case 'uPa/count'
        % convert to counts
        H=  (p.Results.vmax./(2^(p.Results.nbits-1)))*H;  %%Units of uPa/(count);
    case 'uPa/V'
        % do nothing
    otherwise
        error('unrecognized units!')
end
