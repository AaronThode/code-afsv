
%%%%%%%%%%%%TOC_env.m%%%%%%%%%%
%function [lthick,svp,spopt,cmedtop,cmedbot,rho,pwavatten,plottitle,case_ssp,case_bottom]=TOC_env(case_ssp,case_bottom,D)
%load environmental parameters for run
%lthick=[D 30 800];  %Thicknesses of layers in meters
%cmedtop=[1500 1572 1881 5200]; %sound speed at tops of layers in m/sec
%cmedbot=[1500 1593 3245 0];  %sound speed at bottom layers in m/sec
%rho=[1 1.76 2.06 2.66];  %densities of layers in g/cc
%pwavatten=[0 .200 .06 .02]; %Attenuations in dB/wavelength
%	These arrays are Nmedia+1 long because the last elements describe the
%	acoustic halfspace.
%  A second row contains corresponding shear properties...
%
% Output:
%   rms_roughness: roughness of bottom in meters
function [lthick,svp,sspopt,cmedtop,cmedbot,rho, ...
    pwavatten,plottitle,case_ssp,case_bottom,rms_roughness]=TOC_env(case_ssp,case_bottom,D,bath)

rms_roughness=[];
if strcmp(case_ssp,'select_from_menu'),
    menu_chc={'pekeris_kupdemo','pekeris','GOM_822','GOM_mean','GOM_911','GOM_909','GOM_fresh', ...
        'GOM_GDEM','GOM_GDEM_Aaron','Ewing_03032_xt5','Ewing_0303_6_x7d','2015 Cabo San Lucas Pt Ballena','2014 MURI'};
    Ichc=menu('Select a sound speed profile:',menu_chc)
    case_ssp=menu_chc{Ichc};
end

if strcmp(case_bottom,'select_from_menu'),
    menu_chc={'sand_pekeris','sand_kupdemo','sand_shallow','granite_gravellayer_Cabo','sand_sedlayer_Cabo','silt_pekeris','sand_English_channel','wedge_demo'};
    Ichc=menu('Select a bottom profile:',menu_chc)
    case_bottom=menu_chc{Ichc};
end


lthick=D;
cmedtop=[1500 -1];
cmedbot=[1500 0];

sspopt = 'SVW' ;   % sound speed profile options (W = dB/wavelength)
%sspopt = 'SVN' ;   % sound speed profile options (N = Nepers/m)

svp=[];

%%%Assign bottom properties...
%%Check to see if property embedded in name
if ~isempty(findstr(case_bottom,'pekeris'))&&~isnan(str2double(case_bottom((end-3):end)))
    rho = [1 1.9] ;  % specific density from COA/Hamilton
    cmedtop(2) = str2double(case_bottom((end-3):end)) ;
    
    
    
    if cmedtop(2)==1672
        pwavatten = [0 1.02] ; % p-wave attenuation [dB/wavelength] from 2010 inversions
        rho=[1 1.34];  %From 2010 pekeris inversions...
   % elseif cmedtop(2)==1707
   %     pwavatten = [0 1.09] ; % p-wave attenuation [dB/wavelength] from 2010 inversions
    %    rho=[1 1.19];  %From 2010 pekeris inversions...
%     elseif cmedtop(2)==1650
%         cmedtop(1)=1465;
%         cmedbot(1)=1465;
%         pwavatten = [0 0.1] ; % dB/lamda
%          %pwavatten = [0 0.9] ; % dB/lamda
%         rho=[1 1.8];  %
    else
        prompt={'rho','pwavatten (dB/lambda)'};
        def={'1.9','0.8'};
        tittle='Parameters for Pekeris';
        answer=inputdlg(prompt,tittle,1,def);
        rho=[1 eval(answer{1})];
        pwavatten=[0 eval(answer{2})];
        
    end
    
    
else
    switch case_bottom
        case 'sand_pekeris',
            rho=[1 1.5];
            csed =   1700;  % sediment sound velocity
            cmedtop(2)=csed;
            pwavatten =[0 0.2];   % attenuation in db/lam
            pwavatten=[0 0.15]; %attenuation from db/m/Hz to db/wavelength
            pwavatten=pwavatten*cmedtop(2)/1000;
        case 'wedge_demo'
            rho=[1 1.5];
            csed =   1700;  % sediment sound velocity
            cmedtop(2)=csed;
            pwavatten=[0 0.5]; %attenuation from db/m/Hz to db/wavelength
            %rms_roughness=0.4;
            %pwavatten=pwavatten*cmedtop(2)/1000;
        case 'COA_pekeris_demo'
            rho=[1 1.8];
            csed =   1800;  % sediment sound velocity
            cmedtop(2)=csed;
            pwavatten=[0 1.0]; %attenuation from db/m/Hz to db/wavelength
            
        case 'sand_kupdemo'
            rho=[1 1.5];
            csed =   1600;  % sediment sound velocity
            cmedtop(2)=csed;
            pwavatten=[0 0.2]; %attenuation from db/m/Hz to db/wavelength
            %rms_roughness=0.4;
            %pwavatten=pwavatten*cmedtop(2)/1000;
        case 'sand_shallow'
            rho = [1 1.9] ;  % specific density
            csed = 1650;  % sediment p-wave velocity [m/s]
            %csed = 1586;  % sediment p-wave velocity [m/s]
            
            cmedtop(2) = csed ;
            pwavatten = [0 0.8] ; % p-wave attenuation [dB/wavelength]
            %rms_roughness=1.0;
        case 'sand_Arctic_inversion' %derived from R_360s_tiltm6_105013start_airgunonly.jpg
            
            cmedtop=[1535 1684 2526 2621 4007];
            cmedbot=[1523 2526 2621 3156 0];
            rho=[1 1.23 1.23 1.23 2.16];
            %pwavatten=[0 0.2 0.5]; %dB/lamb
            %  rho=[1 1.9 3.0];
            pwavatten=[0 0.8 0.8 0.8 0.817]; %dB/lamb
            lthick=[D 10 10 10];
            
        case 'sand_Arctic_inversionSSP_TwoClosestMFPAvg' %Average of 1.5 and 7 km inversions below,
            %               used for warping benchmarek
             cmedtop=[-1 1519 2086 2888 4211];
            cmedbot=[-1 2086 2888  3506 0];
            rho=[1 [1 1 1]*1.45 2.0];
            pwavatten=[0 0.98*[1 1 1] 0.4]; %dB/lamb my tweak for source level calculations on Oct. 2015
            
            lthick=[D 10 10 10];
        case 'sand_Arctic_inversionSSP_1p5km' %derived from 1.5 km range call on CSDM_Nfft2048_20100820T014239_eigenvector.in
            %Arctic_2010/ShortRangeMFPInversion.dir/15dBSNR_allparams.dir/BestResult.dir
            %Water Depth 57 m
            % ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_flat57m_10to500Hz
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m_fixed.dir/
            %ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_10to500Hz
            cmedtop=[-1 1513 2410 3103 4410];
            cmedbot=[-1 2410 3103 3457 0];
            rho=[1 [1 1 1]*1.31 1.7];
            %pwavatten=[0 0.2 0.5]; %dB/lamb
            %  rho=[1 1.9 3.0];
            %pwavatten=[0 1.08*[1 1 1] 0.4]; %dB/lamb  %Original MFP
            %       inversion
            pwavatten=[0 0.5*[1 1 1] 0.4]; %dB/lamb my tweak for source level calculations on Oct. 2015
            
            lthick=[D 10 10 10];
            
        case 'sand_Arctic_inversion_7kmWhale' %derived from 7 km at CSDM_Nfft2048_20100831T004830_eigenvector.in
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m.dir/
            %ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_flat55m_10to500Hz.mat
            %             cmedtop=[1491 1657 1680 1877 2578];
            %             cmedbot=[1480 1680 1877 2255 0];
            %             rho=[1 1.2 1.2 1.2 2.5];
            %             pwavatten=[0 0.5 0.5 0.5 0.3]; %dB/lamb
            %             lthick=[D 10 10 10];
            
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m_fixed.dir/
            %ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_10to500Hz
            cmedtop=[-1 1524 1761 2674 4012];
            cmedbot=[-1 1761 2674 3555 0];
            rho=[1 [1 1 1]*1.58 2.37];
            %pwavatten=[0 0.2 0.5]; %dB/lamb
            %  rho=[1 1.9 3.0];
            pwavatten=[0 0.5*[1 1 1] 0.817]; %dB/lamb
            lthick=[D 10 10 10];
            
            
        case 'sand_Arctic_inversion_15kmWhale' %derived from 7 km at CSDM_Nfft2048_20100831T004830_eigenvector.in
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m.dir/
            %ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_flat55m_10to500Hz.mat
            %             cmedtop=[1463 1664 1704 1712 3027];
            %             cmedbot=[1437 1704 1712 2672 0];
            %             rho=[1 2.0 2.0 2.0 2.1];
            %             pwavatten=[0 0.5 0.5 0.5 1.0]; %dB/lamb
            %             lthick=[D 10 10 10];
            
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m_fixed.dir/
            %ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_10to500Hz.mat
            
            cmedtop=[-1 1953 2110 2118 2583];
            cmedbot=[-1 2110 2118 2331 0];
            rho=[1 [1 1 1]*2.575 1.787];
            pwavatten=[0 0.628*[1 1 1] .604]; %dB/lamb
            lthick=[D 10 10 10];
        case 'sand_Arctic_inversion_35kmWhale' %derived from 7 km at CSDM_Nfft2048_20100831T004830_eigenvector.in
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m.dir/
            %ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105013_flat55m_10to500Hz.mat
            %             cmedtop=[1487 1532 2524 2697 3288];
            %             cmedbot=[1467 2524 2697 2847 0];
            %             rho=[1 2.575 2.575 2.575 1.7];
            %             %pwavatten=[0 0.2 0.5]; %dB/lamb
            %             %  rho=[1 1.9 3.0];
            %             pwavatten=[0 0.9 0.9 0.9 .9]; %dB/lamb
            %             lthick=[D 10 10 10];
            
            %Projects/Schultz.dir/MATLAB/PropagationModeling.dir/Depth55m_fixed.dir/
            %ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105748_10to500Hz.mat
            
            cmedtop=[-1 1590 1992 2952 3850];
            cmedbot=[-1 1992 2952 3630 0];
            rho=[1 [1 1 1]*1.13 1.96];
            pwavatten=[0 0.69*[1 1 1] .87]; %dB/lamb
            lthick=[D 10 10 10];
        case 'sand_shallow_shear'
            
            rho = [1 1.9] ;  % specific density
            csed = 1650 ;  % sediment p-wave velocity [m/s]
            cmedtop(2) = csed ;
            pwavatten = [0 0.8] ; % p-wave attenuation [dB/wavelength]
            
            cmedtop=[cmedtop; 0 200];
            cmedbot(1,2)=cmedtop(2,2);
            cmedbot(2,2)=cmedtop(2,2);
            pwavatten =[0 0.8;0 2.5];   % attenuation in db/lam
        case 'granite_gravellayer_Cabo',  %Add a sediment layer
            
            %1/Q =(alpha*v/pi*f)=0.05

            %alpha_dB=8.686*alpha*lambda
         %=8.686*lambda*(pi*f/(v*Q))=8.686*pi*(1/Q)
            
            Dsed =  10;  %  thickness of sediment layer in m
            pwavatten =[0 0.6 0.1];   % attenuation in db/lam
           
            
            lthick=[D Dsed] ;
            rho=[1 2.0 2.7];
            csed =   1800;  % sediment sound velocity
            cmedtop(2)=csed;
            sgd =  0;    %gradient
            cmedbot(2)=cmedtop(2)+sgd;
            
            cmedtop(3)=2500; %shear wave of granite
            cmedbot(3)=0;
            
            %pwavatten=[0 0.8 0.1]; %attenuation from db/m/Hz to db/wavelength
            %pwavatten=pwavatten*cmedtop(2)/1000;
            
        case 'sand_sedlayer_Cabo',  %Add a sediment layer
            
            Dsed =  25;  %  thickness of sediment layer in m
            pwavatten =[0 0.2 0.1];   % attenuation in db/lam, punta gorda simulations
            
            Dsed =  5;  %  thickness of sediment layer in m
            pwavatten =[0 0.8 0.1];   % attenuation in db/lam
            
            
  
            lthick=[D Dsed] ;
            rho=[1 1.9 2.3];
            csed =   1650;  % sediment sound velocity
            cmedtop(2)=csed;
            sgd =  0;    %gradient
            cmedbot(2)=cmedtop(2)+sgd;
            
            cmedtop(3)=2700;
            cmedbot(3)=0;
            
            %pwavatten=[0 0.8 0.1]; %attenuation from db/m/Hz to db/wavelength
            pwavatten=pwavatten*cmedtop(2)/1000;
            
        case 'silt_sedlayer',  %Add a sediment layer
            Dsed =  10;  %  thickness of sediment layer in m
            
            lthick=[D Dsed] ;
            rho=[1 1.3 1.5];
            csed =   1550;  % sediment sound velocity
            cmedtop(2)=csed;
            sgd =  10;    %gradient
            cmedbot(2)=cmedtop(2)+sgd;
            
            cmedtop(3)=1800;
            cmedbot(3)=0;
            
            %pwavatten =[0 0.3 0.1];   % attenuation in db/lam
            pwavatten=[0 0.15 0.1]; %attenuation from db/m/Hz to db/wavelength
            pwavatten=pwavatten*cmedtop(2)/1000;
        case 'silt_pekeris',
            rho=[1 1.3];
            csed =   1550;  % sediment sound velocity
            cmedtop(2)=csed;
            pwavatten =[0 0.2];   % attenuation in db/lam
            pwavatten=pwavatten*cmedtop(2)/1000;
        case 'sand_English_channel'
            rho=[1 2.0];
            csed =   1800;  % sediment sound velocity
            cmedtop(2)=csed;
            cmedtop=[cmedtop; 0 400];
            cmedbot=cmedtop;
            pwavatten =[0 0.8;0 2.0];   % attenuation in db/lam
            
            %         cmedtop=[cmedtop; 0 0];
            %         cmedbot=cmedtop;
            %         pwavatten =[0 0.8;0 0];   % attenuation in db/lam
            %         rms_roughness=0.4;
            svp=[0.0     1510
                30.0    1510
                35      1490
                D     1492
                ];
    end
    
    
    
    %lthick=[D 30 800];
    %cmedtop=[1500 1572 1881 5200];
    %cmedbot=[1500 1593 3245 0];
    %rho=[1 1.76 2.06 2.66];
    %pwavatten=[0 .200 .06 .02];
    %	These arrays are Nmedia+1 long because the last elements describe the
    %	acoustic halfspace.
    %
    
    
end  %if pekeris

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Sound speed profiles%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch case_ssp
    case 'pekeris_kupdemo',
        plottitle = 'pekeris';
        lthick(1)=D;
        cmedtop=[1500 -1];
        cmedbot=[1500 0];
        svp=[];
    case 'pekeris',
        plottitle = 'pekeris';
        lthick(1)=D;
        %cmedtop(1)=[1500];
        %cmedbot(1)=[1500];
        svp=[];
   
    case 'arctic_pekeris'
        plottitle = 'arctic_pekeris';
        lthick(1)=D;
        %cmedtop=[1535 -1];  %Original before 7/22/2011
        %cmedbot=[1523 0];   %Original before 7/22/2011
        
        cmedtop(1)=1442.00;  %
        cmedbot(1)=1442.00  ; %
        
        %cmedtop
        
        svp=[];
    case 'arcticSSP_TwoClosestMFPAvg'  %Average of two closest MPF inversions
        %       (1.5 and 7 km), with three point SSP.
        %       Using in warping benchmarking studies.
        
        plottitle = 'arctic';
        lthick(1)=D;
        
                svp=[0 1438;
                    8 1455;
                    D 1445];
    case 'arcticSSP'
        
        %Original inversion with linear profile
        %         plottitle = 'arctic';
        %         lthick=D;
        %         cmedtop=[1465 -1]
        %         cmedbot=[1457 0];
        %         svp=[];
        
        plottitle = 'arctic';
        lthick(1)=D;
        svp=[0 1441;
            D 1452];
%         svp=[0 1441;
%             8 1464;
%             55 1452];
    case 'arcticSSP_7kmWhale'
        %Original inversion with linear profile
        %         plottitle = 'arctic_7kmWhale';
        %         lthick=D;
        %         cmedtop=[1491 -1]
        %         cmedbot=[1480 0];
        %         svp=[];
        
        plottitle = 'arctic';
        lthick(1)=D;
        svp=[0 1434;
            8 1447;
            55 1439];
        
        
    case 'arcticSSP_15kmWhale'
        %Original inversion with linear profile
        
        %         plottitle = 'arctic_15kmWhale';
        %         lthick=D;
        %         cmedtop=[1463 -1]
        %         cmedbot=[1437 0];
        plottitle = 'arctic';
        lthick(1)=D;
        svp=[0 1420;
            8 1447;
            55 1431];
        
    case 'arcticSSP_35kmWhale'
        %Original inversion with linear profile
        
        %         plottitle = 'arctic_15kmWhale';
        %         lthick=D;
        %         cmedtop=[1487 -1]
        %         cmedbot=[1467 0];
        %         svp=[];
        
        plottitle = 'arctic';
        lthick(1)=D;
        svp=[0 1446;
            8 1401;
            55 1383];
        
    case 'GOM_822',
        plottitle = 'Gulf of Mexico';
        load -ascii GOM_822.m
        svp = GOM_822 ; %svp_test;
        Igood=find(svp(:,5)>5);
        svp=svp(Igood,:);
        svp=[svp(:,1) svp(:,5)];
        svp=svp(4:(end-2),:);
        svp=[0 svp(1,2);svp];
    case 'GOM_mean',
        plottitle = 'Gulf of Mexico mean profile';
        L = 4000;     %  number of layers in water column;
        %  interface at L+1
        %  L should be L > 10*D/lam = 10*D*f/c
        load -ascii meansvp.txt
        svp = meansvp ; %svp_test;
        svp=[0 svp(1,2);svp];
        Id=find(svp(:,1)<D);
        Id=Id(end);
        svp=svp(1:Id,:);
        svp=[svp;D svp(end,2)];
        
        
    case 'GOM_911',
        plottitle = 'Gulf of Mexico Sept.11, after storm';
        
        path(path,'..');
        L=4000;
        xbt0911;
        svp=M(:,[1 3]);
        svp=[0 svp(1,2);svp];
    case 'GOM_909',
        plottitle = 'Gulf of Mexico Sept. 9';
        
        path(path,'..');
        L=4000;
        xbt0909;
        svp=M(:,[1 3]);
        svp=[0 svp(1,2);svp];
    case 'GOM_fresh',
        plottitle = 'Gulf of Mexico Sept. 11 with freshwater correction';
        
        path(path,'..');
        L=4000;
        xbt0911;
        svp=M(:,[1 3]);
        svp(1:50,2)=svp(1:50,2)-5;  %29 ppt salinity at surface
        svp(51,2)=0.5*(svp(50,2)+svp(52,2));
        svp=[0 svp(1,2);svp];
    case 'GOM_GDEM',
        plottitle = 'GDEM Gulf of Mexico provided by Del';
        
        
        %svp=[0.0     1542.8836670
        svp=[0.0     1538.8836670
            10.0    1543.0607910
            20.0    1542.9898680
            30.0    1542.3852540
            50.0    1538.5949710
            75.0    1531.7530520
            100.0   1527.0809330
            125.0   1524.6839600
            150.0   1522.5419920
            200.0   1517.4738770
            250.0   1512.1301270
            300.0   1506.8452150
            400.0   1498.1040040
            500.0   1492.9677730
            600.0   1489.9194340
            700.0   1487.9512940
            800.0   1487.0299070
            900.0   1486.8637700
            1000.0  1487.1755370
            1100.0  1487.9349370
            1200.0  1488.9696040
            1300.0  1490.1162110
            1400.0  1491.5003660
            1500.0  1492.9567870
            1750.0  1497.0948490
            2000.0  1501.2102050
            2500.0  1509.0683590
            3000.0  1516.7132570
            3104.0  1518.5915530
            ];
        if isempty(bath),
            lthick=D;
        else
            lthick=max(bath(:,2));
        end
        Igood=find(svp(:,1)<lthick);
        svp_bot=interp1(svp(:,1),svp(:,2),lthick);
        svp=[svp(Igood,:); lthick svp_bot];
    case 'Sitka_measured_2010'
        plottitle = 'Sitka profile measured by Delphine in 2010';
        
        svp=[
            0 1501.699279
            10 1500.129321
            20 1495.362124
            30 1494.693787
            40 1488.814806
            50 1483.318303
            60 1482.89442
            70 1481.439363
            80 1480.371806
            90 1478.965287
            100 1478.34788
            110 1478.19982
            120 1478.034506
            130 1477.489494
            140 1476.463043
            150 1476.166747
            160 1475.594623
            170 1474.776835
            180 1473.970089
            190 1473.648046
            200 1473.586663
            210 1473.245855
            220 1473.067973
            230 1473.070671
            240 1472.678043
            250 1472.91089
            260 1472.847229
            270 1472.832474
            280 1472.900289
            290 1472.836019
            300 1472.920464
            310 1472.822331
            320 1472.723867
            330 1472.741726
            340 1472.809204
            350 1472.810153
            360 1472.894317
            380 1473.146008
            390 1473.313575
            400 1473.381034
            410 1473.548666
            420 1473.7331
            430 1473.78384
            440 1473.95157
            450 1474.119335
            460 1474.287133
            470 1474.429429
            480 1474.533429
            490 1474.637429
            500 1474.741429
            510 1474.845429
            520 1474.949429
            530 1475.053429
            540 1475.157429
            550 1475.261429
            560 1475.365429
            570 1475.469429
            580 1475.573429
            590 1475.677429
            600 1475.781429
            610 1475.885429
            620 1475.989429
            630 1476.093429
            640 1476.197429
            650 1476.301429
            660 1476.405429
            670 1476.509429
            680 1476.613429
            690 1476.717429
            700 1476.821429
            710 1476.925429
            720 1477.029429
            730 1477.133429
            740 1477.237429
            750 1477.341429
            760 1477.445429
            770 1477.549429
            780 1477.653429
            790 1477.757429
            800 1477.861429
            810 1477.965429
            820 1478.069429
            830 1478.173429
            840 1478.277429
            850 1478.381429
            860 1478.485429
            870 1478.589429
            880 1478.693429
            890 1478.797429
            900 1478.901429
            910 1479.005429
            920 1479.109429
            930 1479.213429
            940 1479.317429
            950 1479.421429
            960 1479.525429
            970 1479.629429
            980 1479.733429
            990 1479.837429
            1000 1479.941429
            1010 1480.045429
            1020 1480.149429
            1030 1480.253429
            1040 1480.357429
            1050 1480.461429
            1060 1480.565429
            1070 1480.669429
            1080 1480.773429
            1090 1480.877429
            1100 1480.981429
            1110 1481.085429
            1120 1481.189429
            1130 1481.293429
            1140 1481.397429
            1150 1481.501429
            1160 1481.605429
            1170 1481.709429
            1180 1481.813429
            1190 1481.917429
            1200 1482.021429  ];
        if isempty(bath),
            lthick=D;
        else
            lthick=max(bath(:,2));
        end
        Igood=find(svp(:,1)<lthick);
        svp_bot=interp1(svp(:,1),svp(:,2),lthick);
        svp=[svp(Igood,:); lthick svp_bot];
    case 'GOM_GDEM_Aaron',
        plottitle = 'GDEM Gulf of Mexico provided by Aaron';
        
        
        %svp=[0.0     1542.8836670
        svp=[  1 	 	 0.0 	 	 27.426 	 	 36.146 	 	1541.250
            2 		2.0 		27.404 		36.148 		1541.235
            3 		4.0 		27.383 		36.150 		1541.225
            4 		6.0 		27.351 		36.155 		1541.190
            5 		8.0 		27.308 		36.164 		1541.135
            6 		10.0 		27.263 		36.171 		1541.075
            7 		15.0 		27.128 		36.182 		1540.865
            8 		20.0 		26.945 		36.188 		1540.535
            9 		25.0 		26.663 		36.207 		1539.995
            10 		30.0 		26.316 		36.224 		1539.300
            11 		35.0 		25.913 		36.237 		1538.460
            12 		40.0 		25.498 		36.248 		1537.575
            13 		45.0 		25.082 		36.264 		1536.685
            14 		50.0 		24.694 		36.278 		1535.850
            15 		55.0 		24.297 		36.285 		1534.970
            16 		60.0 		23.935 		36.292 		1534.170
            17 		65.0 		23.595 		36.299 		1533.415
            18 		70.0 		23.266 		36.307 		1532.680
            19 		75.0 		22.953 		36.315 		1531.980
            20 		80.0 		22.675 		36.321 		1531.360
            21 		85.0 		22.417 		36.327 		1530.785
            22 		90.0 		22.157 		36.335 		1530.200
            23 		95.0 		21.907 		36.343 		1529.640
            24 		100.0 		21.664 		36.349 		1529.090
            25 		110.0 		21.166 		36.353 		1527.940
            26 		120.0 		20.683 		36.352 		1526.805
            27 		130.0 		20.213 		36.346 		1525.685
            28 		140.0 		19.724 		36.329 		1524.480
            29 		150.0 		19.252 		36.314 		1523.310
            30 		160.0 		18.806 		36.293 		1522.190
            31 		170.0 		18.368 		36.270 		1521.075
            32 		180.0 		17.930 		36.245 		1519.940
            33 		190.0 		17.500 		36.213 		1518.810
            34 		200.0 		17.088 		36.178 		1517.715
            35 		220.0 		16.309 		36.090 		1515.595
            36 		240.0 		15.561 		36.000 		1513.525
            37 		260.0 		14.898 		35.916 		1511.680
            38 		280.0 		14.294 		35.828 		1509.985
            39 		300.0 		13.738 		35.745 		1508.420
            40 		350.0 		12.443 		35.551 		1504.725
            41 		400.0 		11.301 		35.388 		1501.450
            42 		500.0 		9.414 		35.127 		1496.075
            43 		600.0 		7.922 		34.969 		1491.985
            44 		700.0 		6.872 		34.907 		1489.530
            45 		800.0 		6.049 		34.881 		1487.920
            46 		900.0 		5.481 		34.888 		1487.315
            47 		1000.0 		5.073 		34.907 		1487.345
            48 		1100.0 		4.745 		34.927 		1487.695
            49 		1200.0 		4.527 		34.940 		1488.480
            50 		1300.0 		4.394 		34.949 		1489.610
            51 		1400.0 		4.309 		34.955 		1490.940
            52 		1500.0 		4.255 		34.960 		1492.400
            53 		1600.0 		4.229 		34.961 		1493.970
            54 		1800.0 		4.215 		34.963 		1497.285
            55 		2000.0 		4.212 		34.963 		1500.660
            56 		2200.0 		4.220 		34.964 		1504.090
            57 		2400.0 		4.229 		34.965 		1507.540
            58 		2600.0 		4.242 		34.967 		1511.015
            59 		2800.0 		4.257 		34.969 		1514.515
            60 		3000.0 		4.277 		34.971 		1518.045
            60 		3030.0 		4.277 		34.971 		1518.045
            ];
        svp=svp(:,[2 5]);
        if isempty(bath),
            lthick=D;
        else
            lthick=max(bath(:,2));
        end
        Igood=find(svp(:,1)<lthick);
        svp_bot=interp1(svp(:,1),svp(:,2),lthick);
        svp=[svp(Igood,:); lthick svp_bot];
    case 'Ewing_03032_xt5',
        %27/05/103, 01:39, EW03032, 26 58N   86 56W
        plottitle = 'deep-water profile measured by Ewing May 27 2003';
        load ew0303_6.x7d.txt.mat
        
        lthick=D;
        Igood=find(svp(:,1)<lthick);
        if any(svp(:,1)>lthick),
            svp_bot=interp1(svp(:,1),svp(:,2),lthick);
        else
            svp_bot=svp(end,2);
        end
        svp=[svp(Igood,:); lthick svp_bot];
        
    case 'Ewing_0303_6_x7d',
        plottitle = 'shallow-water profile measured by Ewing June 2 2003';
        load ew03032.xt5.txt.mat
        if isempty(bath),
            lthick=D;
        else
            lthick=max(bath(:,2));
        end
        Igood=find(svp(:,1)<lthick);
        if any(svp(:,1)>lthick),
            svp_bot=interp1(svp(:,1),svp(:,2),lthick);
        else
            svp_bot=svp(end,2);
        end
        svp=[svp(Igood,:); lthick svp_bot];
        
    case 'GOM_GDEM_shallow'
        plottitle = 'shallow-water profile in GOM from GDEM (June)' ;
        load -ascii 'GOM_GDEM_shallow' ;    % depth [m], temp [deg. C], sal [pss-78], sv [m/s]
        svp = [ GOM_GDEM_shallow(:,1) GOM_GDEM_shallow(:,4) ] ;
        if isempty(bath),
            lthick=D;
        else
            lthick=max(bath(:,2));
        end
        Igood=find(svp(:,1)<lthick);
        svp_bot=interp1(svp(:,1),svp(:,2),lthick);
        svp=[svp(Igood,:); lthick svp_bot];
        
    case 'SoCalDeepWater'
        plottitle = 'Southern California deep-water profile from Delphine';
        
        %         svp=[
        %             0     1.506459932041710e+03
        %             1.000000000000000e+01     1.506215615133090e+03
        %             2.000000000000000e+01     1.506029008426910e+03
        %             3.000000000000000e+01     1.505489081047150e+03
        %             5.000000000000000e+01     1.502693324812380e+03
        %             7.500000000000000e+01     1.498126118618880e+03
        %             1.000000000000000e+02     1.494334125573250e+03
        %             1.250000000000000e+02     1.491457660220090e+03
        %             1.500000000000000e+02     1.489333256724740e+03
        %             2.000000000000000e+02     1.486808788058460e+03
        %             2.500000000000000e+02     1.485171015926580e+03
        %             3.000000000000000e+02     1.483764059722830e+03
        %             4.000000000000000e+02     1.481702585807810e+03
        %             5.000000000000000e+02     1.480603297792360e+03
        %             6.000000000000000e+02     1.480237620598360e+03
        %             7.000000000000000e+02     1.480453913455390e+03
        %             8.000000000000000e+02     1.480667903286900e+03
        %             9.000000000000000e+02     1.480985603622360e+03
        %             1.000000000000000e+03     1.481498527886060e+03
        %             1.100000000000000e+03     1.482126831353220e+03
        %             1.200000000000000e+03     1.482895106844280e+03
        %             1.300000000000000e+03     1.483667527863580e+03
        %             1.400000000000000e+03     1.484568151160100e+03
        %             1.500000000000000e+03     1.485482375646880e+03
        %             1.750000000000000e+03     1.488013073450540e+03
        %             2.000000000000000e+03     1.491018425145850e+03
        %             2.500000000000000e+03     1.498285418698310e+03
        %             3.000000000000000e+03     1.506181418056490e+03
        %             3.500000000000000e+03     1.514474457919630e+03
        %             4.000000000000000e+03     1.523113068628910e+03
        %             4.500000000000000e+03     1.532168124749200e+03
        %
        %             ];
        %%From levitcus Santa Catalina basin
        svp=1000*[
            
        0    1.5065
        0.0100    1.5062
        0.0200    1.5060
        0.0300    1.5055
        0.0500    1.5027
        0.0750    1.4981
        0.1000    1.4943
        0.1250    1.4915
        0.1500    1.4893
        0.2000    1.4868
        0.2500    1.4852
        0.3000    1.4838
        0.4000    1.4817
        0.5000    1.4806
        0.6000    1.4802
        0.7000    1.4805
        0.8000    1.4807
        0.9000    1.4810
        1.0000    1.4815
        1.1000    1.4821
        1.2000    1.4829
        1.3000    1.4837
        1.4000    1.4846
        1.5000    1.4855
        1.7500    1.4880
        2.0000    1.4910
        2.5000    1.4983
        3.0000    1.5062
        3.5000    1.5145
        4.0000    1.5231
        4.5000    1.5322];
    
    Igood=find(svp(:,1)<D);
    
    svp=svp(Igood,:);
    svp=[svp;D      svp(end,2) ];
    case '2015 Cabo San Lucas Pt Ballena'
        plottitle = '2015 Cabo San Lucas Pt Ballena';
        svp=1000*[ 
            0 1.5293
    0.0038    1.5293
    0.0085    1.5292
    0.0134    1.5279
    0.0170    1.5278
    0.0178    1.5278
    0.0181    1.5279
    0.0203    1.5277
    0.0241    1.5266
    0.0278    1.5242
    0.0310    1.5227
    0.0347    1.5225
    0.0380    1.5219
    0.0411    1.5219
    0.0450    1.5198
    0.0491    1.5197
    0.0531    1.5189
    0.0542    1.5189
    0.0547    1.5189
    0.0584    1.5182
    0.0612    1.5180
    0.0628    1.5168
    0.0671    1.5150
    0.0696    1.5148
    0.0737    1.5147
    0.0772    1.5145
    0.0796    1.5144
    0.0824    1.5139
    0.0849    1.5100
    0.0881    1.5094
    0.0922    1.5083
    0.0955    1.5079
    0.0982    1.5076
    0.0995    1.5074
    0.1004    1.5076];
    Igood=find(svp(:,1)<D);
    
    svp=svp(Igood,:);
    svp=[svp;D      svp(end,2) ];
     case '2014 MURI'
        plottitle = '2014 MURI experiment deep water';
        svp=[
             0     1.504000376986906e+03
     1.000000000000000e+01     1.503592888669960e+03
     2.000000000000000e+01     1.503446020514418e+03
     3.000000000000000e+01     1.502548336888688e+03
     4.000000000000000e+01     1.501410698053703e+03
     5.000000000000000e+01     1.499965899504160e+03
     6.000000000000000e+01     1.498143298650598e+03
     7.000000000000000e+01     1.496151681342939e+03
     8.000000000000000e+01     1.494342520366701e+03
     9.000000000000000e+01     1.492751720027372e+03
     1.000000000000000e+02     1.491317404530005e+03
     1.100000000000000e+02     1.490159012342119e+03
     1.200000000000000e+02     1.489284031735389e+03
     1.300000000000000e+02     1.488641172912715e+03
     1.400000000000000e+02     1.488117750507746e+03
     1.500000000000000e+02     1.487668746879846e+03
     1.600000000000000e+02     1.487305223137082e+03
     1.700000000000000e+02     1.486982533839291e+03
     1.800000000000000e+02     1.486683615993441e+03
     1.900000000000000e+02     1.486373303013471e+03
     2.000000000000000e+02     1.486042194404160e+03
     2.100000000000000e+02     1.485674065000376e+03
     2.200000000000000e+02     1.485299121029334e+03
     2.300000000000000e+02     1.484918911913654e+03
     2.400000000000000e+02     1.484541284115912e+03
     2.500000000000000e+02     1.484167172695238e+03
     2.600000000000000e+02     1.483788348494932e+03
     2.700000000000000e+02     1.483421217046512e+03
     2.800000000000000e+02     1.483081271464257e+03
     2.900000000000000e+02     1.482754773559175e+03
     3.000000000000000e+02     1.482440084708727e+03
     3.100000000000000e+02     1.482152347954758e+03
     3.200000000000000e+02     1.481882795578505e+03
     3.300000000000000e+02     1.481630487889046e+03
     3.400000000000000e+02     1.481419015085499e+03
     3.500000000000000e+02     1.481262539672060e+03
     3.600000000000000e+02     1.481142490163692e+03
     3.700000000000000e+02     1.481000151516373e+03
     3.800000000000000e+02     1.480858635426472e+03
     3.900000000000000e+02     1.480725701665188e+03
     4.000000000000000e+02     1.480582713086117e+03
     4.100000000000000e+02     1.480458936023132e+03
     4.200000000000000e+02     1.480387465078349e+03
     4.300000000000000e+02     1.480339003851698e+03
     4.400000000000000e+02     1.480292663654683e+03
     4.500000000000000e+02     1.480247805756387e+03
     4.600000000000000e+02     1.480179469004288e+03
     4.700000000000000e+02     1.480114847141252e+03
     4.800000000000000e+02     1.480039177180262e+03
     4.900000000000000e+02     1.479969844643505e+03
     5.000000000000000e+02     1.479924661311473e+03
     5.100000000000000e+02     1.479892494346665e+03
     5.200000000000000e+02     1.479864787866658e+03
     5.300000000000000e+02     1.479852444824867e+03
     5.400000000000000e+02     1.479851395403960e+03
     5.500000000000000e+02     1.479847353919183e+03
     5.600000000000000e+02     1.479847816890158e+03
     5.700000000000000e+02     1.479848779229454e+03
     5.800000000000000e+02     1.479854808236047e+03
     5.900000000000000e+02     1.479861450905341e+03
     6.000000000000000e+02     1.479882629674941e+03
     6.100000000000000e+02     1.479912834725299e+03
     6.200000000000000e+02     1.479950280524080e+03
     6.300000000000000e+02     1.479989396897284e+03
     6.400000000000000e+02     1.480027517189992e+03
     6.500000000000000e+02     1.480059303190049e+03
     6.600000000000000e+02     1.480090263254124e+03
     6.700000000000000e+02     1.480120761044756e+03
     6.800000000000000e+02     1.480157839197660e+03
     6.900000000000000e+02     1.480200238255495e+03
     7.000000000000000e+02     1.480245284153261e+03
     7.100000000000000e+02     1.480296630454750e+03
     7.200000000000000e+02     1.480347542713534e+03
     7.300000000000000e+02     1.480394890340508e+03
     7.400000000000000e+02     1.480439888035717e+03
     7.500000000000000e+02     1.480482190035238e+03
     7.600000000000000e+02     1.480521450185071e+03
     7.700000000000000e+02     1.480568913869696e+03
     7.800000000000000e+02     1.480617556685636e+03
     7.900000000000000e+02     1.480664828325335e+03
     8.000000000000000e+02     1.480715219021054e+03
     8.100000000000000e+02     1.480761184047206e+03
     8.200000000000000e+02     1.480801632088232e+03
     8.300000000000000e+02     1.480843539720138e+03
     8.400000000000000e+02     1.480889127118635e+03
     8.500000000000000e+02     1.480938088870934e+03
     8.600000000000000e+02     1.480990594645177e+03
     8.700000000000000e+02     1.481039841275976e+03
     8.800000000000000e+02     1.481088969068564e+03
     8.900000000000000e+02     1.481136348034993e+03
     9.000000000000000e+02     1.481184077000024e+03
     9.100000000000000e+02     1.481236287514703e+03
     9.200000000000000e+02     1.481290199116178e+03
     9.300000000000000e+02     1.481341028748070e+03
     9.400000000000000e+02     1.481393278881989e+03
     9.500000000000000e+02     1.481448499710275e+03
     9.600000000000000e+02     1.481503382747578e+03
     9.700000000000000e+02     1.481564035241433e+03
     9.800000000000000e+02     1.481630395135556e+03
     9.900000000000000e+02     1.481698708640950e+03
     1.000000000000000e+03     1.481765167090634e+03
     1.010000000000000e+03     1.481831770618472e+03
     1.020000000000000e+03     1.481896570923965e+03
     1.030000000000000e+03     1.481959941201393e+03
     1.040000000000000e+03     1.482025427695281e+03
     1.050000000000000e+03     1.482091594112335e+03
     1.060000000000000e+03     1.482159832932341e+03
     1.070000000000000e+03     1.482230447680397e+03
     1.080000000000000e+03     1.482300638851630e+03
     1.090000000000000e+03     1.482370138965526e+03
     1.100000000000000e+03     1.482438479541455e+03
     1.110000000000000e+03     1.482505610371216e+03
     1.120000000000000e+03     1.482573139510924e+03
     1.130000000000000e+03     1.482644496739684e+03
     1.140000000000000e+03     1.482716688619767e+03
     1.150000000000000e+03     1.482788760094057e+03
     1.160000000000000e+03     1.482860273405186e+03
     1.170000000000000e+03     1.482931287702168e+03
     1.180000000000000e+03     1.482999251300589e+03
     1.190000000000000e+03     1.483068053451198e+03
     1.200000000000000e+03     1.483138937272759e+03
     1.210000000000000e+03     1.483209323228866e+03
     1.220000000000000e+03     1.483279701779062e+03
     1.230000000000000e+03     1.483350492549564e+03
     1.240000000000000e+03     1.483421192532012e+03
     1.250000000000000e+03     1.483493295973730e+03
     1.260000000000000e+03     1.483568537974395e+03
     1.270000000000000e+03     1.483646023580235e+03
     1.280000000000000e+03     1.483725767682395e+03
     1.290000000000000e+03     1.483806306659178e+03
     1.300000000000000e+03     1.483888043374035e+03
     1.310000000000000e+03     1.483970685180973e+03
     1.320000000000000e+03     1.484054649398121e+03
     1.330000000000000e+03     1.484140861867403e+03
     1.340000000000000e+03     1.484229585614799e+03
     1.350000000000000e+03     1.484317902975935e+03
     1.360000000000000e+03     1.484405866258233e+03
     1.370000000000000e+03     1.484493891277377e+03
     1.380000000000000e+03     1.484582107035194e+03
     1.390000000000000e+03     1.484668660047335e+03
     1.400000000000000e+03     1.484757673276852e+03
     1.410000000000000e+03     1.484848209992439e+03
     1.420000000000000e+03     1.484938147702855e+03
     1.430000000000000e+03     1.485028351905397e+03
     1.440000000000000e+03     1.485119705325967e+03
     1.450000000000000e+03     1.485209576766148e+03
     1.460000000000000e+03     1.485300261906210e+03
     1.470000000000000e+03     1.485392908353856e+03
     1.480000000000000e+03     1.485487549741931e+03
     1.490000000000000e+03     1.485584819930969e+03
     1.500000000000000e+03     1.485682282151697e+03
     1.510000000000000e+03     1.485779025757831e+03
     1.520000000000000e+03     1.485875707700199e+03
     1.530000000000000e+03     1.485971414156297e+03
     1.540000000000000e+03     1.486065810852939e+03
     1.550000000000000e+03     1.486160872748129e+03
     1.560000000000000e+03     1.486256995727567e+03
     1.570000000000000e+03     1.486352603341353e+03
     1.580000000000000e+03     1.486448514831198e+03
     1.590000000000000e+03     1.486544187874186e+03
     1.600000000000000e+03     1.486638789322977e+03
     1.610000000000000e+03     1.486734129883486e+03
     1.620000000000000e+03     1.486829549918251e+03
     1.630000000000000e+03     1.486924798212229e+03
     1.640000000000000e+03     1.487020282952709e+03
     1.650000000000000e+03     1.487117402179588e+03
     1.660000000000000e+03     1.487213814930492e+03
     1.670000000000000e+03     1.487310910231714e+03
     1.680000000000000e+03     1.487409761641452e+03
     1.690000000000000e+03     1.487510348581202e+03
     1.700000000000000e+03     1.487611984243111e+03
     1.710000000000000e+03     1.487715044857878e+03
     1.720000000000000e+03     1.487820472885545e+03
     1.730000000000000e+03     1.487923737482222e+03
     1.740000000000000e+03     1.488027821229767e+03
     1.750000000000000e+03     1.488133014468910e+03
     1.760000000000000e+03     1.488238742361399e+03
     1.770000000000000e+03     1.488344151109477e+03
     1.780000000000000e+03     1.488453346688007e+03
     1.790000000000000e+03     1.488564525426230e+03
     1.800000000000000e+03     1.488678076704282e+03
     1.810000000000000e+03     1.488792912107492e+03
     1.820000000000000e+03     1.488908428012250e+03
     1.830000000000000e+03     1.489024001893095e+03
     1.840000000000000e+03     1.489138905669015e+03
     1.850000000000000e+03     1.489253699864437e+03
     1.860000000000000e+03     1.489369922133690e+03
     1.870000000000000e+03     1.489488268229579e+03
     1.880000000000000e+03     1.489608153292664e+03
     1.890000000000000e+03     1.489729325910218e+03
     1.900000000000000e+03     1.489851575259083e+03
     1.910000000000000e+03     1.489973847885715e+03
     1.920000000000000e+03     1.490096526447680e+03
     1.930000000000000e+03     1.490220183690011e+03
     1.940000000000000e+03     1.490345046520354e+03
     1.950000000000000e+03     1.490470999329352e+03
     1.960000000000000e+03     1.490598288234243e+03
     1.970000000000000e+03     1.490725513555371e+03
     1.980000000000000e+03     1.490854091313889e+03
     1.990000000000000e+03     1.490983437693935e+03
     2.000000000000000e+03     1.491113826495653e+03
     2.010000000000000e+03     1.491244465872706e+03
     2.020000000000000e+03     1.491376157095538e+03
     2.030000000000000e+03     1.491507775727561e+03
     2.040000000000000e+03     1.491639663955572e+03
     2.050000000000000e+03     1.491771468602607e+03
     2.060000000000000e+03     1.491904013115601e+03
     2.070000000000000e+03     1.492037308812790e+03
     2.080000000000000e+03     1.492171339389165e+03
     2.090000000000000e+03     1.492306088539716e+03
     2.100000000000000e+03     1.492441539959433e+03
     2.110000000000000e+03     1.492577677343306e+03
     2.120000000000000e+03     1.492714484386326e+03
     2.130000000000000e+03     1.492851944783481e+03
     2.140000000000000e+03     1.492990042229764e+03
     2.150000000000000e+03     1.493128760420164e+03
     2.160000000000000e+03     1.493268083049670e+03
     2.170000000000000e+03     1.493407993813273e+03
     2.180000000000000e+03     1.493548476405965e+03
     2.190000000000000e+03     1.493689514522733e+03
     2.200000000000000e+03     1.493831091858570e+03
     2.210000000000000e+03     1.493973192108463e+03
     2.220000000000000e+03     1.494115798967406e+03
     2.230000000000000e+03     1.494258896130385e+03
     2.240000000000000e+03     1.494402467292393e+03
     2.250000000000000e+03     1.494546496148420e+03
     2.260000000000000e+03     1.494690966393455e+03
     2.270000000000000e+03     1.494835861722489e+03
     2.280000000000000e+03     1.494981165830511e+03
     2.290000000000000e+03     1.495126862412513e+03
     2.300000000000000e+03     1.495272935163484e+03
     2.310000000000000e+03     1.495419367778414e+03
     2.320000000000000e+03     1.495566143952293e+03
     2.330000000000000e+03     1.495713247380113e+03
     2.340000000000000e+03     1.495860661756862e+03
     2.350000000000000e+03     1.496008370777531e+03
     2.360000000000000e+03     1.496156358137109e+03
     2.370000000000000e+03     1.496304607530588e+03
     2.380000000000000e+03     1.496453102652958e+03
     2.390000000000000e+03     1.496601827199209e+03
     2.400000000000000e+03     1.496750764864330e+03
     2.410000000000000e+03     1.496899899343311e+03
     2.420000000000000e+03     1.497049214331144e+03
     2.430000000000000e+03     1.497198693522819e+03
     2.440000000000000e+03     1.497348320613324e+03
     2.450000000000000e+03     1.497498079297651e+03
     2.460000000000000e+03     1.497647953270789e+03
     2.470000000000000e+03     1.497797926227729e+03
     2.480000000000000e+03     1.497947983352187e+03
     2.490000000000000e+03     1.498098154388083e+03
     2.500000000000000e+03     1.498248520522621e+03
     2.510000000000000e+03     1.498399166815233e+03
     2.520000000000000e+03     1.498550179476087e+03
     2.530000000000000e+03     1.498701643135956e+03
     2.540000000000000e+03     1.498853592765969e+03
     2.550000000000000e+03     1.499006019184600e+03
     2.560000000000000e+03     1.499158911725434e+03
     2.570000000000000e+03     1.499312259722056e+03
     2.580000000000000e+03     1.499466052508049e+03
     2.590000000000000e+03     1.499620279417001e+03
     2.600000000000000e+03     1.499774929782495e+03
     2.610000000000000e+03     1.499929992938117e+03
     2.620000000000000e+03     1.500085458217450e+03
     2.630000000000000e+03     1.500241314954080e+03
     2.640000000000000e+03     1.500397552481591e+03
     2.650000000000000e+03     1.500554160133569e+03
     2.660000000000000e+03     1.500711127243600e+03
     2.670000000000000e+03     1.500868443145266e+03
     2.680000000000000e+03     1.501026097172153e+03
     2.690000000000000e+03     1.501184078657846e+03
     2.700000000000000e+03     1.501342376935929e+03
     2.710000000000000e+03     1.501500981339989e+03
     2.720000000000000e+03     1.501659881203609e+03
     2.730000000000000e+03     1.501819065860375e+03
     2.740000000000000e+03     1.501978524643870e+03
     2.750000000000000e+03     1.502138246887681e+03
     2.760000000000000e+03     1.502298221925391e+03
     2.770000000000000e+03     1.502458439090587e+03
     2.780000000000000e+03     1.502618887716852e+03
     2.790000000000000e+03     1.502779557137771e+03
     2.800000000000000e+03     1.502940436686930e+03
     2.810000000000000e+03     1.503101515697913e+03
     2.820000000000000e+03     1.503262783504304e+03
     2.830000000000000e+03     1.503424229439690e+03
     2.840000000000000e+03     1.503585842837655e+03
     2.850000000000000e+03     1.503747613031783e+03
     2.860000000000000e+03     1.503909529355658e+03
     2.870000000000000e+03     1.504071581142868e+03
     2.880000000000000e+03     1.504233757726995e+03
     2.890000000000000e+03     1.504396048441626e+03
     2.900000000000000e+03     1.504558442620344e+03
     2.910000000000000e+03     1.504720929596735e+03
     2.920000000000000e+03     1.504883498704384e+03
     2.930000000000000e+03     1.505046139276874e+03
     2.940000000000000e+03     1.505208840647792e+03
     2.950000000000000e+03     1.505371592150723e+03
     2.960000000000000e+03     1.505534383119249e+03
     2.970000000000000e+03     1.505697202886957e+03
     2.980000000000000e+03     1.505860041437927e+03
     2.990000000000000e+03     1.506022908400895e+03
     3.000000000000000e+03     1.506185836865059e+03
     3.010000000000000e+03     1.506348862969943e+03
     3.020000000000000e+03     1.506512024725560e+03
     3.030000000000000e+03     1.506675360818535e+03
     3.040000000000000e+03     1.506838889374285e+03
     3.050000000000000e+03     1.507002609536036e+03
     3.060000000000000e+03     1.507166519802752e+03
     3.070000000000000e+03     1.507330618673401e+03
     3.080000000000000e+03     1.507494904646948e+03
     3.090000000000000e+03     1.507659376222354e+03
     3.100000000000000e+03     1.507824031898589e+03
     3.110000000000000e+03     1.507988870174615e+03
     3.120000000000000e+03     1.508153889549398e+03
     3.130000000000000e+03     1.508319088521904e+03
     3.140000000000000e+03     1.508484465591097e+03
     3.150000000000000e+03     1.508650019255943e+03
     3.160000000000000e+03     1.508815748015407e+03
     3.170000000000000e+03     1.508981650368453e+03
     3.180000000000000e+03     1.509147724814048e+03
     3.190000000000000e+03     1.509313969851156e+03
     3.200000000000000e+03     1.509480383978742e+03
     3.210000000000000e+03     1.509646965695772e+03
     3.220000000000000e+03     1.509813713501211e+03
     3.230000000000000e+03     1.509980625894023e+03
     3.240000000000000e+03     1.510147701373174e+03
     3.250000000000000e+03     1.510314938437630e+03
     3.260000000000000e+03     1.510482335586354e+03
     3.270000000000000e+03     1.510649891318313e+03
     3.280000000000000e+03     1.510817604132472e+03
     3.290000000000000e+03     1.510985472527795e+03
     3.300000000000000e+03     1.511153495003248e+03
     3.310000000000000e+03     1.511321670057796e+03
     3.320000000000000e+03     1.511489996190403e+03
     3.330000000000000e+03     1.511658471900037e+03
     3.340000000000000e+03     1.511827095685662e+03
     3.350000000000000e+03     1.511995866046241e+03
     3.360000000000000e+03     1.512164781480741e+03
     3.370000000000000e+03     1.512333840488126e+03
     3.380000000000000e+03     1.512503041567364e+03
     3.390000000000000e+03     1.512672383217417e+03
     3.400000000000000e+03     1.512841863937251e+03
     3.410000000000000e+03     1.513011482225832e+03
     3.420000000000000e+03     1.513181236582124e+03
     3.430000000000000e+03     1.513351125505093e+03
     3.440000000000000e+03     1.513521147493704e+03
     3.450000000000000e+03     1.513691301046922e+03
     3.460000000000000e+03     1.513861584663713e+03
     3.470000000000000e+03     1.514031996843040e+03
     3.480000000000000e+03     1.514202536375537e+03
     3.490000000000000e+03     1.514373210692675e+03
     3.500000000000000e+03     1.514544036800972e+03
     3.510000000000000e+03     1.514715031769298e+03
     3.520000000000000e+03     1.514886212190919e+03
     3.530000000000000e+03     1.515057593649098e+03
     3.540000000000000e+03     1.515229181414859e+03
     3.550000000000000e+03     1.515400971949944e+03
     3.560000000000000e+03     1.515572961422843e+03
     3.570000000000000e+03     1.515745146002045e+03
     3.580000000000000e+03     1.515917521856038e+03
     3.590000000000000e+03     1.516090085153313e+03
     3.600000000000000e+03     1.516262832062358e+03
     3.610000000000000e+03     1.516435758751662e+03
     3.620000000000000e+03     1.516608861389715e+03
     3.630000000000000e+03     1.516782136145006e+03
     3.640000000000000e+03     1.516955579186025e+03
     3.650000000000000e+03     1.517129186681260e+03
     3.660000000000000e+03     1.517302954799201e+03
     3.670000000000000e+03     1.517476879708337e+03
     3.680000000000000e+03     1.517650957577158e+03
     3.690000000000000e+03     1.517825184574151e+03
     3.700000000000000e+03     1.517999556867808e+03
     3.710000000000000e+03     1.518174070626617e+03
     3.720000000000000e+03     1.518348722019067e+03
     3.730000000000000e+03     1.518523507213647e+03
     3.740000000000000e+03     1.518698422378847e+03
     3.750000000000000e+03     1.518873463683156e+03
     3.760000000000000e+03     1.519048627295062e+03
     3.770000000000000e+03     1.519223909383058e+03
     3.780000000000000e+03     1.519399306115629e+03
     3.790000000000000e+03     1.519574813661266e+03
     3.800000000000000e+03     1.519750428188460e+03
     3.810000000000000e+03     1.519926145865697e+03
     3.820000000000000e+03     1.520101962861467e+03
     3.830000000000000e+03     1.520277875344261e+03
     3.840000000000000e+03     1.520453879482567e+03
     3.850000000000000e+03     1.520629971444873e+03
     3.860000000000000e+03     1.520806147399671e+03
     3.870000000000000e+03     1.520982403515448e+03
     3.880000000000000e+03     1.521158735960695e+03
     3.890000000000000e+03     1.521335140903901e+03
     3.900000000000000e+03     1.521511614513552e+03
     3.910000000000000e+03     1.521688152958142e+03
     3.920000000000000e+03     1.521864752406158e+03
     3.930000000000000e+03     1.522041409026089e+03
     3.940000000000000e+03     1.522218118986423e+03
     3.950000000000000e+03     1.522394878455653e+03
     3.960000000000000e+03     1.522571683602266e+03
     3.970000000000000e+03     1.522748530594750e+03
     3.980000000000000e+03     1.522925405346740e+03
     3.990000000000000e+03     1.523102310004365e+03
     4.000000000000000e+03     1.523279257812500e+03];
    
       
            
        
    Igood=find(svp(:,1)<D);
    
    svp=svp(Igood,:);
    svp=[svp;D      svp(end,2) ];
    
    
    
end  %case_ssp
