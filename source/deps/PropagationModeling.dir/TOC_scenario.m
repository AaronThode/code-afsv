%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%TOC_scenario.m%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interactive menus to select propagation scenario
% Output:
%       savename: string of suggested filename
%       code_chc: name of propagation model, 'KRAKEN','RAM','bellhop'
%       case_bath: string describing bathymetry model
%       case_ssp: string describing sound speed profile
%       case_bottom: string describing bottom type
%       bellhop:  ??
%       ro: ranges in meters for received field.  Can be vector
%       rd: receiver depths in meters
%       D: bottom depth in meters
%       fl,fu,fs:  lower and upper frequency, sampling frequency (Hz)
%       Tmin, Tmax: min and max times to compute for a time series synthesis

function [case_output,savename,code_chc,case_bath,case_ssp,case_bottom,bellhop,ro,sd,rd,D,fl,fu,fs,Tmin,Tmax,tilt]=TOC_scenario(Iscenario, case_output, flag_suppress_io)

if nargin<3
    flag_suppress_io=0;  %Use menus
end

if nargin<2
    case_output='time_series';
end
scenario_list={'CedricArctic55mModel';'ShallowBeaufortCloseWhaleInversionSSP';'ShallowBeaufortMidWhaleInversionSSP';'ShallowBeaufortFarWhaleInversionSSP'; ...
    'ShallowBeaufortUltraWhaleInversionSSP'; 'SoCalDeepWater';'Sitka2010';...
    'Computational_Ocean_Acoustics_Pekeris_demo';'bellhop_test';'ram_test_2km'; ...
    'bellhop_test_5km';'bellhop_slope_test2'; ...
    'shallow_water'; 'shallow_water_Ewing'; ...
    'deep_water_Ewing';'freespace_demo';'pekeris_demo'; ...
    'English Channel optimum frequency';'ShallowBeaufortPekeris';'ShallowBeringPekeris'; ...
    'PuntaPiedraLSI';'PuntaPiedraLSIburied'; ...
    'ShallowBeaufortCloseWhaleInversion';'classic_wedge_test';'MURI 2014';'Cabo_granite'; ...
    };
menu_flag=(nargin<1||isempty(Iscenario)||Iscenario<0);

if menu_flag||~flag_suppress_io
    Iscenario=menu('Select a scenario',scenario_list);
end
case_scenario=scenario_list{Iscenario};
fprintf('Case scenario is %s\n',case_scenario);


prompt1={'Output format (group_velocity,TL_C,TL_I,TL_S,time_series,ray)','Propagation code(KRAKEN, BELLHOP, RAM)','water depth (m)','default sound speed keyword-keep blank to select from menu:', ...
    'default bottom type-keep blank to select from menu:','default number of BELLHOP beams', ...
    'default bathymetry type-keep blank to select from menu:'};
prompt2={'ranges (m)','source depth (m)','receiver depths (m)--use D for water depth', '[lower upper] frequencies (Hz)', ...
    'Sampling frequency (Hz)','minimum time modeled (sec)-use ro for range','maximum time modeled (sec)--use Tmin for min time','Tilt (m offset)'};

dlgTitle1='Parameters for scenario selection';
dlgTitle2='Parameters for scenario geometry';

switch case_scenario
    case 'CedricArctic55mModel'
        def1={case_output,'KRAKEN','D','arcticSSP','sand_Arctic_inversionSSP','0','flat'};
        % def1={case_output,'KRAKEN','D','pekeris','sand_Arctic_inversionSSP','0','flat'};
        %def1={case_output,'KRAKEN','D','arctic_pekeris','pekeris1707','0','flat'};
        [rd,D]=get_VA_2010_depths_and_offset(2010,[],case_scenario);
        switch case_output
            case 'time_series'
                def2={'1000','40','rd','[10 500]','6250','ro/1600','Tmin+2','5'};  %Second item is source depth...
            case 'group_velocity'
                def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','4','5'};
            case 'TL_C'  %For source level modeling...
                def2={'[50:25:7000]','54','0:1:55','[10 500]','6250','zeros(size(ro))','4','5'};
        end
        
    case 'ShallowBeaufortCloseWhaleInversionSSP'
        
        def1={case_output,'KRAKEN','D','arcticSSP','sand_Arctic_inversionSSP','0','flat'};
        % def1={case_output,'KRAKEN','D','pekeris','sand_Arctic_inversionSSP','0','flat'};
        %def1={case_output,'KRAKEN','D','arctic_pekeris','pekeris1707','0','flat'};
        
        [rd,D]=get_VA_2010_depths_and_offset(2010,[],case_scenario);
         D=53;
        D=47;
        D=37;
        
        switch case_output
            case 'time_series'
                def2={'1000','40','rd','[10 500]','6250','ro/1600','Tmin+2','5'};  %Second item is source depth...
            case 'group_velocity'
                def2={'[50:100:40000]','2','0:1:55','[10 500]','6250','zeros(size(ro))','4','5'};
            case 'TL_C'  %For source level modeling...
               def2={'[50:25:20000]','D-0.5','0:1:D','[10 300]','6250','zeros(size(ro))','4','5'};
        end
       
                
    case 'ShallowBeaufortMidWhaleInversionSSP'
        
        [rd,D]=get_VA_2010_depths_and_offset(2010,[],case_scenario);
        
        def1={case_output,'KRAKEN','D','arcticSSP_7kmWhale','sand_Arctic_inversion_7kmWhale','5','flat'};
        % def1={case_output,'KRAKEN','D','pekeris','sand_Arctic_inversion_7kmWhale','5','flat'};
        %def1={case_output,'KRAKEN','D','arctic_pekeris','pekeris1707','5','flat'};
         D=53;
        D=47;
        D=37;
        switch case_output
            case 'time_series'
                def2={'7000','40','rd','[10 500]','6250','ro/1600','Tmin+2','5'};
                
            case 'group_velocity'
                def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','4','5'};
            case 'TL_C'  %For source level modeling...
                def2={'[50:25:20000]','D-0.5','0:1:D','[10 300]','6250','zeros(size(ro))','4','5'};
                
        end
        
    case 'ShallowBeaufortFarWhaleInversionSSP'
        
        def1={case_output,'KRAKEN','D','arcticSSP_15kmWhale','sand_Arctic_inversion_15kmWhale','5','flat'};
        %def1={case_output,'KRAKEN','D','pekeris','sand_Arctic_inversion_15kmWhale','5','flat'};
        %def1={case_output,'KRAKEN','D','arctic_pekeris','pekeris1707','5','flat'};
        [rd,D]=get_VA_2010_depths_and_offset(2010,[],case_scenario);
        
        switch case_output
            case 'time_series'
                def2={'17000','40','rd','[10 500]','6250','ro/1600','Tmin+3','5'};
                
            case 'group_velocity'
                def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','4','5'};
            case 'TL_C'  %For source level modeling...
                def2={'[50:25:7000]','54','0:1:55','[10 500]','6250','zeros(size(ro))','4','5'}; %For source level est
                
        end
        
    case 'ShallowBeaufortUltraWhaleInversionSSP'
        
        def1={case_output,'KRAKEN','D','arcticSSP_35kmWhale','sand_Arctic_inversion_35kmWhale','5','flat'};
        %def1={case_output,'KRAKEN','D','pekeris','sand_Arctic_inversion_35kmWhale','5','flat'};
        %def1={case_output,'KRAKEN','D','arctic_pekeris','pekeris1707','5','flat'};
        [rd,D]=get_VA_2010_depths_and_offset(2010,[],case_scenario);
        
        switch case_output
            case 'time_series'
                def2={'35000','40','rd','[10 500]','6250','ro/1600','Tmin+4','5'};
                
            case 'group_velocity'
                def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','4','5'};
            case 'TL_C'  %For source level modeling...
                def2={'[50:25:7000]','54','0:1:55','[10 500]','6250','zeros(size(ro))','4','5'}; %For source level est
                
        end
    case 'SoCalDeepWater'
        def1={case_output,'BELLHOP','1000','SoCalDeepWater','silt_pekeris','200','flat'};
        def2={'[10:20:20000]','20','10:10:D','[150 150]','2000','(ro)/1600','Tmin+0.001','0'};
        
    case 'Sitka2010'
        def1={case_output,'BELLHOP','1000','Sitka_measured_2010','silt_pekeris','2','flat'};
        def2={'[10:20:20000]','287','10:50:D','[300 300]','2000','(ro)/1600','Tmin+0.001','0'};
        
    case 'MURI 2014'
        def1={case_output,'BELLHOP','4000','2014 MURI','silt_pekeris','3','flat'};
        def2={'[10:20:20000]','330','100:50:D','[4000 4000]','25000','(ro)/1600','Tmin+0.001','0'};
        
        %angles of interest: -5 -2 0.75
        
    case 'PuntaPiedraLSIburied'
        %         def1={'group_velocity','KRAKEN','10','pekeris','sand_sedlayer','0','flat'};
        %         def2={'[10:20:2000]','0.5','11','[50 500]','12500','(ro)/1600','Tmin+0.1'};
        %
        %Highest attenuation to date..
        def1={case_output,'KRAKEN','10','pekeris','sand_sedlayer','0','flat'};
        def2={'[100:100:40000]','0.5','10.5','[100 500]','12500','(ro)/1600','Tmin+0.01','0'};
        
        %def1={'group_velocity','KRAKEN','10','pekeris','silt_pekeris','0','flat'};
        %def2={'[10:20:2000]','0.5','11','[50 500]','12500','(ro)/1600','Tmin+0.1'};
        
        
    case 'Computational_Ocean_Acoustics_Pekeris_demo',  %Automatically produces TL plot of Fig 2.30 in COA.  Lossy bottom
        def1={case_output,'KRAKEN','100','pekeris','COA_pekeris_demo','0','flat'};  %Original textbook example
        def2={'1000*[0.01:0.01:3]','36','46','[20 20]','200','ro/1600','Tmin+2','0'};  %Original textbook example
        
        %Bellhop test
        %def1={case_output,'bellhop','100','pekeris','COA_pekeris_demo','0','flat'};
        %def2={'1000*[0.01:0.01:3]','36','2:2:98','[20 20]','200','ro/1600','Tmin+2','0'};
        
    case 'classic_wedge_test'  %COA p.452
        def1={case_output,'RAM','200','pekeris','wedge_demo','750','wedge_demo'};
        %def2={'100:10:2000','D/2','[D/4 D/2]','[25 75]','200','0*ones(size(ro))','Tmin+2','0'};
        def2={'100:10:15000','180','1:10:600','[25 25]','200','0*ones(size(ro))','Tmin+2','0'};
        
    case 'bellhop_test',
        def1={'group_velocity','bellhop','500+500*sqrt(3)/2','pekeris','sand_kupdemo','750','flat'};
        def2={'500*[1 2]','D/2','[D/4 D/2]','[25 75]','200','zeros(size(ro))','Tmin+3','0'};
        
    case 'ram_test_2km',
        def1={case_output,'RAM','100','pekeris','sand_kupdemo','750','flat'};
        %def2={'100:10:2000','D/2','[D/4 D/2]','[25 75]','200','0*ones(size(ro))','Tmin+2','0'};
        def2={'100:10:10000','D/2','1:2:D','[100 200]','200','0*ones(size(ro))','Tmin+2','0'};
        
    case 'bellhop_test_5km',
        def1={case_output,'bellhop','100','pekeris','sand_kupdemo','750','flat'};
        def2={'[ 2.5 5]*1000','D/2','[D/4 D/2]','[25 75]','200','3*ones(size(ro))','Tmin+3','0'};
        
    case 'bellhop_slope_test2',
        def1={case_output,'bellhop','500+500*sqrt(3)/2','pekeris','sand_kupdemo','750','wedge'};
        def2={'1000:50:12500','D/2','1:10:D','[25 75]','200','0*ones(size(ro))','Tmin+2','0'};
        
    case 'bellhop_slope_test',
        def1={case_output,'bellhop','100','pekeris','sand_kupdemo','750','wedge'};
        def2={'100:50:12500','D/2','[D/4 D/2]','[25 75]','200','0*ones(size(ro))','Tmin+3','0'};
        
    case 'pekeris_demo',
        def1={'group_velocity','KRAKEN','10','pekeris','sand_kupdemo','0','flat'};
        def2={'[50:10:5000]','5','[1:1:(D-1)]','[10 500]','4000','(ro)/1600','Tmin+0.1','0'};
        
    case 'shallow_water',
        def1={'group_velocity','bellhop','52.3','GOM_GDEM_shallow','sand_shallow','750','slope'};
        def2={'[1000 2000:2000:10000]','18','7.5','[150 200]','400','ro/1600','Tmin./cos(75*pi/180)','0'};
        
    case 'shallow_water_Ewing',
        def1={'group_velocity','bellhop','30','GOM_GDEM_shallow','sand_shallow','750','Ewing_shallow_slope'};
        def2={'[1000 2000:2000:10000]','7.5','18','[20 200]','400','ro/1600','Tmin./cos(75*pi/180)','0'};
    case 'deep_water_Ewing',
        def1={'group_velocity','KRAKEN','3030','GOM_GDEM_Aaron','silt_pekeris','750','flat'};
        def2={'[250:250:10000]','7.5','[20 330 500]','[10 200]','1500','ro/1600','Tmin./cos(60*pi/180)','0'};
        
    case 'freespace_demo',
        def1={'group_velocity','KRAKEN','3000','pekeris','sand_kupdemo','0','flat'};
        def2={'500','1500','1500','[25 75]','200','0','Tmin+15','0'};
        
    case 'English Channel optimum frequency',
        def1={'group_velocity','KRAKEN','125','pekeris','sand_English_channel','0','flat'};
        def2={'[50:100:80000]','50','59','[50 1600]','4000','(ro)/1600','Tmin+0.1','0'};
        
    case 'ShallowBeaufortPekeris'
        
        %'group velocity' settings
        def1={'group_velocity','KRAKEN','D','arctic','pekeris1638','0','flat'};
        def2={'[50:100:40000]','2','0:1:53','[10 300]','6250','zeros(size(ro))','2','0'};
        
        [rd,D]=get_VA_2010_depths_and_offset(2010,[],case_scenario);
        D=53;
        %D=47;
        D=37;
      
        def1={'group_velocity','KRAKEN',num2str(D),'arctic','pekeris1707','0','flat'};
        def2={'17000','40','rd','[75 125]','6250','ro/1450','Tmin+4','0'};
        
         switch case_output
            case 'time_series'
                def2={'1000','40','rd','[10 500]','6250','ro/1600','Tmin+2','5'};  %Second item is source depth...
            case 'group_velocity'
                def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','4','5'};
            case 'TL_C'  %For source level modeling...
                def1={'TL_C','KRAKEN','D','pekeris','pekeris1707','0','flat'};
                def2={'[50:25:8000]','D-0.5','0:1:D','[10 300]','6250','zeros(size(ro))','4','5'};
        end
    case 'ShallowBeringPekeris'
        %For Bering Sea
        def1={'time_series','KRAKEN','51','pekeris','pekeris1650','0','flat'};
        def2={'9000','44','2:10:D','[30 250]','10000','3','Tmin+4','0'};
        
    case 'Cabo_granite',
        
        def1={'group_velocity','KRAKENC','100','2015 Cabo San Lucas Pt Ballena','granite_gravellayer_Cabo','0','flat'};
        def2={'[50:100:40000]','10','0:0.25:100','[50 1000]','2000','zeros(size(ro))','2','0'};
        
        
        
    case 'ShallowBeaufortCloseWhaleInversion',
        
        def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','2','0'};
        
        
        [rd,D]=get_VA_2010_depths_and_offset(2010);
        def1={'group_velocity','KRAKEN','D','arctic','sand_Arctic_inversion','0','flat'};
        
        
    case 'PuntaPiedraLSI'
        def1={case_output,'KRAKEN','10','pekeris','sand_pekeris','0','flat'};
        def1={case_output,'bellhop','30','pekeris','sand_pekeris','0','slope3.0km20m'};
        
        
        %def2={'[50:100:40000]','2','0:0.25:55','[10 500]','6250','zeros(size(ro))','4','5'};
        def2={'[100:20:2000]','0.5','[ 5 D-0.2]','[100 500]','12500','ro/1600','Tmin+0.01','0'};  %Shallow panga
        def2={'[100:20:2000]','0.55','1:1:D','[100 500]','12500','ro/1600','Tmin+0.01','0'};  %Shallow panga, single frequency
        def2={'[100:20:3000]','5','1:1:D','[100 500]','12500','ro/1600','Tmin+0.01','0'};  %Shallow whale
        
        
        
end  % switch case_scenario
if menu_flag||~flag_suppress_io
    answer=inputdlg(prompt1,dlgTitle1,1,def1);
else
    answer=def1;
end
case_output=answer{1};
code_chc=answer{2};
D=eval(answer{3});  %depth of water [m]

if ~isempty(deblank(answer{4})),
    case_ssp=answer{4};  %'Ewing_0303_6_x7d','Ewing_03032_xt5', pekeris, pekeris_kupdemo, GOM_GDEM_shallow
else
    case_ssp='select_from_menu';
end

if ~isempty(deblank(answer{5})),
    case_bottom=answer{5};
else
    case_bottom='select_from_menu';
end

bellhop.nbeams=str2num(answer{6});  %Number of rays to use with bellhop
if ~isempty(deblank(answer{7})),
    case_bath=answer{7};
else
    case_bath='select_from_menu';
end

if menu_flag||~flag_suppress_io
    answer2=inputdlg(prompt2,dlgTitle2,1,def2);
else
    answer2=def2;
end
ro = eval(answer2{1});   %Plot distances

sd =eval(answer2{2}); %source depth 8 for airguns --> using buoy hydrophone depth (reciprocity) [m]
rd=eval(answer2{3});
if strcmpi(code_chc,'ram')&&length(rd)==1
    errordlg('When using RAM only you need to specify at least two receiver depths--sorry');
    %stop
end
ff=eval(answer2{4});
if strcmpi(code_chc,'ram')&&length(ff)<2
    errordlg('When using RAM need two frequencies minimum');
    stop
end
fl=ff(1);           %lower frequency limit for Kuperman problem
fu=ff(2);           %  upper frequency limit for Kuperman problem
fs=eval(answer2{5});  %Sampling rate
%Tmin=max(0,max(ro)/1600); %Minimum time,sec , changed from 1700 on
Tmin=eval(answer2{6});
if length(Tmin)~=length(ro),
    error('Your vector for Tmin must be same length as ro vector--try using ones(size(ro)) in dialog box'   );
end
%Tmax=max(ro)/1300; % Maximum time, sec
Tmax=eval(answer2{7});  %Based on maximum expected arrival angle

disp(sprintf('Tmin: %6.2f Tmax: %6.2f',Tmin,Tmax));
tilt=eval(answer2{8});

prompt3={'Suggested file name for output'};
dlgTitle3='Parameters for final file output';
def3={sprintf('%s_%s_%s_%s_%s%im_%ito%iHz_%s',case_scenario,case_ssp,case_bottom,code_chc,case_bath,D,fl,fu,case_output)};
if menu_flag||~flag_suppress_io
    answer=inputdlg(prompt3,dlgTitle3,1,def3);
else
    answer=def3;
end
savename=answer{1};
