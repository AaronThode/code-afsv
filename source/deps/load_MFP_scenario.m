function [MFP_replica_file_out,range_str_out,default_tilt_out,default_SNR_out]=load_MFP_scenario

%%MFP replicas can be divided into Pekeris models, and more complex models...

if strcmp(getenv('USER'),'thode')
    model_dir='/Users/thode/Projects/Arctic_2010/PropagationModeling.dir/';
elseif strcmp(getenv('USER'),'shabadi');
    model_dir='/Users/shabadi/Documents/Shima-SIO/PropagationModeling.dir/';
elseif strcmp(getenv('USERNAME'),'bonnelju');
    model_dir=fullfile('C:','Users','bonnelju','Desktop','SanDiego','PropagationModeling.dir');
    %model_dir='C:\Users\bonnelju\Desktop\SanDiego\PropagationModeling.dir\';
end

%%%%%%%%Short-Range Estimate%%%%%%%%%%%%%%%%
MFP_replica_file{2}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_10to500Hz.mat');
%MFP_replica_file{2}{1}=[model_dir 'Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_ShortRange_CSDM_Nfft2048_20100820T014239_eigenvector_10to500Hz.mat'];
default_tilt{2}{1}='1.94';
range_str{2}{1}='100:5:5000';
default_SNR{2}{1}=[45.78 76.29 88.5 170.9 210.6 253.3];

MFP_replica_file{2}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_ShortRange_arctic_pekeris1707_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{2}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_ShortRange_arctic_pekeris1707_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{2}{2}='2.41';
range_str{2}{2}=range_str{2}{1};
default_SNR{2}{2}=default_SNR{2}{1};

%%%%%%%%%%%%%Medium range 7 km est%%%%%%%%%%%%%%%%
MFP_replica_file{3}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_10to500Hz.mat');
%MFP_replica_file{3}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_MidRange_CSDM_Nfft2048_20100831T004830_eigenvector_10to500Hz.mat'];
%/Users/thode/Projects/Arctic_2010/MidRangeMFPInversion.dir/Thermocline_11dB_eigenvector_inversion.dir/RunAA
%SSP EOF included, 40 pop, thermocline depth fixed, depth restricted to 55 m, gradient greater than -50 m/s
%receiver depth limited to +/1 m/s.
default_tilt{3}{1}='1.35';
range_str{3}{1}='5000:100:10000';
default_SNR{3}{1}=[112.9 116 119 122.1 164.8 167.8 170.9 174 177 180.1 235 238];  %Original 12 frequencies
%used for inversion
%default_SNR{3}{1}=[116 119 122.1 164.8 167.8 170.9 174 177 180.1 ]; %Plotting results

MFP_replica_file{3}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{3}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{3}{2}='1.585';
range_str{3}{2}=range_str{3}{1};
default_SNR{3}{2}=default_SNR{3}{1};

MFP_replica_file{3}{3}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{3}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_MidRange_arctic_pekeris1638_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{3}{3}='1.585';
range_str{3}{3}=range_str{3}{1};
default_SNR{3}{3}=default_SNR{3}{1};

%%%%%%%%%%%%%%%%%Long range estimate (17 km)%%%%%%%%%%%%%%%%%%%%%%%%
MFP_replica_file{4}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_10to500Hz.mat');
%MFP_replica_file{4}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_FarRange_CSDM_Nfft8192_20100831T105013_10to500Hz.mat'];
default_tilt{4}{1}='1.35';
range_str{4}{1}='5000:100:20000';
%default_SNR=12;
default_SNR{4}{1}=[74.01 78.58 85.45 88.5 90.03 92.32 93.08 95.37 99.95 103 107.6 109.9];

MFP_replica_file{4}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_FarRange_arctic_pekeris1692_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{4}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_FarRange_arctic_pekeris1692_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{4}{2}='1.585';
range_str{4}{2}=range_str{4}{1};
default_SNR{4}{2}=default_SNR{4}{1};

%%%%%%%%%%%%%%%%%35 km range estimate%%%%%%%%%%%%%%%%
MFP_replica_file{5}{1}=fullfile(model_dir,'Depth55m_Fixed.dir','ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105748_10to500Hz.mat');
%MFP_replica_file{5}{1}=[model_dir '/Depth55m_Fixed.dir/ShallowBeaufortWhaleInversionSSP_UltraRange_CSDM_Nfft8192_20100831T105748_10to500Hz.mat'];

default_tilt{5}{1}='3.93';
range_str{5}{1}='30000:100:40000';
range_str{5}{1}='100:100:40000';
default_SNR{5}{1}=[99.1800  102.2000  105.3000  108.3000  111.4000  114.4000  117.5000  120.5000 ...
    123.6000  126.6000  129.7000  132.8000  135.8000  138.9000  141.9000  145.0000 ...
    148.0000  151.1000  154.1000  157.2000  160.2000 ];
%default_SNR=default_SNR(1:(end-1));
% default_SNR=default_SNR((end-3):end)

MFP_replica_file{5}{2}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_UltraRange_arctic_pekeris1553_KRAKEN_flat55m_10to500Hz.mat');
%MFP_replica_file{5}{2}=[model_dir '/Depth_55m_fixed_Pekeris.dir/ShallowBeaufortPekeris_UltraRange_arctic_pekeris1553_KRAKEN_flat55m_10to500Hz.mat'];
default_tilt{5}{2}='3.34';
range_str{5}{2}=range_str{5}{1};
default_SNR{5}{2}=default_SNR{5}{1};


%Averaged Pekeris model...
for I=2:5
    %%Averaged Pekeris model
    MFP_replica_file{I}{3}=fullfile(model_dir,'Depth_55m_fixed_Pekeris.dir','ShallowBeaufortPekeris_arctic_pekeris1673_KRAKEN_flat55m_10to500Hz.mat');
    default_tilt{I}{3}=default_tilt{I}{2};
    range_str{I}{3}=range_str{I}{2};
    default_SNR{I}{3}=default_SNR{I}{2};
    
end

yes=menu('Which inversion?','Blank...','Whale Close Range inversion', ...
    'Whale 7km Range inversion','Whale 15km Range inversion','Whale 35km Range inversion');

MFP_replica_file_out=MFP_replica_file{yes};
default_tilt_out=default_tilt{yes};
range_str_out=range_str{yes};
default_SNR_out=default_SNR{yes};

yes=menu('Which model?','Detailed model','Inverted Pekeris','Average Pekeris','All');

if yes<4
    MFP_replica_file_out=MFP_replica_file_out(yes);
    default_tilt_out=default_tilt_out(yes);
    range_str_out=range_str_out(yes);
    default_SNR_out=default_SNR_out(yes);
    
end



end  %load_MFP_scenario
