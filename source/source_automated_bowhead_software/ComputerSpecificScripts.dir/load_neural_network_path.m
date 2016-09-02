function [net_param]=load_neural_network_path(param_chc)
%function [net_param]=load_neural_network_path(param_chc)
%net_param fields are [ndir,name,extra,threshold]
% provide path and filename to neural network for feature filtering...
[s,local_machine]=unix('hostname');
local_machine=deblank(local_machine);
extra{1}=[];
threshold=[];

%%Search for keyword to use updated (two step) neural network..

if ~isempty(findstr(lower(param_chc),'thresholdadjust'))
    disp('Loading manually_set thresholds');
    %threshold=[-0.8 0];
    Istart=findstr(lower(param_chc),'thresholdadjust');
    threshold=sscanf(lower(param_chc(Istart:end)),'thresholdadjust%fn%f');
    disp(sprintf('Threshold is %s',mat2str(threshold)));
end

if ~isempty(findstr(lower(param_chc),'neuralnetupdated'))||~isempty(findstr(lower(param_chc),'dawngrebner'))
    param_chc='Shell09_NeuralNetupdated';
end

    
switch param_chc
    
    %Original neural network used in the 2008 analysis--single network
    %only..
    case {'Shell08_final_parameters','Shell08_March_optimized','Shell09_new_interval_test'}
        if strcmp(local_machine,'macmussel.ucsd.edu')
            ndir='/Volumes/macmussel1/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
        elseif ~isempty(findstr(local_machine,'KatsMacPro')) %KatMacGSI
            ndir='/Users/thode/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
        else ~isempty(findstr(local_machine,'thode')) %KatMacGSI
            ndir='/Users/thode/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
            
        end
      
    case {'Shell09_initial','Shell07_initial'}
        if strcmp(local_machine,'macmussel.ucsd.edu')
            ndir='/Volumes/macmussel1/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
            extra{1}.dir='/Users/Shared/Projects/2009_Arctic_Analysis/TSV_files_Shell09/Shell09_PinnipedAnalysis';
        elseif ~isempty(findstr(local_machine,'KatsMacPro')) %KatMacGSI
            
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
            ndir='/Users/thode/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
            extra{1}.dir='/Users/thode/2009_Arctic_Analysis/TSV_files_Shell09/Shell09_PinnipedAnalysis';
        end
    
         %This was the neural net used for the second run of the 2009 analysis
      
    case{'Shell09_NeuralNetupdated'} 
        if strcmp(local_machine,'macmussel.ucsd.edu')
            ndir='/Volumes/macmussel1/Arctic_2009/Processed/Site_05/Shell09_Site5_NeuraNetJan22/Neural_Networks.dir';
            %name{1}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_NonBioVsBio.mat';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
            name{2}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_WhaleVsPinniped.mat';
        elseif ~isempty(findstr(local_machine,'KatsMacPro')) %KatMacGSI
            ndir='/Users/thode/2009_Arctic_Analysis/Processed2009/Site_05/Shell09_Site5_NeuraNetJan22/Neural_Networks.dir';
            %name{1}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_NonBioVsBio.mat';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
            name{2}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_WhaleVsPinniped.mat';
        else ~isempty(findstr(local_machine,'thode'))
            %ndir='/Users/thode/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
            ndir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/NeuralNets.dir';
            name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
            name{2}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_WhaleVsPinniped.mat';
        end
        
       
%     case{'Shell09_NeuralNetDoublePass'}
%         if strcmp(local_machine,'macmussel.ucsd.edu')
%             ndir='/Volumes/macmussel1/Arctic_2009/Processed/Site_05/Shell09_Site5_NeuraNetJan22/Neural_Networks.dir';
%             name{1}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_NonBioVsBio.mat';
%             name{2}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_WhaleVsPinniped.mat';
%         elseif ~isempty(findstr(local_machine,'KatsMacPro')) %KatMacGSI
%             ndir='/Users/thode/2009_Arctic_Analysis/Processed2009/Site_05/Shell09_Site5_NeuraNetJan22/Neural_Networks.dir';
%             name{1}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_NonBioVsBio.mat';
%             %name{1}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_WhaleVsNotWhale.mat';
%             name{2}='Nnetwork_all_Sites1to5_Shell09_initial_10cells_20090827_20090930_WhaleVsPinniped.mat';
%         elseif ~isempty(findstr(local_machine,'thode'))
%             %ndir='/Users/thode/Arctic_2008/Processed/Site_05/Shell08_Site5_PeakNeuralNetTrain/Neural_Networks.dir';
%             %name{1}='Nnetwork_all_Sites1to5_Shell08_March_optimized_10cells_20080821_20080929.mat';
%             
%         end
        
        
end

net_param.name=name;
net_param.dir=ndir;
net_param.extra=extra;
net_param.threshold=threshold;

