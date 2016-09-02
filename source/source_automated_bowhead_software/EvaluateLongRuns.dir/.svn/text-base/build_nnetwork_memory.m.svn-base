%%%%%%%%build_nnetwork.m%%%%%%%
%Input:  whale: [Isite Iday] input cell matrix of features associated with desired category.  whale{I,J}
%                   contains a series of fields that each describe a feature.  The field value can be a scalar
%                   or a vector
%        Nwhale: number of training examples
%        other: [Nfeature Nother] input matrix of features associated with false examples
%        miss_fraction: target miss fraction (missed true examples/ total true examples), or 1-"recall"
%        nweights: number of hidden weights to use in multilayer-perceptron
%        run_options: structure with the following fields:
%                       'type':  backpropagation performance criteria.  Default is mean-square-error 'mse'
%                       'train_chc': backpropagation training scheme.  Default is 'trainlm', batch-mode
%                       levenberg-marquardt
%                       'skip_feature':{'duration_debug','ctime','ctime_min','ctime_debug','miss'}--fields in whale{I,J} to skip
%                       'save_memory':  if one, don't write a history showing the origin of each pattern vector in terms of the whale cell matrix
%                       'pca': conduct principal component analysis on the inputs
%                       'debug_plot': plot debug information
%  Outputs:
%        pfinal:  [Nfeature Nexamples] matrix of training data where the first Nwhale columns are the desired
%               target, and the remainder are false examples.  Each row (feature) has zero mean and std=1;
%        target_est: trained network output for each pattern vector
%       
%        historry: [3 Nexamples]  maps feature vectors to original location in whale{} input cell matrix
%        net: neural network structure
%        threshold: scalar that acheives desired miss_fraction
%        pattern_scaling:  output of mapstd;  the command ptest=mapstd('apply',pattern,pattern_scaling)
%           converts the pattern into ptest, which has a mean of zero and standard deviation of one for each row (feature)
%        pattern_scaling_pca: output of processpca, which finds PCA of normalized pattern data.  Used as follows ptest = processpca('apply',ptest,pattern_scaling_pca);
%               Note: this transformation should only be applied after mapstd!
%        feature_names: cell matrix of features used by the network.  Useful to confirm input feature vector
%               constructed correctly
%[pattern,target,historry,net,threshold,preprocess_info,preprocess_info_pca,feature_names
%
%
%  Also check out %NNSAMPLE Sample Training Session

function [pfinal, target_est,historyy,net,threshold,pattern_scaling,pattern_scaling_pca,feature_names]=build_nnetwork_memory(whale,Nwhale,other,Nother,miss_fraction,nweights,run_options)

if strcmp(run_options.type,'mse')
    false_value=-1;
else
    false_value=0;
end
pattern_scaling=[];
pattern_scaling_pca=[];
threshold=[];
historyy=[];

%target_w=[];
%pattern_w=[];
%history_w=[];

if ~isfield(run_options,'train_chc')
    run_options.train_chc='trainlm';
end
%feature_names=fieldnames(whale{1});
preallocate=0;
Iindex=0;
for K=1:size(whale,1)
    for SS=1:size(whale,2)
        if ~isempty(whale{K,SS})
            target_pattern=[1;false_value];  %First row is positive for whale...
            [pattern0,target0,feature_names]=build_inputs(whale{K,SS},target_pattern,run_options.skip_feature);
            if preallocate==0
               pattern=zeros(size(pattern0,1),Nwhale+Nother);
               target=zeros(2,Nwhale+Nother);
               if run_options.save_memory==0
                   historry=zeros(3,Nwhale+Nother);
               end
               preallocate=1;
            end
            %pattern_w=[pattern_w pattern0];
            %target_w=[target_w target0];
            pattern(:,(Iindex+(1:size(pattern0,2))))=pattern0;
            target(:,(Iindex+(1:size(pattern0,2))))=target0;
            Iindex=Iindex+size(pattern0,2);
            
            if run_options.save_memory==0 %history tells you where every pattern feature vector came from in the whale structure
                historyy(:,(Iindex+(1:size(pattern0,2))))=[[1;K;SS]*ones(1,size(pattern0,2)); whale{K,SS}.ctime_min ;whale{K,SS}.duration];
            end
        end
    end
end



pattern_o=[];
target_o=[];
history_o=[];

for K=1:size(other,1)
    for SS=1:size(other,2)
        if ~isempty(other{K,SS})
            target_pattern=[false_value;1];
            [pattern0,target0]=build_inputs(other{K,SS},target_pattern,run_options.skip_feature);
            %pattern_o=[pattern_o pattern0];
            %target_o=[target_o target0];
            
            pattern(:,(Iindex+(1:size(pattern0,2))))=pattern0;
            target(:,(Iindex+(1:size(pattern0,2))))=target0;
            Iindex=Iindex+size(pattern0,2);
            
            if run_options.save_memory==0
                historyy(:,(Iindex+(1:size(pattern0,2))))=[[0;K;SS]*ones(1,size(pattern0,2)); other{K,SS}.ctime_min ;other{K,SS}.duration];
                %history_o=[history_o [[0;K;SS]*ones(1,size(pattern0,2)); other{K,SS}.ctime_min; other{K,SS}.duration]];
            end

        end
    end
end


disp(sprintf('%i whale patterns and %i other patterns',Nwhale,Nother));

Ibad=find(isinf(pattern));
disp(sprintf('%i Inf values in pattern matrix',length(Ibad)));
pattern(Ibad)=30;

%%Delphine, these commands can be used directly%%
%[pattern_norm,pattern_scaling] = mapminmax(pattern);
[pattern_norm,pattern_scaling] = mapstd(pattern);


%%Final test pattern
%ptest=mapminmax('apply',[pattern_w pattern_o],pattern_scaling);
ptest=mapstd('apply',pattern,pattern_scaling);
%ttest=[target_w target_o];

%[pattern_norm,pattern_scaling] = mapstd(pattern);
%[pn,ps1] = mapstd(p);
if run_options.pca==1
    [pfinal,pattern_scaling_pca] = processpca(pattern_norm,0.01);
    ptest = processpca('apply',ptest,pattern_scaling_pca);
    disp(sprintf('PCA applied, %i pattern vectors remain',size(pfinal,1)));
else
    pfinal=ptest;

end
%target_org=target;

% if run_options.debug_plot==1,
%     figure
%     subplot(2,1,1);
%     imagesc(pattern_norm);
%     title('normalized pattern');
%     subplot(2,1,2);
%     imagesc(pfinal);
% end

%randomize
Nsamples=size(pattern,2);
Irand=randperm(Nsamples);

clear pattern_norm other whale

%Prepare to initialize neural network with a small amount of data...

for KK=1:length(nweights)
    if strcmp(run_options.type,'mse')
        %net =
        %newff(pfinal(:,Irand(1:Initial_max)),target(:,Irand(1:Initial_max)), ...
        %  nweights(KK)*[1],{},run_options.train_chc);  %Original mistaken code up to Dec 23, 2011
        
        %%  From email 	Date: 	June 3, 2009 12:42:21 PM PDT from Sankar at Mathworks.
        %%  A mistake one makes in newff is making sure that appropriate
        %%  data ranges are incroporated when initializing network, without
        %%  entering all input data.  Below is the **wrong** way to do it...
        %     net=newff(nn.pattern(:,1:10),nn.target(:,1:10),nn.nhidden,{},train_chc);
        %
        %     %%Adjust network statistics for input ranges...
        %     net.inputs{1}.range(:,1)=min(nn.pattern');
        %     net.inputs{1}.range(:,2)=max(nn.pattern');
        % And here is the correct way:
        inputs2 = [minmax(pfinal) pfinal(:,1)];
        targets2 = [minmax(target) target(:,1)];
        net=newff(inputs2,targets2,nweights(KK)*[1],{},run_options.train_chc);
        
         
    
    else
        %net = newff(pfinal(:,Irand(1:Initial_max)),target(:,Irand(1:Initial_max)),[nweights(KK)],{'tansig','logsig'},run_options.train_chc);

        %net.layers{2}.transferFcn = 'purelin';
        %net.layers{3}.transferFcn = 'logsig';
        
        inputs2 = [minmax(pfinal) pfinal(:,1)];
        targets2 = [minmax(target) target(:,1)];
        net=newff(inputs2,targets2,nweights(KK),{'tansig','logsig'},run_options.train_chc);
        net.performFcn='cross_entropy';
        
    end
    
    %%Adjust network statistics for input ranges...  part of mistaken initialization of newff up to Dec 23,
    %%    2011
    %net.inputs{1}.range(:,1)=min(pfinal(:,Irand)');
    %net.inputs{1}.range(:,2)=max(pfinal(:,Irand)');


    net.trainParam.epochs = 200;
    %net.trainParam.goal=1e-10;
    net = train(net,pfinal(:,Irand),target(:,Irand));
   % netsave{KK}=net;


    target_est = sim(net,pfinal);

    threshold=plot_results;
    if run_options.debug_plot==1,

        plot_weights;
        gtext(sprintf('num hidden units=%i',nweights(KK)));

    end

end


%keyboard;
    function plot_weights

        if run_options.pca==0
            figure
            subplot(2,1,1);
            imagesc(net.Iw{1});
            colormap(flipud(gray));colorbar;
            set(gca,'fontweight','bold','fontsize',14);
            grid on;
            xlabel('Feature Index');
            ylabel('Hidden unit');
            subplot(2,1,2);
            imagesc(net.LW{2}'); set(gca,'fontweight','bold','fontsize',14);
            grid on;
            ylabel('Hidden unit');
            xlabel('Output unit');
            gtext('A negative value in top plot means a high feature value indicates whale signal, and vice versa.')


        else
            figure

            subplot(3,1,1);
            imagesc(pattern_scaling_pca.transform');
            colormap(flipud(gray));colorbar;
            set(gca,'fontweight','bold','fontsize',14);
            grid on;
            ylabel('Original feature index');
            xlabel('PCA component');
            subplot(3,1,2);
            imagesc(net.Iw{1});
            colormap(flipud(gray));colorbar;
            set(gca,'fontweight','bold','fontsize',14);
            grid on;
            xlabel('PCA component');
            ylabel('Hidden unit');
            subplot(3,1,3);
            imagesc(net.LW{2}'); set(gca,'fontweight','bold','fontsize',14);
            grid on;
            ylabel('Hidden unit');
            xlabel('Output unit');
        end

        % keyboard;
    end
%%%%%%%%%%%%%%%%%%%%%%%plot_results.m%%%%%%%%%%%%

    function threshold=plot_results

        Iwhale=find(target(1,:)>0);

        figure
        for I=1:2,
            if I==1,
                indexx=Iwhale;
            else
                indexx=((max(Iwhale)+1):size(target,2));
            end


            for J=1:2,
                [NN{I,J},X{I,J}]=hist(target_est(J,indexx),(-2:0.05:2));
                cum{I,J}=cumsum(NN{I,J})/sum(NN{I,J});

                if run_options.debug_plot==1,
                    subplot(2,2,2*I-2+J)

                    bar(X{I,J},NN{I,J});
                    xlim([-2 2]);

                    title(sprintf('Category %i, weight index %i',I,J));
                end
            end

        end

        %if run_options.debug_plot==1,
        figure;
        subplot(3,1,1)
        plot(X{1,1},cum{1,1},X{2,1},cum{2,1});
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Cumulative percent of target output');
        ylabel('Fraction rejected');
        legend('whale','not whale')
        grid on
        subplot(3,1,2);

        plot(X{1,1},Nwhale*(1-cum{1,1}),'ro',X{2,1},Nother*(1-cum{2,1}),'ko');
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Output layer threshold');
        ylabel('Detectons remaining'); grid on
        legend('whale','not whale')

        subplot(3,1,3);
        false_alarm_ratio=(1-cum{2,1})*Nother./((1-cum{1,1})*Nwhale);
        plot(cum{1,1},false_alarm_ratio,'ko');
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Missed call fraction');
        ylabel('False alarm ratio');grid on;
        ylim([0 2])
        xlim([0 0.5]);
        orient tall
        disp(sprintf('Miss fraction desired is %6.2f',miss_fraction));
        %tmp=ginput(1);
        %miss_ratio=tmp(1);
        %end
        [junk,Iwant]=min(abs(cum{1,1}-miss_fraction));
        threshold=X{1,1}(Iwant);
        disp(sprintf('threshold value is %6.4f',threshold));

        for I=1:2,
            if I==1,
                indexx=Iwhale;
            else
                indexx=((max(Iwhale)+1):size(target,2));
            end
            Ipass{I}=find(target_est(1,indexx)>threshold);
            Npass(I)=length(Ipass{I});
        end

        disp(sprintf('%i manual calls: Before filtering, %i true and %i false alarms',Nwhale,Nwhale,Nother));
        disp(sprintf('%i manual calls: After filtering, %i true and %i false alarms',Nwhale,Npass(1),Npass(2)));


    end

end



function [pattern,target,feature_names]=build_inputs(feature_cell,target_pattern,skip_feature)
names=fieldnames(feature_cell);
Nsamples=length(feature_cell.ctime_min);
Nfeature=length(names);
%pattern=zeros(Nfeature,Nsamples);
%skip_feature={'duration_debug','ctime','ctime_min','ctime_debug','miss'};
J=1;
for I=1:length(names),
    switch lower(names{I})
        case skip_feature
            continue
        otherwise,
            pattern(J,:)=feature_cell.(names{I});
            feature_names{J}=names{I};
            J=J+1;
    end
end
target=target_pattern*ones(1,Nsamples);
end
