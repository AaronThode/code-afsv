%%%%%%%%build_nnetwork.m%%%%%%%
function [ptest, target_est,historyy,net,threshold,pattern_scaling,pattern_scaling_pca,feature_names]=build_nnetwork(whale,other,miss_fraction,nweights,run_options)

if strcmp(run_options.type,'mse')
    false_value=-1;
else
    false_value=0;
end
pattern_scaling=[];
pattern_scaling_pca=[];
threshold=[];

target_w=[];
pattern_w=[];
history_w=[];

if ~isfield(run_options,'train_chc')
    run_options.train_chc='trainlm';
end
%feature_names=fieldnames(whale{1});

for K=1:size(whale,1)
    for SS=1:size(whale,2)
        if ~isempty(whale{K,SS})
            target_pattern=[1;false_value];  %First row is positive for whale...
            [pattern0,target0,feature_names]=build_inputs(whale{K,SS},target_pattern,run_options.skip_feature);
            pattern_w=[pattern_w pattern0];
            if run_options.save_memory==0
            target_w=[target_w target0];
            history_w=[history_w [[1;K;SS]*ones(1,size(pattern0,2)); whale{K,SS}.ctime_min ;whale{K,SS}.duration]];
            end
        end
    end
end

Nwhale=size(target_w,2);
%pattern=[pattern_w pattern_w];
%target=[target_w target_w];
%pattern=[pattern pattern];
%target=[target target];
pattern=pattern_w;
target=target_w;
historyy=history_w;

pattern_o=[];
target_o=[];
history_o=[];

for K=1:size(other,1)
    for SS=1:size(other,2)
        if ~isempty(other{K,SS})
            target_pattern=[false_value;1];
            [pattern0,target0]=build_inputs(other{K,SS},target_pattern,run_options.skip_feature);
            pattern_o=[pattern_o pattern0];
            if run_options.save_memory==0
            target_o=[target_o target0];
            history_o=[history_o [[0;K;SS]*ones(1,size(pattern0,2)); other{K,SS}.ctime_min; other{K,SS}.duration]];
            end

        end
    end
end

Nother=size(target_o,2);

pattern=[pattern pattern_o];
target=[target target_o];
if run_options.save_memory==0
historyy=[historyy history_o];
end

%pattern_org=pattern;
%target_org=target;

%pattern=[pattern_w pattern_w pattern_w pattern_w pattern_w pattern];
%target=[target_w target_w target_w target_w target_w target];

disp(sprintf('%i whale patterns and %i other patterns',Nwhale,Nother));

Ibad=find(isinf(pattern));
disp(sprintf('%i Inf values in pattern matrix',length(Ibad)));
pattern(Ibad)=30;

%%Delphine, these commands can be used directly%%
%[pattern_norm,pattern_scaling] = mapminmax(pattern);
[pattern_norm,pattern_scaling] = mapstd(pattern);


%%Final test pattern
%ptest=mapminmax('apply',[pattern_w pattern_o],pattern_scaling);
ptest=mapstd('apply',[pattern_w pattern_o],pattern_scaling);
ttest=[target_w target_o];

%[pattern_norm,pattern_scaling] = mapstd(pattern);
%[pn,ps1] = mapstd(p);
if run_options.pca==1,
    [pfinal,pattern_scaling_pca] = processpca(pattern_norm,0.001);
    ptest = processpca('apply',ptest,pattern_scaling_pca);
    disp(sprintf('PCA applied, %i pattern vectors remain',size(pfinal,1)));
else
    pfinal=pattern_norm;

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

%pattern_train=pfinal(:,Irand);
%target_train=target(:,Irand);

%nweights=[3  ];

for KK=1:length(nweights),
    if strcmp(run_options.type,'mse')
        net = newff(pfinal(:,Irand),target(:,Irand),nweights(KK)*[1],{},run_options.train_chc);
    else
        net = newff(pfinal(:,Irand),target(:,Irand),[nweights(KK)],{'tansig','logsig'},run_options.train_chc);

        net.performFcn='cross_entropy';
        %net.layers{2}.transferFcn = 'purelin';
        %net.layers{3}.transferFcn = 'logsig';
    end

    net.trainParam.epochs = 500;
    %net.trainParam.goal=1e-10;
    net = train(net,pfinal(:,Irand),target(:,Irand));
    netsave{KK}=net;


    target_est = sim(net,ptest);

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

        Iwhale=find(ttest(1,:)>0);

        figure
        for I=1:2,
            if I==1,
                indexx=Iwhale;
            else
                indexx=((max(Iwhale)+1):size(ttest,2));
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
        grid on
        subplot(3,1,2);

        plot(X{1,1},Nwhale*(1-cum{1,1}),'ro',X{2,1},Nother*(1-cum{2,1}),'ko');
        set(gca,'fontweight','bold','fontsize',14);
        xlabel('Output layer threshold');
        ylabel('Detectons remaining'); grid on

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
                indexx=((max(Iwhale)+1):size(ttest,2));
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
