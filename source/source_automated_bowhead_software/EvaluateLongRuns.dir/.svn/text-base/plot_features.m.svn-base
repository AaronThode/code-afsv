%%%plot and compare detection features against each other%%%%
% function plot_features(best_calls,I_false_alarms,I_auto_match,I_missed_calls);

function attributes=plot_features(best_calls,I_false_alarms,I_auto_match,I_missed_calls);

features=best_calls.features;


%%%Obtain feature names and secondary index if needed..
[feature_name, feature_index]=get_feature_names(features);

%%Extract feature vector to level desired...
%keyboard;
Nshapes=3;  %Number of shapes per detection to evaluate...
attributes.false=extract_feature_vector(features(I_false_alarms),feature_name,feature_index,Nshapes);

attributes.match=extract_feature_vector(features(I_auto_match),feature_name,feature_index,Nshapes);

attributes.missed=extract_feature_vector(features(I_missed_calls),feature_name,feature_index,Nshapes);



attributes.Nfeatures=length(feature_name);
attributes.feature_name=feature_name;

%%plot comparisons of features...
Nplots=attributes.Nfeatures-1;
colorstr='kgcym';
plotstr='oxs';
for Iplot=1:Nplots,
    for Ishape=1:Nshapes,
        %         figure(Ishape)
        %         ps=[colorstr(Ishape) 'o'];
        %
        %         subplot(Nplots,1,Iplot);
        %         plot(attributes.false{1}(:,Ishape),attributes.false{Iplot+1}(:,Ishape),ps);hold on;
        %         plot(attributes.match{1}(:,Ishape),attributes.match{Iplot+1}(:,Ishape),['r' plotstr(Ishape)]);
        %         grid on
        %         xlabel(feature_name{1});ylabel(feature_name{Iplot+1});

        %Make contour plot
        figure(Ishape);
        subplot(2,1,2);
        contour_2D_histogram(attributes.false{1}(:,Ishape),attributes.false{Iplot+1}(:,Ishape),50,50);
        xlabel(feature_name{1});ylabel(feature_name{Iplot+1});
        title(sprintf('Shape %i, false alarms',Ishape));
        subplot(2,1,1);
        contour_2D_histogram(attributes.match{1}(:,Ishape),attributes.match{Iplot+1}(:,Ishape),50,50);
       
        xlabel(feature_name{1});ylabel(feature_name{Iplot+1});
        title('true calls');
        % keyboard;
    end
    hold off
end

figure
for Iplot=1:(Nplots+1),
    subplot(Nplots+1,2,2*Iplot-1);
    [NN,XX]=hist(attributes.false{Iplot},200);
    bar(XX,NN);
    grid on
    title(feature_name{Iplot});
    ylabel('False detections');

    subplot(Nplots+1,2,2*Iplot);
    hist(attributes.match{Iplot},XX);

    grid on
    title(feature_name{Iplot});
    ylabel('True detections');

    hold off
end

figure
if attributes.Nfeatures>2,
    plot3(attributes.false{1}(:),attributes.false{2}(:),attributes.false{3}(:),'ko', ...
        attributes.match{1}(:),attributes.match{2}(:),attributes.match{3}(:),'ro');grid on
    xlabel(feature_name{1});ylabel(feature_name{2});zlabel(feature_name{3});
end

winnow_chc=menu('How set discrimination?','Manually','Optimization','Skip');
if winnow_chc~=3,
    if winnow_chc==1,
        optvec=manual_select_thresholds(attributes.feature_name);
        optvec=reshape(optvec,2,[])';
    else
        disp('Insert default selection criterica here...');
    end

    %%Create the evaluation function...
    eval_str={'first','all'};
    debug=1;
    for I=1:2,

        disp(sprintf('Using evaluation criteria: %s',eval_str{I}));
        disp('False alarm filter...');
        [attributes.false_filter, ctimes_false,attributes.false_pass_I,attributes.false_fail_I]= ...
            filter_feature_vector(attributes.false,feature_name,best_calls.ctime(I_false_alarms), ...
            optvec,eval_str{I},'and',debug);


        disp('True call filter...');
        [attributes.match_filter, ctimes_match,attributes.match_pass_I,attributes.match_fail_I]= ...
            filter_feature_vector(attributes.match,feature_name,best_calls.ctime(I_auto_match), ...
            optvec,eval_str{I},'and',debug);

        Nlost=length(attributes.match{1})-length(ctimes_match);
        Nmissed_new=length(attributes.missed{1})+Nlost;
        Nfalse_new=length(ctimes_false);
        Nmatch_new=length(ctimes_match);

        disp(sprintf('After discrimination: %i matched calls %i missed calls, %i false alarms',Nmatch_new,Nmissed_new,Nfalse_new));
        disp(sprintf('%i true calls lost, %i false alarms stripped\n\n\n',Nlost, length(attributes.false{1})-Nfalse_new));

        missed_fraction=Nmissed_new/(Nmissed_new+Nmatch_new);
        false_ratio=Nfalse_new/Nmatch_new;
        fit=1*missed_fraction+1*false_ratio;
    end

    keyboard;
    % lost_missed=features(attributes.match_fail_I{1});
    %  figure;
    %
    % for J=1:length(lost_missed),
    %     imshow(lost_missed{J}.Image);
    %     pause;
    % end
    % keyboard;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TEXT TO DEMONSTRATE POLYNOMIAL FITS%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %I_false_alarms =Ix_nomatch;
%         %I_auto_match =Ix_match;
%         %I_manual_match =Iy_match;
%         %I_missed_calls =Iy_nomatch;
%
%         %Compare magnitudes...
%         %         magnitude_calls=10*log10(eps+sum(best_calls .magnitude(1:2,I_auto_match )));
%         %         magnitude_false=10*log10(eps+sum(best_calls .magnitude(1:2,I_false_alarms )));
%         %
%         Ifreq=1:10;
%         magnitude_calls=10*log10(eps+sum(best_calls .magnitude(Ifreq,I_auto_match )));
%         magnitude_false=10*log10(eps+sum(best_calls .magnitude(Ifreq,I_false_alarms )));
%
%
%         figure;
%         subplot(2,1,1);plot(10*log10(eps+sum(magnitude_calls)),'ro');hold on;grid on;ylim([50 140]);title('peak SEL');
%         subplot(2,1,2);plot(10*log10(eps+sum(magnitude_false)),'ko');grid on; ylim([50 140]);
%
%
%         duration_calls=(best_calls .duration(1,I_auto_match ));
%         duration_false=(best_calls .duration(1,I_false_alarms ));
%
%         figure;
%         subplot(2,1,1);plot(duration_calls,'ro');hold on;grid on;title('total duration');
%         subplot(2,1,2);plot(duration_false,'ko');grid on
%
%         freq_calls=(best_calls .median_freq(I_auto_match ));
%         freq_false=(best_calls .median_freq(I_false_alarms ));
%
%         figure;
%         subplot(2,1,1);plot(freq_calls,magnitude_calls,'ro');ylim([40 130]);
%         hold off;grid on;xlabel('median frequency');ylabel('magnitude(dB SEL)');
%
%
%         subplot(2,1,2);plot(freq_false,magnitude_false,'ko');ylim([40 130]);
%         grid on;xlabel('median frequency');ylabel('magnitude(dB SEL)');
%
%
%         %Plot contour polynomial fit
%         for III=1:3, %Contour level..
%             figure;
%
%             subplot(2,2,1);
%             plotcc='.x*sdvphox*sdvph';
%             %for III=1:length(best_calls .p),
%             p_calls=best_calls .p{III}(:,I_auto_match );
%             plot(p_calls(end-1,:),p_calls(end-2,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .p{III}(:,I_false_alarms );
%             plot(p_false(end-1,:),20+p_false(end-2,:),[plotcc(2) 'k']);hold on
%             title('scaled slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('slope (Hz/sec)');ylabel('Curvature');
%             grid on
%
%             subplot(2,2,2);
%             p_calls=best_calls .p{III}(:,I_auto_match );
%             plot(p_calls(end,:),p_calls(end-1,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .p{III}(:,I_false_alarms );
%             plot(p_false(end,:),20+p_false(end-1,:),[plotcc(2) 'k']);hold on
%             title('scaled slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('mean F');ylabel('Slope (Hz/sec)');
%             grid on
%
%             subplot(2,2,3);
%
%             p_calls=best_calls .praw{III}(:,I_auto_match );
%             plot(p_calls(end-1,:),p_calls(end-2,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .praw{III}(:,I_false_alarms );
%             plot(p_false(end-1,:),5000+p_false(end-2,:),[plotcc(2) 'k']);hold on
%             title('raw slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('slope (Hz/sec)');ylabel('Curvature');
%             grid on
%             legend('black: airgun','red: calls')
%
%             subplot(2,2,4);
%
%             p_calls=best_calls .praw{III}(:,I_auto_match );
%             plot(p_calls(end,:),p_calls(end-1,:),[plotcc(1) 'r']);hold on
%
%             p_false=best_calls .praw{III}(:,I_false_alarms );
%             plot(p_false(end,:),5000+p_false(end-1,:),[plotcc(2) 'k']);hold on
%             title('raw slopes');
%             set(gca,'fontweight','bold','fontsize',14);
%             xlabel('mean F');ylabel('Slope (Hz/sec)');
%             grid on
%             legend('black: airgun','red: calls')
%         end