  function Ipass=crude_station_feature_filter(stations,param,Idebug)

        filter_names={'duration','Centroid_freq','time_band_product','robust_fmin','local_bandwidth','robust_fmax','Orientation','Eccentricity'};

        if Idebug>1,
            disp(sprintf(' %i candidates to test',length(feature)));
        end
        Ipass=1:length(feature);
        for I=1:length(filter_names),
            Ipass1=find([feature(Ipass).(filter_names{I})]<=param.morph.(filter_names{I}).max & [feature(Ipass).(filter_names{I})]>=param.morph.(filter_names{I}).min);
            if isempty(Ipass1),
                if Idebug>1,
                    disp(sprintf('morph: No candidates pass %s test,: min: %10.6f max: %10.6f', ...
                        filter_names{I},param.morph.(filter_names{I}).min,param.morph.(filter_names{I}).max));
                    disp(sprintf('Actual values for features: %10.6f \n',feature(Ipass).(filter_names{I})));
                    %figure;imagesc(labeled_raw);pause;close;

                end
                Ipass=[];

                return;

            end
            if Idebug>1,
                disp(sprintf('%i candidates pass %s test',length(Ipass1),filter_names{I}));
            end
            Ipass=Ipass(Ipass1);
        end

    end