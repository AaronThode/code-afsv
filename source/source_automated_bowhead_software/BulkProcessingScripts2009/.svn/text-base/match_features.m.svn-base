%function [Ifit,Ishape_anchor,dt,best_mismatch]=match_features(station_anchor,Igood,station_other,params,Nshape);
%Input:
%   params.names={'ctime','fmin','Centroid','BoundingBox'};
%   params.weight=[0 1 0 0];
%   params.tol=[0 25 0 0 ];
%  params.time_tol=distance/1500;
%   params.debug=0;  %If 1, display intermediate output..

%   station{Istation}.features{feature_type}=[Ncall by Nshape] matrix
%   station{Istation}.SNR{1}=[Ncall by Nshape] matrix of PSD-based SNR
%   station{Istation}.ctimes=best_ctimes;
%   Nshape: maximum number of shapes per detection to test..
%Output:
%    Ifit: [Ncall by 2 matrix] of indicies of station_other that
%           match station_anchor indicies.  If no match value of -1.
%           First column: indicie of station_other.features.
%           Second column:  shape number of station_other.features(Ifit(1));
%    Ishape: shape number of station_anchor that yielded best fit.
%
% dt: vertical vector of arrival_time(anchor)-arrival_time(other).

function [Ifit,Ishape_anchor,dt,best_mismatch]=match_features(station_anchor,Igood,station_other,params,Nshape);

ITIME=1;

Ncalls=length(station_anchor.ctimes);
Nshape_max=min([Nshape size(station_anchor.SNR{1},2)]);

Ifit=-1*ones(length(Igood),2);
dt=Ifit(:,1);
Ishape_anchor=dt;
best_mismatch=dt;
for Icall=1:Ncalls,
    if rem(Icall, 500)==0, disp(sprintf('%6.2f percent complete...',Icall/Ncalls));end
    best_mismatch(Icall)=Inf;  %lowest mismatch for each "other" shape number

    for Ishape=1:Nshape_max,
        anchor_time=station_anchor.features{ITIME}(Icall,Ishape);
        for Ishape_other=1:Nshape_max,

            %Are detections within allowable time delay?
            dt_test=abs(anchor_time-station_other.features{ITIME}(:,Ishape_other));
            Icandidate=find(dt_test<=params.time_tol);

            if isempty(Icandidate)
                continue;
            end

            %Compare feature list
            mismatch_err=zeros(length(Icandidate),1);
            for If=2:length(params.names), %Skip time part, cycle through features

                %%Does feature mismatch exceed tolerance?
                anchor_feature=station_anchor.features{If}(Icall,Ishape);
                test_feature=station_other.features{If}(Icandidate,Ishape_other);
                tol_test=abs(test_feature-anchor_feature)./anchor_feature;

                Igood=find(tol_test<params.tol(If));
                Ibad=setdiff(1:length(Icandidate),Igood);

                mismatch_err(Igood)=mismatch_err(Igood)+params.weight(If)*tol_test(Igood).^2;
                mismatch_err(Ibad)=Inf;

            end
            [low_mismatch,Ibest]=min(mismatch_err);
            if ~isinf(low_mismatch)&&low_mismatch<best_mismatch(Icall),

                Ifit(Icall,2)=Ishape_other;
                Ifit(Icall,1)=Icandidate(Ibest);
                best_mismatch(Icall)=low_mismatch;
                dt(Icall)=dt_test(Ifit(Icall,1));
                Ishape_anchor(Icall)=Ishape;

                if params.debug==1,
                    disp(sprintf('New best match: Call %i: Shape %i, Shape other %i',Icall,Ishape,Ishape_other));

                    for If=2:length(params.names)
                        anchor_feature=station_anchor.features{If}(Icall,Ishape);
                        test_feature=station_other.features{If}(Ifit(Icall,1),Ishape_other);
                        tol_test=abs(test_feature-anchor_feature)./anchor_feature;

                        disp(sprintf(' Feature %10s:\t weight %6.2f: anchor value: %6.2f, best match value: %6.2f err: %6.4f', ...
                            params.names{If},params.weight(If),anchor_feature,test_feature,tol_test));


                    end
                    disp(sprintf('dt: %10.6f',dt(Icall)));
                    disp(' ');
                    pause
                end

            end

        end %Ishape_other

    end %Ishape
end
