%%extract_optimization_parameters.m%%%%%%%
%  function [param, opt_vec,LB,UB]=extract_parameters(param,opt_vec,action,keyword),
%%DELPHINE, I'm sorry not commented very well.
%% This program does several things:
%% (1) It will take a set of values in opt_vec and assign them as fields
%% to param, if action-'assign'
%%  (2) It will take fields in param and assign them as element in opt_vec
%%  if action='extract'l
%% (3) It also creates upper and lower limits of allowed parameter values.
%%     Thus UB are upper limits to the values in opt_vec that the
%%     optimization cannot exceeed, while LB are lower limits.
%%
%%  The string 'keyword' contains a phrase like 'MorphContourEnergy'.
%%  If the program detects the presence of a particular word in the phrase
%% (e.g. 'Energy') then it will assign certain vallues of the contents of
%% opt_vec into the param.energy fields.
%%
%% The one confusing thing about this particular program is that
%% I find I get better performance if I scale, or 'normalize', the parameters in opt_vec so
%% that the values of opt_vec always lie between 0 and 1. I'll point these
%% steps out below.
function [param, opt_vec,LB,UB]=extract_parameters(param,opt_vec,action,keyword)



LB=[];UB=[];

%%DELPHINE, if you want to restrict an optimization value to a set of
%%discrete possibilities, instead of a continuous range of values, you
%% have to do something similar to what I do here.  FOr example, the
%% value of possible spectrogram overlap percents should only be 0, 25%, 50%
%% 75% and 7/8.  Follow how the ovlap_vec is used below...
ovlap_vec=[0 0.25 .5 .75 7/8];
I=1;

if ~isempty(findstr(keyword,'Morph'))
    LB(I)=5;UB(I)=12;I=I+1;  %min dB level, for thresholding background noise
    LB(I)=1; UB(I)=10;I=I+1;
    LB(I)=1*1000/128; UB(I)=20*1000/128;I=I+1;
    LB(I)=0.01; UB(I)=0.2;I=I+1;

    LB(I)=0.05; UB(I)=5;I=I+1;
    LB(I)=0.05; UB(I)=5;I=I+1;

    %LB(I)=0.1; UB(I)=3;I=I+1;
    %LB(I)=4; UB(I)=1000;I=I+1;

    %LB(I)=5; UB(I)=75;I=I+1;
    %LB(I)=200; UB(I)=400;I=I+1;

    %LB(I)=20; UB(I)=300;I=I+1;
    %LB(I)=10; UB(I)=300;I=I+1;

    %LB(I)=0; UB(I)=20;I=I+1;
    %LB(I)=10; UB(I)=400;I=I+1;

    %LB(I)=30; UB(I)=90;I=I+1;
    %LB(I)=-90; UB(I)=-30;I=I+1;

    %LB(I)=0; UB(I)=0.95;I=I+1;

end

%scaling parameters
switch action
    case 'assign'

        %%DELPHINE, this command assumes that the values of opt_vec lie
        %%between 0 and 1.  These values are converted into 'real'
        %%parameter values using UB and LB.  Thus a value of '0' is
        %%converted into the lower bound (LB) for that parameter.

        opt_vec=LB+opt_vec.*(UB-LB);
        I=1;
        param.eq_frac_update=1;

        %%%DELPHINE, example of how to assign opt_vec values to param
        %%%fields.  These parameters are in the same order as 'UB' and 'LB'
        %%%above.
        if ~isempty(findstr(keyword,'Morph')),
            param.morph.SNRmin=opt_vec(I);I=I+1;  %min dB level, for thresholding background noise
            param.morph.dynamic_range=opt_vec(I);I=I+1;
            param.morph.gap_f=opt_vec(I);I=I+1;
            param.morph.gap_t=opt_vec(I);I=I+1;

             param.morph.duration.min=opt_vec(I)-opt_vec(I+1);
             param.morph.duration.max=opt_vec(I)+opt_vec(I+1);I=I+2;
% 
%             %param.morph.time_band_product.min=9.529000000000001e-01;
%             param.morph.time_band_product.min=opt_vec(I);
%             param.morph.time_band_product.max=opt_vec(I);I=I+2;
% 
%             param.morph.robust_fmin.min=opt_vec(I)-opt_vec(I+1);
%             param.morph.robust_fmin.max=opt_vec(I)+opt_vec(I+1);I=I+2;
% 
%             param.morph.robust_fmax.min=opt_vec(I)-opt_vec(I+1);
%             param.morph.robust_fmax.max=opt_vec(I)+opt_vec(I+1);I=I+2;
% 
%             param.morph.local_bandwidth.min=opt_vec(I)-opt_vec(I+1);
%             param.morph.local_bandwidth.max=opt_vec(I)+opt_vec(I+1);I=I+2;
% 
%             param.morph.Orientation.max=opt_vec(I);I=I+1;
%             param.morph.Orientation.min=opt_vec(I);I=I+1;
% 
%             %param.morph.Centroid_freq.min=param.morph.robust_fmin.min+opt_vec(I);LB(I)=20; UB(I)=400;I=I+1;
%             %param.morph.Centroid_freq.max=param.morph.Centroid_freq.min+opt_vec(I);LB(I)=20; UB(I)=400;I=I+1;
% 
%             param.morph.Eccentricity.min=opt_vec(I);I=I+1;

        end

        if ~isempty(findstr(keyword,'Energy')),

            param.energy.threshold=opt_vec(I);I=I+1;
            param.energy.eq_time=opt_vec(I);I=I+1;
            param.energy.bandwidth=round(opt_vec(I));I=I+1;
            param.energy.MinTime=opt_vec(I);I=I+1;
            param.energy.Nfft=2^round(opt_vec(I));I=I+1;
            param.energy.ovlap=round(param.energy.Nfft*ovlap_vec(round(opt_vec(I))));I=I+1;
            param.energy.f_low=round(opt_vec(I));I=I+1;
            param.energy.MaxTime=opt_vec(I);I=I+1;
        end

        if ~isempty(findstr(keyword,'Interval')),


            param.interval_remove.Ndet=opt_vec(I);I=I+1;
            param.interval_remove.ICItol=opt_vec(I);I=I+1;
        end
        if ~isempty(findstr(keyword,'Median')),
            param.median.on=1;
            param.median.size(1)=opt_vec(I);I=I+1;
            param.median.size(2)=opt_vec(I);I=I+1;

        end

        if ~isempty(findstr(keyword,'Filter')),
            param.filter.on=1;
            % if param.filter.on==1,
            %     disp('Gaussian filter activated');
            %  end
            param.filter.size(1)=opt_vec(I);I=I+1;
            param.filter.size(2)=opt_vec(I);I=I+1;
            param.filter.sigma=opt_vec(I);I=I+1;
        end

        if ~isempty(findstr(keyword,'Feature')),
            param.feature.optvec(1,2)=opt_vec(I);I=I+1;
            param.feature.optvec(2,1)=opt_vec(I);I=I+1;
            param.feature.optvec(2,2)=opt_vec(I);I=I+1;
        end

        %param.median_filter=opt_vec(15);

    case 'extract'
        I=1;

        %%DELPHINE, example of how opt_vec is constructed from a set of
        %%fields in param.  Note order of parameters follows same sequence
        %% as 'UB' and the 'assign' statements for 'Morph'
         if ~isempty(findstr(keyword,'Morph')),
            opt_vec(I)=param.morph.SNRmin;I=I+1;  %min dB level, for thresholding background noise
            opt_vec(I)=param.morph.dynamic_range;I=I+1;
            opt_vec(I)=param.morph.gap_f;I=I+1;
            opt_vec(I)=param.morph.gap_t;I=I+1;

             opt_vec(I)=0.5*(param.morph.duration.max+param.morph.duration.min);I=I+1;
             opt_vec(I)=0.5*(param.morph.duration.max-param.morph.duration.min);I=I+1;
% 
%             opt_vec(I)=param.morph.time_band_product.min;I=I+1;
%             opt_vec(I)=param.morph.time_band_product.max;I=I+1;
% 
%             opt_vec(I)=0.5*(param.morph.robust_fmin.min+param.morph.robust_fmin.max);I=I+1;
%             opt_vec(I)=0.5*(param.morph.robust_fmin.max-param.morph.robust_fmin.min);I=I+1;

%             opt_vec(I)=0.5*(param.morph.robust_fmax.min+param.morph.robust_fmax.max);I=I+1;
%             opt_vec(I)=0.5*(param.morph.robust_fmax.max-param.morph.robust_fmax.min);I=I+1;
%            
%             opt_vec(I)=0.5*(param.morph.local_bandwidth.min+param.morph.local_bandwidth.max);I=I+1;
%             opt_vec(I)=0.5*(param.morph.local_bandwidth.max-param.morph.local_bandwidth.min);I=I+1;
%          
%             opt_vec(I)=param.morph.Orientation.max;I=I+1;
%             opt_vec(I)=param.morph.Orientation.min;I=I+1;
% 
%             opt_vec(I)=param.morph.Eccentricity.min;I=I+1;

        end

       

        if ~isempty(findstr(keyword,'Energy')),


            opt_vec(I)=param.energy.threshold;I=I+1;
            opt_vec(I)=param.energy.eq_time;I=I+1;
            opt_vec(I)=param.energy.bandwidth;I=I+1;
            opt_vec(I)=param.energy.MinTime;I=I+1;
            opt_vec(I)=round(log10(param.energy.Nfft)/log10(2));I=I+1;
            opt_vec(I)=3;I=I+1;
            opt_vec(I)=param.energy.f_low;I=I+1;
            opt_vec(I)=param.energy.MaxTime;I=I+1;

        end

        if ~isempty(findstr(keyword,'Interval')),

            opt_vec(I)=param.interval_remove.Ndet;I=I+1;
            opt_vec(I)=param.interval_remove.ICItol;I=I+1;
        end
        if ~isempty(findstr(keyword,'Median')),
            opt_vec(I)=param.median.size(1);I=I+1;
            opt_vec(I)=param.median.size(2);I=I+1;

        end
        if ~isempty(findstr(keyword,'Filter')),
            %  opt_vec(I)=param.filter.on;I=I+1;
            opt_vec(I)=param.filter.size(1);I=I+1;
            opt_vec(I)=param.filter.size(2);I=I+1;
            opt_vec(I)=param.filter.sigma;I=I+1;
        end
        if ~isempty(findstr(keyword,'Feature')),
            opt_vec(I)=param.feature.optvec(1,2);I=I+1;
            opt_vec(I)=param.feature.optvec(2,1);I=I+1;
            opt_vec(I)=param.feature.optvec(2,2);I=I+1;
        end

        %%DELPHINE, this step normalizes opt_vec
        %% vector so values lie between 0 and 1.
        opt_vec=(opt_vec-LB)./(UB-LB);

end

param.no_energy_run_needed=0;

if isempty(findstr(keyword,'Energy')),
    param.no_energy_run_needed=1;
    disp('No energy run needed');
end

%%DELPHINE, since I have normalized opt_vec,
%% the upper and lower bounds given to the optimization
%% routine are simply 0 and 1.
LB=zeros(size(LB));
UB=ones(size(UB));


end

