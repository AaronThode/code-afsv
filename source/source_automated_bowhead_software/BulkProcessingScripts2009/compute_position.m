%function [locations,confidence_ellipse,Igood]=compute_position(locations,goodFile_loc,goodFile,param,Icase,Isite,run_options)
%%%%Given a 'locations' structure, add a field 'positions' that uses
%%%%geographic information and bearing information to track a call.
%% Inputs:
%%      locations: structure of locations returned by cross_matching and compute_bearing
%%      goodFile: cell array of strings that contain DASAR file names and pathnames that correspond with columns in
%%              locations{1}.ctime
%%      goodFile_loc: string containing letters of DASARS with non-zero weight
%%      param:  parameter structure:
%%              param.localize.min_weight;
%%              param.localize.Ikeep:  an array the length of goodFile, indicating initial weights.  Used to
%%                          eliminate bad DASARs right off the bat.
%%                  plus other fields if debug plotting is desired
%%      Icase: Used to identify year and thus DASAR locations (load_DASAR_coords)
%%      Isite: scalar site number (load_DASAR_coords)
%%      run_options:  parameter structure:
%%              run_options.localization_alg='medianHuber','repeated_median','Huber','Andrews'
%%              run_options.filter_chc-'min_weight','dt';
%%              run_options.plot_locations

function [locations,confidence_ellipse,Igood]=compute_position(locations,goodFile_loc,goodFile,param,Icase,Isite,run_options)

A=[];
B=[];
ANG=[];

%%Set weights to zero for DASARS designated as bad...
w_loc=-1*ones(length(goodFile),1);
for I=1:length(goodFile)
    for J=1:length(goodFile_loc)
        if strcmp(goodFile{I},goodFile_loc{J})
            w_loc(I)=1;
        end
    end
    
end
DASAR_coords=load_DASAR_coords(Icase,Isite,goodFile);
Nunits=length(goodFile);
cwater=1490;
mean_coords=mean(DASAR_coords);
CRITVAL=4.60517; %chi2inv(0.90,2);
%CRITVAL=1; %chi2inv(0.90,2);

%if strcmp(run_options.localization_alg,'HuberCalibratedKappa')
%    dinfo=load('../DASARlocations/dsumryS508/dsumryS508.mat');
%    kappa = sd2kappa(dinfo.dstd);
%end

min_weight=param.localize.min_weight;
Igood=0;
confidence_ellipse=-1*ones(1,length(locations));
for I=1:length(locations)
    if (rem(I,500)==0), disp(I);end
    disp('');
    bearings=locations{I}.bearing;
    w=ones(length(bearings),1);
    
    if length(w)~=length(w_loc)
        disp('WARNING!  localization weight check not consistent: compute_position line 59');
    else
        w=w_loc;
    end
    
    %Remove kappa values of zero (bearing uncertainity so high
    %   Van Mises is uniform distribution
    
    Ikeep=find(~isnan(bearings)&(w>=min_weight)&(locations{I}.kappa>0));
    if isempty(Ikeep)
        outcome='Not enough good bearings';
        
    end
    Ikeep_last=[];
    while ((length(Ikeep)~=length(Ikeep_last))&&~isempty(Ikeep))
        Ikeep_last=Ikeep;
        if strcmp(run_options.localization_alg,'medianHuber')
            [VMm,Ikeepm,outcome]=repeated_median(bearings,Ikeep,DASAR_coords);
            w(setdiff(1:length(bearings),Ikeepm))=0;
            if ~isempty(Ikeepm)
                [VM,Qhat,w0,outcome] = vmmle_r(bearings(Ikeepm),DASAR_coords(Ikeepm,:),'h');
                w(Ikeepm)=w0';
                
            end
            
        elseif strcmp(run_options.localization_alg,'repeated_median')
            [VM,Ikeep,outcome]=repeated_median(bearings,Ikeep,DASAR_coords);
            Qhat=[];
        elseif strcmp(run_options.localization_alg,'Huber')
            if isfield(locations{I},'kappa')
                %Ikeep2=find(locations{I}.kappa(Ikeep)>0);
                %Ikeep=Ikeep(Ikeep2);
                [VM,Qhat,w0,outcome] = vmmle_r(bearings(Ikeep),DASAR_coords(Ikeep,:),'h',locations{I}.kappa(Ikeep));
            else
                [VM,Qhat,w0,outcome] = vmmle_r(bearings(Ikeep),DASAR_coords(Ikeep,:),'h');
            end
            w(Ikeep)=w0';
            if isempty(findstr(outcome,'successful')),
                
                
                
                dt_err=NaN*ones(length(Ikeep),1);
                dt_true=dt_err;
                dt_est=dt_true;
            end
            if strcmp(run_options.filter_chc,'min_weight')
                Ikeep=find(~isnan(bearings)&w>=min_weight&(locations{I}.kappa>0));
                
            elseif strcmp(run_options.filter_chc,'dt_error');
                %%Compute the relative arrival times of the signal at all DASARS,
                %%assuming source is at VM
                dt_true=locations{I}.dt(Ikeep);
                [dt_err,dt_est]=dt_error_check(VM,DASAR_coords(Ikeep,:),dt_true);
                
                dt_true_xcorr=locations{I}.dt_xcorr(Ikeep);
                [dt_err2,dt_est2]=dt_error_check(VM,DASAR_coords(Ikeep,:),dt_true_xcorr);
                
                
                if all(isinf(dt_err))
                    Ikeep=[];
                    locations{I}.position.outcome='anchor station has bearing failure';
                    outcome='anchor station has bearing failure';
                    %keyboard;
                    continue;
                end
                if run_options.plot_locations==1
                    disp('%%%%%%%%%%%%%%%');
                    disp(sprintf('Location %i outcome %s',I,outcome));
                    %disp(sprintf('Time:%s, Position: %s\n',ctime2str(locations{I}.position.ctime),mat2str(VM,8)));
                    disp(sprintf('Indicies:%s, weights: %s',mat2str(Ikeep),mat2str(w0,4)));
                    disp(sprintf('Measured dt: %s\nModeled dt:  %s\nabs diff: %s',mat2str(dt_true,4),mat2str(dt_est,4), ...
                        mat2str(dt_err,4)));
                end
                Ikeep=Ikeep(find(abs(dt_err)<0.75|dt_err==0));
                
            end
        elseif strcmp(run_options.localization_alg,'HuberCalibratedKappa')
            [VM,Qhat,w0,outcome] = vmmle_r(bearings(Ikeep),DASAR_coords(Ikeep,:),'h',kappa(Ikeep));
            w(Ikeep)=w0';
            %Ikeep=find(~isnan(bearings)&w>=min_weight);
            Ikeep=find(~isnan(bearings)&(w>=min_weight)&(locations{I}.kappa>0));
            
        elseif strcmp(run_options.localization_alg,'Andrews')
            [VM,Qhat,w0,outcome] = vmmle_r(bearings(Ikeep),DASAR_coords(Ikeep,:),'a');
            w(Ikeep)=w0';
            Ikeep=find(~isnan(bearings)&(w>=min_weight)&(locations{I}.kappa>0));
        end
    end  %While Ikeep is changing
    
    if ~isempty(findstr(outcome,'successful')),
        
        Igood=Igood+1;
        if ~strcmp(run_options.localization_alg,'repeated_median')
            [AREA,A,B,ANG,Baxis] = ellipsparms(Qhat,CRITVAL,mean_coords,VM);
            confidence_ellipse(I)=AREA;
        end
        
        %%Estimate actual time of call, by protecting back in time from
        %%physically closest DASAR
        dist=sqrt(sum((ones(length(Ikeep),1)*VM-DASAR_coords(Ikeep,:)).^2,2));
        [mindist,Iwant]=min(dist);
        locations{I}.position.ctime=locations{I}.ctime_min(Ikeep(Iwant))-mindist/cwater;
        
        locations{I}.position.outcome=outcome;
        locations{I}.position.Ikeep=Ikeep;
        locations{I}.position.location=VM;
        locations{I}.position.w=w(Ikeep);
        locations{I}.position.center_dist=sqrt(sum((VM-mean_coords).^2));
        %locations{I}.position.dt_err=dt_err;
        
        if ~strcmp(run_options.localization_alg,'repeated_median')
            locations{I}.position.Area=AREA;
            locations{I}.position.major=A;
            locations{I}.position.minor=B;
            locations{I}.position.ellipse_ang=ANG;
            locations{I}.position.Baxis=Baxis;
            locations{I}.position.Qhat=Qhat(:)';
            if run_options.plot_locations==1,
                disp(sprintf('Error ellipse size, %10.4e km^2, major, minor axis: %6.2f %6.2f km\n',AREA/1e6,A/1e3,B/1e3));
            end
        end
        
        locations{I}.position.Nused=sum(~isnan(locations{I}.dt));
        
        
        if run_options.plot_locations==1
            disp('%%%%%%%%%%%%%%%');
            disp(sprintf('Location %i',I));
            disp(sprintf('Time:%s, Position: %s\n',ctime2str(locations{I}.position.ctime),mat2str(VM,8)));
            disp(sprintf('Indicies:%s, weights: %s',mat2str(Ikeep),mat2str(w0,4)));
            % disp(sprintf('Measured dt: %s\nModeled dt:  %s\nAbs diff: %s',mat2str(dt_true,4),mat2str(dt_est,4), ...
            %    mat2str(dt_err,4)));
            %disp(sprintf('\n\nMeasured dt_xcorr: %s\nModeled dt_xcorr:  %s\nAbs diff: %s',mat2str(dt_true_xcorr,4),mat2str(dt_est2,4), ...
            %    mat2str(dt_err2,4)));
            display_automated_crosslink_data(locations,I,param,goodFile);
            
            if strcmp(run_options.localization_alg,'repeated_median')
                plot_location(DASAR_coords,bearings,Ikeep,VM);
            elseif strcmp(run_options.localization_alg,'medianHuber')
                disp(sprintf('bearings into median: %i, survived median: %i',length(Ikeep),length(Ikeepm)));
                plot_location(DASAR_coords,bearings,Ikeep,VMm);title('Median, first round');
                plot_location(DASAR_coords,bearings,Ikeepm,VM,A,B,ANG);
                
            else
                figure;
                plot_location(DASAR_coords,bearings,Ikeep,VM,A,B,ANG);
            end
            title(sprintf('Location %i, %s algorithm',I,run_options.localization_alg));
            
            pause;
            close all;
        end
        
        
        %
    else  %if not succesful localization.
        locations{I}.position.outcome=outcome;
        if run_options.plot_locations==1
            Igood=find(~isnan(bearings));
            display_automated_crosslink_data(locations,I,param,goodFile);
            
            if strcmp(run_options.localization_alg,'medianHuber')
                disp(sprintf('bearings into median: %i, survived median: %i',length(Ikeep),length(Ikeepm)));
                plot_location(DASAR_coords,bearings,Ikeep,VMm);title('Median, first round');
                if ~isempty(Ikeepm)
                    plot_location(DASAR_coords,bearings,Ikeepm);
                end
            else
                plot_location(DASAR_coords,bearings,Igood);
                if ~isempty(Ikeep)
                    plot_location(DASAR_coords,bearings,Ikeep);
                end
            end
            title(sprintf('Location %i, , %s algorithm,outcome: %s',I,run_options.localization_alg,outcome));
            
            pause
            close all;
        end  %plot_locations
    end %if succesful localization
    
end %I: locations



end

function [VM,Ikeep,outcome]=repeated_median(bearings,Ikeep,DASAR_coords);
maxdist=35000;
[xm,ym,ints] = rptdmed(pi/2-bearings(Ikeep)*pi/180,DASAR_coords(Ikeep,:),maxdist);
VM=[xm ym];
Ikeep=Ikeep(ints>0);
if all(ints)==0,
    outcome='fail';
else
    outcome='successful';
end
end
