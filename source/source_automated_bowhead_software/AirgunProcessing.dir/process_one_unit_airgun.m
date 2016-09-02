%%%%%%process_one_unit_airgun.m%%%%%%
%% Given a file name and run parameters, extract signals from file and
%% derive level and duration characteristics...
% function[raw_detections,interval_detections,airgun_shots]=process_one_unit_airgun(fname,short_fname,param,run_options)
%Input:
%   fname: Complete filename, including pathname, to be processed (gsi
%       file, aka Greeneridge format file)
%   short_fname: short filename, without complete pathname or extension.
%       Used to name CFAR JAVA files.  The 'goodNames' output of
%       find_DASAR_dates.m is a typical input...
%   param: structure of automated detection parameters, with a param.energy and param.interval_remove and param.airgun substructure.
%       See TOC_param for details.
%   run_options:  Option flags that include
%       force_energy_detector_run:  If 1, force JAVA energy detector to
%           produce a new 'detsum' file
%       Ncalls_to_sample: number of detections to read into memory for
%           detailed processing at any given time.
%       filtering_stage: 'morph','contour','both','none'
%       debug: includes fields
%
%Output:
%   raw_detections: First stage output from energy detector containing
%       .ctime
%       .duration
%   interval_detection: Second stage output, same structure as
%       'raw_detections'
%  	airgun_shots: Output of extract_airgun_levels.m--see that file for
%       details.  Also includes interval and bearing...

%Revised Sept. 6, 2008
% Further ICI processing Revised April 7, 2010

function [raw_detections,interval_detections,airgun_shots]=process_one_unit_airgun(fname,short_fname,param,run_options)

raw_detections=[];interval_detections=[];airgun_shots=[];
debug_params=run_options.debug;
airgun_shots=[];
if isempty(fname),
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Stage one:  run CFAR cell detector, aka 'energy' detector using JAVA %%
%%%program, appending bearing calculation at end                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_options.force_energy_detector_run==1,
    tic
    run_energy_detector(fname,param.energy);
    disp(sprintf('JAVA program run in %6.2f seconds',toc));
    
end
tic

%%Upload all raw detections from JAVA *detsum files into memory%%
detsum=dir([param.energy.dir_out '/' short_fname '*detsum']);
[data_all,head]=readEnergySummary([param.energy.dir_out '/' detsum.name], Inf);
thet=-ones(1,size(data_all.features,2));  %Set default using number of energy detections
%%%Compute bearing of all detections using readEnergySnips
if isempty(strfind(short_fname,'AURAL'))  %Single-channel AURAL files cannot measure bearing, so set to a constant
    thet=extract_snips_bearing(data_all,run_options,param,short_fname);
end

%%Add bearing information to data_all and update data_all.names
data_all.features=[data_all.features; thet];
Nfeature=size(data_all.features,1);
data_all.names{Nfeature}='bearing';

%Provide an opportunity to plot the bearing, ICI, and frequency range of
%raw detections
if run_options.debug.interval==2
    plot_energy_result(1:length(data_all.ctime),'x');
end



if isempty(data_all.ctime)
    disp('No events detected');
    return
end

tmin=datestr(datenum(1970,1,1,0,0,data_all.ctime(1)));
tmax=datestr(datenum(1970,1,1,0,0,data_all.ctime(end)));

disp(sprintf('Raw data loaded and bearings estimated in %6.2f seconds, min time %s max time %s',toc, tmin, tmax));

%%Create raw_detections output variable..
raw_ctimes=data_all.ctime;
raw_detections.ctime=raw_ctimes;
raw_detections.duration=data_all.npt(1,:)/param.Fs;

%%Check if debug ctimes were detected by energy detector...
Imatch_raw=debug_raw_times;
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Stage 2:  Remove (or keep) regular intervals****
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Prepare debug parameters if desired.
if run_options.debug.interval==1;
    interval_debug.fname=fname;
    interval_debug.names=param.interval_remove.names;
    interval_debug.index=param.interval_remove.debug_index;
    if isfield(param.interval_remove,'debug_time')
        interval_debug.time=param.interval_remove.debug_time;
    end
    
else
    interval_debug=[];
end

%%Identify what features to use when filtering ICI
Iwant=[];
for JJ=1:length(param.interval_remove.names)
    for KK=1:length(data_all.names),
        if strcmp(data_all.names{KK},param.interval_remove.names{JJ})
            Iwant=[Iwant KK];
        end
    end
end

%%Compute ICI (regular intervals) from raw data, using selected features
%%for assistance.
try
ICI=compute_ici_bothways_feature(datenum(1970,1,1,0,0,data_all.ctime),param.interval_remove.ICI_range, ...
    data_all.features(Iwant,:), param.interval_remove.names,...
    param.interval_remove.Ndet,param.interval_remove.Nmiss,...
    param.interval_remove.tol_feature,param.interval_remove.ICItol,interval_debug);
catch
    disp('Check your MATLAB path: ICI_detectors may be missing');
   keyboard 
end
%%Added April 5, 2010:  Sometimes walrus sequences and heavy bowhead whale
%%sequences can have an ICI, but the ICI is inconsistent between calls.
%  Thus here we march through each ICI detection and check whether detections
%  nearby share the same ICI.

Iguns=find(ICI>0);

if run_options.debug.interval==2
    hold on
    plot_energy_result(Iguns,'go');
end
ICI_score=ones(size(ICI));
for I=1:length(Iguns)
    current_time=data_all.ctime(Iguns(I));
    current_ICI=ICI(Iguns(I));
    Itest=find(abs(current_time-data_all.ctime(Iguns))<=0.5*param.interval_remove.Ndet*current_ICI);
    Ipass=0;
    for J=1:2  %Harmonic loop:  checks for possibility that a 10 s ICI may have been assigned a 20 s ICI.
       Ipass=Ipass+length(find(abs(ICI(Iguns(Itest))/J-current_ICI)<=param.interval_remove.ICI_std));
    end
    
    %%Are there enough matching ICIs close to the value of the current ICI?
    if (Ipass-1)<param.interval_remove.Nstd
       ICI_score(Iguns(I))=0; 
    else
        %disp('good');
    end
    
end


%%Assign airgun labels to signals with ICIs that are consistent with other
%%ICI signals...

Iguns=find(ICI_score.*ICI>0);
Igood=setdiff(1:length(ICI),Iguns);

%%Assign output variable interval_detections..
interval_ctimes=raw_ctimes(Igood);
interval_detections.ctime=interval_ctimes;
interval_detections.duration=data_all.npt(1,Igood)/param.Fs;

disp(sprintf('After %6.2f sec interval check %i regular signals, %i out of %i detections remain',toc, ...
    length(Iguns),length(Igood),length(data_all.ctime)));
Igood=debug_interval_times;

if run_options.debug.interval==2
    plot_energy_result(Iguns,'rv');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Calculate airgun characteristics%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_options.short_fname=short_fname;
run_options.calibration_keyword=param.calibration_keyword;
%bufferTime=1;
bufferTime=param.energy.bufferTime;
%fprintf('Buffertime is now %6.2f sec \n',bufferTime);
tic
airgun_shots=extract_airgun_levels(Iguns,data_all,head.Fs, ...
    bufferTime,param.airgun.bandwidth,param.energy.Nfft,param.energy.ovlap, ...
	param.airgun.filter_bandwidth,run_options,fname,param.energy.dir_out);
disp(sprintf('%i Airgun signals processed in in %6.2f seconds, %6.2f per sec',length(Iguns),toc,length(Iguns)/toc));
tic


%%Append additional info about detections
airgun_shots.ICI=ICI(Iguns);
airgun_shots.bearing=thet(Iguns);

%eval(sprintf('!mv *detsum %s/',detsumdir));
%!rm *snips *detsum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Inner functions%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%plot_energy_result.m%%%%%%%%%%
    function plot_energy_result(Irange,ssymbol)
        tabs=datenum(1970,1,1,0,0,data_all.ctime(Irange));
        figure(1)
        subplot(4,1,1)
        plot(tabs,thet(Irange),ssymbol);hold on
        xlabel('Time');ylabel('bearing');
        datetick('x',14)
        grid on
        Iwant=[1 5];
        F=param.energy.f_low:(param.energy.bandwidth/2):250;
        Th=(0:10:360);
        
        %     for JJ=1:length(Iwant);
        %         subplot(3,1,JJ+1)
        %         HH=zeros(length(Th)-1,length(F)-1);
        %         for It=1:(length(Th)-1)
        %             Igood=find(thet>=Th(It)&thet<=Th(It+1));
        %             for If=1:(length(F)-1)
        %                 Igood2=find(data_all.features(Iwant(JJ),Igood)>=F(If)&data_all.features(Iwant(JJ),Igood)<=F(If+1));
        %                 HH(It,If)=log10(length(Igood2));
        %             end
        %         end
        %         imagesc(F,Th,HH);
        %         axis('xy');
        %
        %     end
        subplot(4,1,2)
        plot(tabs,data_all.features(1,Irange),ssymbol);ylabel('Min freq (Hz)');xlabel('Time');hold on
        grid on;datetick('x',14)
        subplot(4,1,3)
        plot(tabs, data_all.features(5,Irange),ssymbol);ylabel('Peak freq (Hz)');xlabel('Time');hold on
        grid on;datetick('x',14)
        if exist('ICI')
            subplot(4,1,4)
            plot(tabs, ICI(Irange),ssymbol);ylabel('ICI(sec)');xlabel('Time');hold on
            grid on;datetick('x',14)
        end
        
        %keyboard
    end

%%%%%%%debug_interval_times.m%%%%%%%%%%
    function Igood_out=debug_interval_times
        
        if ~isempty(interval_debug)
            figure;
            tabs=datenum(1970,1,1,0,0,data_all.ctime);
            
            subplot(2,1,1)
            plot(tabs,thet,'x');hold on;plot(tabs(Iguns),thet(Iguns),'ro')

            datetick('x',14,'keeplimits');
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Local time');
            ylabel('bearing(deg)');
            grid on;
            
             subplot(2,1,2)
            plot(tabs,ICI,'x');

            datetick('x',14,'keeplimits');
            set(gca,'fontweight','bold','fontsize',14);
            xlabel('Local time');
            ylabel('interval (sec)');
            grid on;
            ylim([0 40]);
            keyboard;
            
        end
        
        Igood_out=Igood;
        if debug_params.raw>0&&~isempty(debug_params.ctimes_debug),
            disp('');
            for JJ=1:length(Imatch_raw),
                disp(sprintf('ICI for raw detection %i is %6.2f',Imatch_raw(JJ),ICI(Imatch_raw(JJ))));
                %
            end
            
            [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match]=find_similar_elements(debug_params.ctimes_debug, ...
                data_all.ctime(Igood),debug_params.tol);
            disp(sprintf('%i ctimes present after ICI stripping',length(Iy_match)));
            
            
            if length(Iy_match)==0,
                
                yes=input('Activating debug ctime to see if it survives further processing?');
                
                if ~isempty(yes)
                    %                     ICI_temp=compute_ici_bothways_feature(datenum(1970,1,1,0,0,data_all.ctime),[5 40],data_all.features(Iwant,:), ...
                    %                         param.interval_remove.Ndet,param.interval_remove.Nmiss,...
                    %                         param.interval_remove.tol_feature,param.interval_remove.ICItol,interval_debug);
                    
                    ICI(Imatch_raw(JJ))=-1;  %permit this call to be checked in future
                    Igood_out=find(ICI<0);
                end
                
            end
        end
        
    end

%%%%%%%%%%debug_raw_times.m%%%%%%%%%%
    function Imatch_raw=debug_raw_times
        Imatch_raw=[];
        if debug_params.raw>0&&~isempty(debug_params.ctimes_debug),
 %%%           [~,~,~,Iy_match]=find_similar_elements(debug_params.ctimes_debug, ...
            [dum1,dum2,dum3,Iy_match]=find_similar_elements(debug_params.ctimes_debug, ...
                raw_ctimes,debug_params.tol);
            disp(sprintf('%i ctimes present out of %i in raw energy data',length(Iy_match),length(debug_params.ctimes_debug)));
            Imatch_raw=Iy_match;
            for JJ=1:length(Iy_match),
                disp(sprintf('Raw index %i at time %s',Iy_match(JJ),ctime2str(raw_ctimes(Iy_match(JJ)))));
            end
            %keyboard;
        end
    end
end

