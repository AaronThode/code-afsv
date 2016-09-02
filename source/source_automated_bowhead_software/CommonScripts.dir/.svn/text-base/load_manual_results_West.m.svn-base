%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load_manual_results_West.m  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function localized_out=load_manual_results_West(localized, individual, manualdir,goodDASAR,date_range,max_ctime,run_options)),
%%  Aaron Thode
%%  Uses WEST results for making list of manual detections...
%
%% Input:
%               manualdir: string containing directory of GSI TSV files,
%               full path name
%               goodDASAR: cell matrix [Ndays x 1] each cell itself a cell
%                   array containing name strings like 'D07s4e', e.g., stations
%                   that are recording on that particular day.
%               date_range: matrix of strings of dates that cover desire
%                       time period.  Each row is a different string for a distinct
%                       day.  String format is 'YYYYMMDD'
%               max_ctime:  If only a partial day's results is desired, use
%                       this...
%               run_options:
%                       .min_stations:  Minimum number of stations to count
%                           as a detection...
%                       .success_only:  include a manual result only if
%                       succesfully localized
%                       .west_results: include only if WEST has localized
%                           to less than 10 km
%  Output:
%       whale calls that exist on at least (run_options.minStation)
%       stations...
%
%      localized: cell matrix [Ndays x 1], each cell a structure containing
%      the following [Ncall x 1] fields:
%               ctev: [Ncall x1 double]-absolute c-time of a localized call
%               atev: [2860x1 double]
%               utmx: [2860x1 double]
%               utmy: [2860x1 double]
%               wctype: [2860x1 double]
%               area: [2860x1 double]
%               axmajor: [2860x1 double]
%               axminor: [2860x1 double]
%               Baxis: [2860x1 double]
%               Ang: [2860x1 double]
%               Nused: [2860x1 double]-number of stations used in result...
%               Outcome: {2860x1 cell}
%               comment: {2860x1 cell}
%               OperatorID: {2860x1 cell}
%               DateTimeProc: {2860x1 cell}
%
%      individual: cell matrix [Ndays x 1] , each cell conatining a
%               structure
%               ctime: [2860x7 double]
%               wctype: [2860x7 double]  call type
%               sigdb: [2860x7 double]
%               stndb: [2860x7 double]  signal-to-noise ratio
%               flo: [2860x7 double]  lowest frequency
%               fhi: [2860x7 double]
%               duration: [2860x7 double] duration in seconds
%
%               maxSNR: cell matrix [Ndays x1] containing max SNR detected for call
%               across array.
%
%               anchor_station: cell matrix [Ndays x 1] containing station number
%               yielding max SNR.
%       localized_other,individual_other:  same as above, for for signals
%       other than whale calls.

function [localized_in,individual]=load_manual_results_West(localized_in,individual,Istrip,manualdir,goodDASAR,date_range,max_ctime,run_options)

%duplicate_time_tol=run_options.duplicate_time_tol;  %Remove times within this sec of previous time

mydir=pwd;
cd(manualdir);
tsv=load('wc0821S508');
cd('WestEstimates.dir');
west=load('internal_vm2008_site5');
cd(mydir);
for J=1:size(date_range,1),
    Nstations=size(individual{J}.ctime,2);
    Ncalls=size(individual{J}.ctime,1);

    %Restrict data to date in question
    labels=num2str(west.obsid);
    Igood=[];
    for I=1:length(west.ctev),
        if (strcmp(date_range(J,:),labels(I,1:8))==1);
            Igood=[Igood I];
        end
    end

    %Limit to whale calls and songs
%     wctype=floor(west.wctype(Igood));
%     Igood2=find(wctype>0&wctype<=7);
%     Igood=Igood(Igood2);
% 
% 
%     Igood2=find(west.nused(Igood)>=2);
%     Igood=Igood(Igood2);
% 
     west=trim_struct(west,Igood);

    %Igood2=find(west.distok(Igood)==1);
    %Igood=Igood(Igood2);
    
    
    
%     west=trim_struct(west,Istrip);
%     
% 
%    

    local_west=trim_struct(west,Istrip);


    %Replace and update output fields

    localized_in{J}.distok=local_west.distok;
    localized_in{J}.utmx=local_west.VMh(:,1);
    localized_in{J}.utmy=local_west.VMh(:,2);
    localized_in{J}.ctev=local_west.ctev;
    localized_in{J}.atev=ctime2str(local_west.ctev);
    localized_in{J}.area=local_west.area_h;
    localized_in{J}.axmajor=local_west.a_h;
    localized_in{J}.axminor=local_west.b_h;
    localized_in{J}.Baxis=local_west.angax_h*180/pi;
    localized_in{J}.Ang=local_west.ang_h*180/pi;
    localized_in{J}.mindist=local_west.mindist;
    localized_in{J}.Outcome_vmh=local_west.outcome_vmh;
    localized_in{J}.wgt_vmh=local_west.wgt_vmh;

     Igood=[];
    for It=1:length(local_west.did),
        if strcmp(local_west.outcome_vmh{It},'successful')
            Igood=[Igood It];
        end
    end
    
     local_west=trim_struct(local_west,Igood);


    %Convert Huber weight vector into a more convenient form.
    % Now every DASAR has an element in wgt, even if not used.

    localized_in{J}.wgt=individual{J}.wgt;

    for I=1:length(localized_in{J}.ctev),
        Iused=find(localized_in{J}.wgt(I,:)==1);
        if ~isempty(Iused)
            localized_in{J}.wgt(I,Iused)=local_west.wgt_vmh{I};
        end

    end

    if run_options.west_results==1,
        Igood=find(local_west.distok==1);
        localized_in{J}=trim_struct(localized_in{J},Igood);
        individual{J}=trim_struct(individual{J},Igood);
        disp(sprintf('Out of %i calls, %i pass west test, %6.2f percent pass',Ncalls,length(Igood),100*length(Igood)/Ncalls));
    end

end


end

function struct_trim=trim_struct(struct_in,Igood);
names=fieldnames(struct_in);
for I=1:length(names),
    %keyboard;
    switch(names{I})
        case 'Qh'
            struct_trim.Qh=struct_in.Qh(:,:,Igood);
        case {'did','outcome_vmh','wgt_vmh'},
            struct_trim.(names{I})=struct_in.(names{I})(Igood);
        case 'DASARstr'
            struct_trim.(names{I})=struct_in.(names{I});
        otherwise
            try,
                struct_trim.(names{I})=struct_in.(names{I})(Igood,:);
            catch,
                struct_trim.(names{I})=struct_in.(names{I})(:,Igood);

            end
    end
end


end




