%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load_manual_results_TSV.m  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  [localized,individual,Istrip,localized_false,individual_false, ...
%    localized_uncertain,individual_uncertain]=load_manual_results_TSV(manualdir,goodDASAR,date_range,max_ctime,run_options,d_coords,Isite)
%%  Aaron Thode
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
%                           succesfully localized
%                       .center_dist_limit:  distance from center of site
%                           in km, beyond which manual calls are rejected.
%                           Only used if .success_only==1
%                       .all_biologics:  
%                           If zero, load whales only as true, pinnipeds as false.
%                           If one, load all whales,pinnipeds and unknown marine mammal as true,
%                                everything else as false.
%                           If two, use 2008 option, load whales only as true, everything else as
%                                   false...
%                       .min_weight:  minimum user-assigned weight to
%                           accept a manual call.
%                       .plot_manual_statistics
%               d_coords:  coordinates of DASARs, used for
%                   restricting distance.
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
%       Istrip: indicies of surviving calls (compared with raw TSV file).
%       localized_false,individual_false:  same as above, for for signals
%       other than whale calls.

function [localized,individual,Istrip,localized_false,individual_false, ...
    localized_uncertain,individual_uncertain]=load_manual_results_TSV(manualdir,goodDASAR,date_range,max_ctime,run_options,d_coords,Isite)

%duplicate_time_tol=run_options.duplicate_time_tol;  %Remove times within this sec of previous time

mydir=pwd;
cd(manualdir);
localized=[];
individual=[];
Istrip=[];
localized_false=[];
individual_false=[];
localized_uncertain=[];
individual_uncertain=[];

for J=1:size(date_range,1)
    
    [site_number,DASARstr,file_templates]=lookup_files(J,date_range,goodDASAR,manualdir);
    
    try
        fnames=dir(file_templates.mat);  %Site number in DASAR_list
        if ~isempty(fnames),
            disp(['Loading ' fnames(1).name]);
            load(fnames(1).name);
            
        else
            fprintf('Searching for %s\n',file_templates.tsv);
            fnames=dir(file_templates.tsv);  %Site number in DASAR_list
            
            disp(['Importing ' fnames(1).name]);
            [individual_stored,localized_stored]=read_tsv(fnames(1).name,0,Inf,goodDASAR{J});
            Idot=findstr('.tsv',fnames(1).name)-1;
            save_name=[fnames(1).name(1:Idot) '.mat'];
            save(save_name,'individual_stored','localized_stored');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%
        %%Trim times to less than requested max_ctime
        %%%%%%%%%%%%%%%%%%%%%
        
        if isempty(individual_stored)
            cd(mydir);
            return;
        else
            Igood=[];
            for I=1:size(individual_stored.ctime,1)
                Ipass=find(individual_stored.ctime(I,:)>0);
                if min(individual_stored.ctime(I,Ipass))<=max_ctime
                    Igood=[Igood I];
                end
            end
            %Igood=find(min(individual_stored.ctime')<=max_ctime);
            individual{J}=trim_struct(individual_stored,Igood);
            localized{J}=trim_struct(localized_stored,Igood);
            
            Istrip=1:length(localized{J}.ctev);
            disp(sprintf('Maximum time loaded is %s',ctime2str(max(max(individual{J}.ctime)))));
        end
    catch
        disp(sprintf('File %s failed to load',['*wc' date_range(J,(end-3):end) '*s' site_number '*.tsv']));
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make sure calls are on minimum of two stations...
    % Originally, derived this information from the 'individual ctime' structure,
    % but learned that analysts can log a call without including in
    % localization..
    
    num_good_stations=localized{J}.Nused;
    
    if length(run_options.min_stations)==1
        Igood=find(num_good_stations>=run_options.min_stations);
        fprintf('Out of %i  manual whale locations, %i used %i or greater DASARS. \n',length(num_good_stations),length(Igood),run_options.min_stations);
        
    else
        Igood=find(num_good_stations>=run_options.min_stations(1)&num_good_stations<=run_options.min_stations(2));
        fprintf('Out of %i  manual whale locations, %i used between %i and %i DASARS. \n',length(num_good_stations),length(Igood), ...
            run_options.min_stations(1),run_options.min_stations(2));
        
        
    end
    localized{J}=trim_struct(localized{J},Igood);
    individual{J}=trim_struct(individual{J},Igood);
    Istrip=Istrip(Igood);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %If succesful localizations only desired, trim..
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if run_options.success_only==1
        restrict_manual_data_by_position;
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Keep whale calls only--asssuming 2008 classification scheme below
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----------------------------------------------------------------------
    % 	wctype		Call Classification
    % 	----------------------------------------------------------------------
    % 	1.0			Up Sweep(/)
    % 	2.0			Down Sweep(\)
    % 	3.0			Near Constant tone(-)
    % 	4.0			Undulation(U)
    % 	5.0			Undulation(^)
    % 	6.0			Undulation(other)
    % 	7.0			Complex Call
    % 	8.0			Bearded Seal
    % 	9.0			Walrus
    % 	10.0		Man-made
    % 	11.0		Uncertain(Marine Mammal)
    % 	12.0		Unknown(Not Whale)
    % 	1.x			Call Sequence: Up Sweep(/)
    % 	2.x			Call Sequence: Down Sweep(\)
    % 	3.x			Call Sequence: Near Constant tone(-)
    % 	4.x			Call Sequence: Undulation(U)
    % 	5.x			Call Sequence: Undulation(^)
    % 	6.x			Call Sequence: Undulation(other)
    % 	7.x			Call Sequence: Complex Call
    % 	8.x			Call Sequence: Bearded Seal
    % 	9.x			Call Sequence: Walrus
    % 	10.x			Call Sequence: Man-made
    % 	11.x			Call Sequence: Uncertain,Marine Mammal
    % 	12.x			Call Sequence: Unknown(Not Whale)
    %
    % 	where x denotes an operator's estimate of the number of calls.
    %
    
    wctype=floor(localized{J}.wctype);
    nsequence=localized{J}.nseq;
   % nseq
    
    if run_options.all_biologics==0
        %Divide between whales and pinnipeds.  Iwhale are wahles and Ifalse
        %are pinnipeds
        %
        
        Iwhale=find(wctype>0&wctype<=7);
        Iuncertain=find(wctype==12);
        Iman=find(wctype==10);  %Man made
        Ifalse=setdiff(1:length(wctype),unique([Iwhale; Iuncertain; Iman]));  %%Ifalse are pinnipeds
        
        fprintf('Out of %i  manual locations, %i are whales, %i are pinnipeds, and %i are uncertain. \n', ...
            length(localized{J}.ctev),length(Iwhale), length(Ifalse),length(Iuncertain));
    elseif run_options.all_biologics==1
        %Divide between biologics as Iwhale (included pinnipeds) and everything else as false
        Iwhale=find(wctype>0&wctype<=9);
        Iwhale=[Iwhale; find(wctype==11)];
        Iuncertain=find(wctype==12);
        Ifalse=setdiff(1:length(wctype),unique([Iwhale; Iuncertain]));
        %Ifalse=find((wctype>9)&(wctype~=11));
        fprintf('Out of %i  manual locations, %i are biologics, %i are not biologics, and %i are uncertain non-biologics. \n', ...
            length(localized{J}.ctev),length(Iwhale), length(Ifalse),length(Iuncertain));
    elseif run_options.all_biologics==2
        %Divide between whales and everything else.  Iwhale are whales, and
        %everything else (including pinnipeds) are false.  2008 option.
        %
        
        Iwhale=find(wctype>0&wctype<=7);
        Iuncertain=find(wctype==11);
        Ifalse=setdiff(1:length(wctype),[Iwhale; Iuncertain]);  %%Ifalse are everything else
        
        fprintf('Out of %i  manual locations, %i are whales, %i are not whale, and %i are uncertain. \n', ...
            length(localized{J}.ctev),length(Iwhale), length(Ifalse),length(Iuncertain));
       
    end
    
    
    Iseq=find(nsequence>0);
    
    fprintf('%i sequences detected \n',length(Iseq));
    %for KK=1:length(Iseq),
    %    disp(sprintf('Sequence %i: %i instances of call type %i',KK,nsequence(Iseq(KK)), wctype(Iseq(KK))));
    %end
    
    Istrip=Istrip(Iwhale);
    
    localized_false{J}=trim_struct(localized{J},Ifalse);
    individual_false{J}=trim_struct(individual{J},Ifalse);
    
    localized_uncertain{J}=trim_struct(localized{J},Iuncertain);
    individual_uncertain{J}=trim_struct(individual{J},Iuncertain);
    
    localized{J}=trim_struct(localized{J},Iwhale);
    individual{J}=trim_struct(individual{J},Iwhale);
    
   
    %Pull out index of maximum SNR
    [individual{J}.maxSNR,individual{J}.anchor_station]=max(individual{J}.stndb');
    
    individual{J}.DASARstr=DASARstr;
    
end

fprintf('Finished loading manual results.\n\n');

cd(mydir);

if run_options.plot_manual_statistics==1,
    
    plot_call_statistics;
    disp('');
end

    function restrict_manual_data_by_position
        
        %%Important!  If DASAR not used in localization, then we aren't
        %%interested...
        
        fprintf('Want succesfull localization only; begin with %i localizations. \n',length(Igood));
        
        %%Remove bearings in detection that are never used in
        %%localization...
        disp('Am keeping all times, even those that are not used in localization');
        % individual{J}.ctime=individual{J}.ctime.*sign(individual{J}.wgt);   
        
        Igood=[];
        for It=1:length(localized{J}.ctev),
            if strcmp(localized{J}.Outcome{It}','successful')
                Igood=[Igood It];
            end
        end
        fprintf('Out of %i remaining manual locations, %i were successfully tracked. \n',It,length(Igood));
        
        localized{J}=trim_struct(localized{J},Igood);
        individual{J}=trim_struct(individual{J},Igood);
        Istrip=Istrip(Igood);
        
        %%%Restrict distances to less than a certain range...
        % center_dist_limit=1000*input('Enter maximum distance from DASAR center in km (default 11 km + 10 km= 21 km)');
        center_dist_limit=run_options.center_dist_limit;
        if isempty(center_dist_limit)
            center_dist_limit=100000;
        end
        fprintf('Localizations less than %i km used. \n',center_dist_limit/1000);
        
        %center_distance=sqrt((localized{J}.utmx-mean_coords(1)).^2+(localized{J}.utmy-mean_coords(2)).^2);
        %center_distance=sqrt((localized{J}.utmx-mean_coords(1)).^2);
        p = [localized{J}.utmx localized{J}.utmy];                     % Whale call locations
        
        if Isite~=1
            p1 = mean(d_coords(1:2,:));  % Get endpoints of a line segment through
            p2 = mean(d_coords(6:7,:));  % Get endpoints of a line segment through
            center_distance = segment_points_dist_2d(p1,p2,p);
            Igood=find(center_distance<=center_dist_limit);
            
        else
            
            p1 = mean(d_coords(1:2,:));  % Get endpoints of a line segment through
            p2 = mean(d_coords(6:7,:));  % Get endpoints of a line segment through
            center_distance1 = segment_points_dist_2d(p1,p2,p);
            
            p1 = mean(d_coords(8:9,:));  % Get endpoints of a line segment through
            p2 = mean(d_coords(10:11,:));  % Get endpoints of a line segment through
            center_distance2 = segment_points_dist_2d(p1,p2,p);
            center_distance=min([center_distance1 center_distance2],[],2);
            Igood=find(center_distance<=center_dist_limit);
            
        end
        
        %   the approximate center of the array.
        
        
        disp(sprintf('Out of %i remaining manual locations, %i were less than %i km from center.', ...
            length(Istrip),length(Igood),center_dist_limit/1000));
        
        localized{J}=trim_struct(localized{J},Igood);
        individual{J}=trim_struct(individual{J},Igood);
        localized{J}.distance=center_distance(Igood);
        Istrip=Istrip(Igood);
        
        %%%finally, remove individual times form consideration, if assigned Huber weight is too low.
        Nlocs_zerowgt=[0 0];
        for It=1:length(localized{J}.ctev)
            Ibad=find(individual{J}.wgt(It,:)<run_options.min_weight);
            
            individual{J}.ctime(It,Ibad)=0;
            individual{J}.wgt(It,Ibad)=0;
            localized{J}.Nused(It)=sum(sign(individual{J}.wgt(It,:)));
            Nlocs_zerowgt(1)=Nlocs_zerowgt(1)+length(Ibad);
            Nlocs_zerowgt(2)=Nlocs_zerowgt(2)+1;
        end
        fprintf('%i manual locations had a total of %i bearings with a min Huber weight below %6.2f',Nlocs_zerowgt(2),Nlocs_zerowgt(1),run_options.min_weight);
    end

    function plot_call_statistics
        
        names={'flo','fhi','duration'};
        unitts={'Hz','Hz','sec'};
        maxval=[500 500 3;
            100 100 1];
        
        for Idate=1:length(individual),
            figure
            clear tmp
            for Iname=1:length(names),
                for I=1:size(individual{Idate}.ctime,1),
                    Igood=find(individual{Idate}.ctime(I,:)>0);
                    tmp.(names{Iname}).mean(I)=mean(individual{Idate}.(names{Iname})(I,Igood));
                    tmp.(names{Iname}).std(I)=std(individual{Idate}.(names{Iname})(I,Igood));
                    
                end
                
                subplot(length(names),2,2*Iname-1)
                [N,F]=hist(tmp.(names{Iname}).mean,linspace(0,maxval(1,Iname),250));
                if Idate==1,
                    Ntot{Iname}=N;
                else
                    Ntot{Iname}=Ntot{Iname}+N;
                end
                cum=cumsum(N)/sum(N);
                [junk,F95]=min(abs(0.95-cum));
                [junk,F75]=min(abs(0.75-cum));
                
                bar(F,N,'k');
                xlim([0 maxval(1,Iname)]);
                
                set(gca,'fontweight','bold','fontsize',14);
                ylabel('Number events');xlabel(unitts{Iname});grid on
                title(sprintf('mean value of %s on %s, %6.2f to %6.2f covers 75-95%% range', ...
                    names{Iname},date_range(Idate,:),F(F75),F(F95)));
                
                subplot(length(names),2,2*Iname)
                [N,F]=hist(tmp.(names{Iname}).std,linspace(0,maxval(2,Iname),250));
                if Idate==1,
                    Ntot_std{Iname}=N;
                else
                    Ntot_std{Iname}=Ntot_std{Iname}+N;
                end
                
                cum=cumsum(N)/sum(N);
                [junk,F95]=min(abs(0.95-cum));
                [junk,F75]=min(abs(0.75-cum));
                
                bar(F,N,'k');
                
                set(gca,'fontweight','bold','fontsize',14);
                xlim([0 maxval(2,Iname)]);
                
                ylabel('Number events');xlabel(unitts{Iname});grid on
                title(sprintf('std of %s on %s, %6.2f to %6.2f covers 75-95%% range', ...
                    names{Iname},date_range(Idate,:),F(F75),F(F95)));
                
                
                
                
            end
            keyboard;
            orient landscape
            print('-djpeg',sprintf('CallStats_%s',date_range(Idate,:)));
            
            
        end
        
        figure
        for Iname=1:length(names),
            subplot(length(names),2,2*Iname-1)
            cum=cumsum(Ntot{Iname})/sum(Ntot{Iname});
            [junk,F95]=min(abs(0.95-cum));
            [junk,F75]=min(abs(0.75-cum));
            F=linspace(0,maxval(1,Iname),250);
            bar(F,Ntot{Iname},'k');ylabel('Number events');xlabel(unitts{Iname});grid on
            title(sprintf('mean value of %s over all dates, %6.2f to %6.2f covers 75-95%% range', ...
                names{Iname},F(F75),F(F95)));
            
            cum=cumsum(Ntot_std{Iname})/sum(Ntot_std{Iname});
            [junk,F95]=min(abs(0.95-cum));
            [junk,F75]=min(abs(0.75-cum));
            subplot(length(names),2,2*Iname);
            F=linspace(0,maxval(2,Iname),250);
            bar(F,Ntot_std{Iname},'k');ylabel('Number events');xlabel(unitts{Iname});grid on
            title(sprintf('std value of %s over all dates, %6.2f to %6.2f covers 75-95%% range', ...
                names{Iname},F(F75),F(F95)));
            
        end
        
        orient landscape
        print('-djpeg',sprintf('CallStats_total'));
        
        
    end
end

function struct_trim=trim_struct(struct_in,Igood)
names=fieldnames(struct_in);
for I=1:length(names),
    struct_trim.(names{I})=struct_in.(names{I})(Igood,:);
end
end

function [site_number,DASARstr,file_templates]=lookup_files(J,date_range,goodDASAR,manualdir)
%if ~isempty(findstr(manualdir,'Shell07'))
   %  site_number=goodDASAR{J}{1}(2);  %DASAR a of this site
    %DASARstr=[];
   %  for K=1:length(goodDASAR{J}),
   %      DASARstr=[DASARstr goodDASAR{J}{K}(end-1)];
   %  end
   %  file_templates.base=['*wc' date_range(J,(end-3):end) '*s' site_number ];
     
 %else
    site_number=goodDASAR{J}{1}(2);  %DASAR a of this site
    file_templates.base=['*wc' date_range(J,(end-3):end) '*S' site_number ];
    DASARstr=[];
    for K=1:length(goodDASAR{J}),
        DASARstr=[DASARstr goodDASAR{J}{K}(end-1)];
    end
%end

file_templates.tsv=[file_templates.base '*.tsv'];
file_templates.mat=[file_templates.base '*.mat'];



end

function d = segment_points_dist_2d ( p1, p2, p )

%% SEGMENT_POINTS_DIST_2D: distance ( line segment, point ) in 2D.
%
%  Discussion:
%
%    A line segment is the finite portion of a line that lies between
%    two points.
%
%    The nearest point will satisfy the condition
%
%      PN = (1-T) * P1 + T * P2.
%
%    T will always be between 0 and 1.
%
%  Modified:
%
%    03 May 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real P1(2), P2(2), the endpoints of the line segment.
%
%    Input, real P(2), the point whose nearest neighbor on the line
%    segment is to be determined.
%
%    Output, real D, the distance from the point to the line segment.
%
%  Modified 19 February 2009 from Burkardt's SEGMENT_POINT_DIST_2D
%    to handle multiple points in an m-by-2 matrix (each row is a point).
%    Still assumes a single line segment.         Chris Nations

p1 = p1(:)';
p2 = p2(:)';
[m,n] = size(p);
if n~=2,
    error('p must be m-by-2')
end
one = ones(m,1);
p1 = one*p1;
p2 = one*p2;

%  If the line segment is actually a point, then the answer is easy.

if p1==p2,
    t = zeros(m,2);
else
    bot = sum((p2 - p1).^2,2);
    t = sum((p - p1).*(p2 - p1),2) ./ bot;
    t = max(t,0);
    t = min(t,1);
    t = t*[1 1];
end

pn = p1 + t.*(p2 - p1);
d = sqrt(sum((pn - p).^2,2));
end


