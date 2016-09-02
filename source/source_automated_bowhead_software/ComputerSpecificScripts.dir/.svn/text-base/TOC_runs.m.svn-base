%%%TOC_runs.m%
%[Site_vector,date_str,DASAR_str,keyword,detection_stage,datatype]=TOC_runs(Icase),
%
%% Select scope of time and space for processing DASAR sites
%% Note that these values can be overridden by load_local_params.m, but generally dates for manual analyses and
%%          DASARs used are stored here
%
% function [Site_vector,date_str,DASAR_str]=TOC_runs(Icase),
%  Input: Icase-keyword eg. 'Site_4_short.morph.crosschecked'
%         Format is 'scenario.algorithm.process_stage', divided by periods.
%                algorithm:  %'contour','morph','both','none':
%               process_stage:  String of concatenated phrases like
%                   'RawCrosscheckedRandomClassified'.
%  Output:
%       Site_vector:  vector of indicies between 1 and 5, representing a site
%       date_str:  Two element cell array with strings of year/month/day; e.g. '20090902'
%           is 2009 Sept. 2.  A set of multiple entries, e.g.
%               date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};
%           defines a set of non-contiguous dates
%       DASAR_str='a' through 'g', or '*' for all DASARS at a site...
%       keyword:  structure with three fields, 'scenario','algorithm', and
%               'stage', to be used for outputing exact automated output
%               results
%       detection_stage:  cell array of detection stage strings, for use in
%           master_evaluation.m..

function [Site_vector,date_str,DASAR_str,keyword,detection_stage]=TOC_runs(Icase)

Site_vector=[];
date_str=[];
DASAR_str=[];
detection_stage=[];

%datatype='short';
%%Divide keyword into three parts.
keyword.scenario=[];
keyword.algorithm=[];
keyword.stage=[];

Idots=findstr(Icase,'.');
if isempty(Idots),
    keyword.scenario=Icase;
elseif length(Idots)==1,
    keyword.scenario=Icase(1:(Idots(1)-1));
    keyword.algorithm=Icase((Idots(1)+1):end);
elseif length(Idots)==2,
    keyword.scenario=Icase(1:(Idots(1)-1));
    keyword.algorithm=Icase((Idots(1)+1):(Idots(2)-1));
    keyword.stage=Icase((Idots(2)+1):end);
elseif length(Idots)==3,
    keyword.scenario=Icase(1:(Idots(1)-1));
    keyword.algorithm=Icase((Idots(1)+1):(Idots(2)-1));
    keyword.stage=Icase((Idots(2)+1):(Idots(3)-1));
end

%%Create a 'case string' from the Icase input for use in switch statement
Nsites=5;
years=7:30;  %We'll worry about 2031 when we come to it.
for I=years
    if ~isempty(findstr(Icase,sprintf('Shell%02i',I)))
        for J=0:Nsites
            if ~isempty(findstr(Icase,sprintf('Site%i',J)))
                case_word=sprintf('Shell%02i_Site%i',I,J);
                
                if ~isempty(findstr(Icase,'airgun'))
                    case_word=[case_word '_airgun'];
                end
            end
        end
    end
    
    if ~isempty(findstr(Icase,sprintf('AURAL%02i',I)))
        for J=1:Nsites
            if ~isempty(findstr(Icase,sprintf('Site%i',J)))
                case_word=sprintf('AURAL%02i_Site%i',I,J);
                
                if ~isempty(findstr(Icase,'airgun'))
                    case_word=[case_word '_airgun'];
                end
            end
        end
    end
    
    if ~isempty(findstr(Icase,sprintf('BP%02i',I)))
        for J=1:Nsites
            if ~isempty(findstr(Icase,sprintf('Site%i',J)))
                case_word=sprintf('BP%02i_Site%i',I,J);
                
                if ~isempty(findstr(Icase,'airgun'))&&J==1
                    case_word=[case_word '_airgun'];
                end
            end
        end
    end
    
end

disp(case_word)

switch case_word
    
    %%%%%2011%%%%%
    case {'Shell11_Site1','Shell11_Site1_airgun'}
        Site_vector=[1];
        date_str{1}={'20110725'};
        date_str{2}={'20111010'};
        %DASAR_str='*';
        DASAR_str='abcdef';  %No g
    case {'Shell11_Site2','Shell11_Site2_airgun'}
        Site_vector=[2];
        date_str{1}={'20110725'};
        date_str{2}={'20111010'};
        DASAR_str='*';
    case {'Shell11_Site3','Shell11_Site3_airgun'}
        Site_vector=[3];
        date_str{1}={'20110725'};
        date_str{2}={'20111010'};
        DASAR_str='acdef';  %No bg
        %DASAR_str='af';
    case {'Shell11_Site4','Shell11_Site4_airgun'}
        Site_vector=[4];
        date_str{1}={'20110725'};
        date_str{2}={'20111010'};
        DASAR_str='*';
    case {'Shell11_Site5','Shell11_Site5_airgun'}
        Site_vector=[5];
        date_str{1}={'20110725'};
        date_str{2}={'20111010'};
        DASAR_str='*';
        
        %%%%%%%2010%%%%%%%%
        
    case {'Shell10_Site1','Shell10_Site1_airgun'}
        Site_vector=[1];
        date_str{1}={'20100801'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20101006'};  %End date
        DASAR_str='*';
        DASAR_str='*';
        
        
    case {'Shell10_Site2','Shell10_Site2_airgun'}
        Site_vector=[2];
        date_str{1}={'20100801'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20101006'};  %End date
        DASAR_str='*';
    case {'Shell10_Site3','Shell10_Site3_airgun'}
        Site_vector=[3];
        date_str{1}={'20100801'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20101006'};  %End date
        DASAR_str='bcdeg';
    case {'Shell10_Site4','Shell10_Site4_airgun'}
        Site_vector=[4];
        date_str{1}={'20100801'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20101006'};  %End date
        DASAR_str='*';
    case {'Shell10_Site5','Shell10_Site5_airgun'}
        Site_vector=[5];
        date_str{1}={'20100801'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20101006'};  %End date
        DASAR_str='*';
    case {'Shell07_Site1','Shell07_Site1_airgun'},
        Site_vector=[ 1];
        %date_str{1}={'20070820'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{1}={'20070910'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20071015'};  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'Shell08_Site1',
        Site_vector=[ 1];
        date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %End date
        
        date_str{1}={'20080812'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{1}={'20080918'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN *re-run*
        date_str{2}={'20081008'};  %End date
        
        %date_str{1}={'20080916'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN *re-run*
        %date_str{2}={'20080916'};  %End date
        
        
        %DASAR_str='aceghijk';
        %DASAR_str='d';
        DASAR_str='acdeghijk';  %%I added 1F back into mix June 30, 2010
    case {'Shell09_Site1','Shell09_Site1_airgun'}
        Site_vector=[ 1];
        date_str{1}={'20090823'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20091005'};  %End date
        %date_str{1}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{1}={'20090907'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090907'};  %End date
        
        %DASAR_str='acdeghijk';
        DASAR_str='abcdef';
    case 'Shell08_Site1_airgun',
        Site_vector=[ 1];
        
        date_str{1}={'20080821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081005'};  %End date
        
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'BP08_Site1_airgun',
        Site_vector=[ 1];
        
        date_str{1}={'20080826'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080924'};  %End date
        
        DASAR_str='j'; %Use '*' to process all DASARs
    case 'BP09_Site1_airgun', %Kat use this
        Site_vector=[ 1];
        
        %%%date_str{1}={'20090826'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{1}={'20090925'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20090929'};  %End date
        
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'BP10_Site1_airgun', %Kat use this
        Site_vector=[ 1];
        
        date_str{1}={'20100808'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20101015'};  %End date
        
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'BP12_Site1_airgun', %Kat use this
        Site_vector=[ 1];
        
        date_str{1}={'20120826'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20121007'};  %End date
        
        DASAR_str='*'; %Use '*' to process all DASARs       
    case {'Shell07_Site2','Shell07_Site2_airgun'}
        Site_vector=[ 2];
        %date_str{1}={'20070820'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{1}={'20070910'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20071015'};  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'Shell08_Site2',
        Site_vector=[ 2];
        date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %End date
        
        %date_str{1}={'20080816'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20081007'};  %End date
        %date_str{2}={'20080901'};  %End date
        
        DASAR_str='*';
    case 'Shell08_Site2_airgun',
        Site_vector=[ 2];
        date_str{1}={'20080824'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081005'};  %End date
        DASAR_str='f'; %Use '*' to process all DASARs
        
    case {'Shell09_Site2','Shell09_Site2_airgun'}
        Site_vector=[ 2];
        
        
        date_str{1}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        
        date_str{1}={'20090820'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20091006'};  %End date
        
        %date_str{1}={'20090912'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090912'};  %End date
        
        DASAR_str='*';
        
    case {'Shell07_Site3','Shell07_Site3_airgun'}
        Site_vector=[ 3];
        date_str{1}={'20070820'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20071015'};  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
        
    case 'Shell08_Site3'
        Site_vector=[ 3];
        date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %End date
        
        date_str{1}={'20080821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081006'};  %End date
        
        DASAR_str='abcdfg'; %Use '*' to process all DASARs
    case 'Shell08_Site3_airgun',
        Site_vector=[ 3];
        date_str{1}={'20080824'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081005'};  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
        
    case {'Shell09_Site3','Shell09_Site3_airgun'}
        Site_vector=[ 3];
        
        date_str{1}={'20090821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20091001'};  %End date
        
        
        %date_str{1}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        DASAR_str='abcdfg'; %Use '*' to process all DASARs
        DASAR_str='*';
        DASAR_str='abcdef';
        
    case {'Shell07_Site4','Shell07_Site4_airgun'}
        Site_vector=[ 4];
        date_str{1}={'20070820'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20071015'};  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'Shell08_Site4'
        Site_vector=[ 4];
        date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %End date
        
        date_str{1}={'20080820'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081005'};  %End date
        
        DASAR_str='bcefg'; %Use '*' to process all DASARs
    case 'Shell08_Site4_airgun',
        Site_vector=[ 4];
        date_str{1}={'20080821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081005'};  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
        
    case {'Shell09_Site4','Shell09_Site4_airgun'}
        Site_vector=[ 4];
        
        
        date_str{1}={'20090817'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20091004'};  %End date
        
        DASAR_str='*'; %Use '*' to process all DASARs
        %DASAR_str='e'; %Use '*' to process all DASARs
        
        %date_str{1}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        
    case 'Shell08_Site5_airgun'
        Site_vector=[ 5];
        
        date_str{1}={'20080819'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20080913'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        date_str{2}={'20081010'};  %End date
        
        %date_str{1}={'20080831'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20080831'};  %End date
        
        DASAR_str='g'; %Use '*' to process all DASARs
    case {'Shell07_Site5','Shell07_Site5_airgun'}
        Site_vector=[ 5];
        date_str{1}={'20070821','20070828','20070906','20070913','20070921','20070929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20070821','20070828','20070906','20070913','20070921','20070929'};  %End date
        
        date_str{1}={'20070821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20071010'};  %End date
        
        DASAR_str='*';
    case 'Shell08_Site5',
        Site_vector=[5];
        date_str{1}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080821','20080828','20080906','20080913','20080921','20080929'};  %End date
        
        date_str{1}={'20080819'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20081003'};  %End date
        
        date_str{1}={'20080831'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20080831'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'Shell09_Site5',
        Site_vector=[ 5];
        
        date_str{1}={'20090819'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20091005'};  %End date
        
        
        %date_str{1}={'20090912'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090912'};  %End date
        
        
        %date_str{1}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090827','20090903','20090908''20090912','20090914','20090918','20090925','20090930'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        DASAR_str='*'; %Use '*' to process all DASARs
        %       DASAR_str='e';
        
    case {'Shell09_Site5_airgun'}
        Site_vector=[ 5];
        
        date_str{1}={'20090821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20091004'};  %End date
        
        %date_str{1}={'20090831'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090831'};  %End date
        
        DASAR_str='*'; %Use '*' to process all DASARs
        
    case 'Shell09_Final_February',
        Site_vector=[ 5];
        date_str{1}={'20090821','20090828','20090906','20090913','20090921','20090929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}={'20090821','20090828','20090906','20090913','20090921','20090929'};  %End date
        
        %date_str{1}={'20090828','20090906','20090913','20090921','20090929'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090828','20090906','20090913','20090921','20090929'};  %End date
        
        %date_str{1}={'20090821'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090821'};  %End date
        
        %date_str{1}={'20090905'};  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}={'20090905'};  %End date
        
        %datatype='short';  %Since DASAR-A are 53-B sonobuoys, have a frequency-dependent calibration.
        DASAR_str='*'; %Use '*' to process all DASARs
        
        
        
    case 'NorthStar08',
        Site_vector=0;
        date_str{1}='20090827';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20090925';  %End date
        %datatype='short';  %Since DASAR-A are 53-B sonobuoys, have a frequency-dependent calibration.
        DASAR_str='*'; %Use '*' to process all DASARs
        
        %Site_vector=;
        %date_str{1}='20090907';  %Two distant airgun surveys
        %date_str{2}='20090907';  %End date
        
        %date_str{1}='20090920';  %Strong airguns
        %date_str{2}='20090920';  %End date
        
        %datatype='double';  %Since DASAR-A are 53-B sonobuoys, have a frequency-dependent calibration.
        %DASAR_str='j'; %Use '*' to process all DASARs
        
    case 'Liberty08_Aug3',
        Site_vector=5;
        date_str{1}='20090803';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20090803';  %End date
        %datatype='short';  %Since DASAR-A are 53-B sonobuoys, have a frequency-dependent calibration.
        DASAR_str='a'; %Use '*' to process all DASARs
        
        %Site_vector=;
        %date_str{1}='20090820';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{2}='20090820';  %End date
        datatype='short';  %Since DASAR-A are 53-B sonobuoys, have a frequency-dependent calibration.
        %DASAR_str='a'; %Use '*' to process all DASARs
    case 'Liberty08_all',
        
        datatype='short';  %Since DASAR-A are 53-B sonobuoys, have a frequency-dependent calibration.
        Site_vector=[5 6 8];
        date_str{1}='20090801';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20090826';  %End date
        DASAR_str='a'; %Use '*' to process all DASARs
    case 'Site 1 all',
        Site_vector=1;
        date_str{1}='20070824';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20071012';  %End date
        DASAR_str='ac'; %Use '*' to process all DASARs
        
    case 'Site 2 all',
        Site_vector=2;
        date_str{1}='20070823';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20071011';  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'all sites',
        Site_vector=2:5;
        date_str{1}='20070823';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{1}='20070918';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RU
        date_str{1}='20070918';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20070918';  %End date
        DASAR_str='beg'; %Use '*' to process all DASARs
    case 'Site_4_short',
        Site_vector=4;
        %         date_str{1}='20070913';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %         date_str{2}='20070916';  %End date
        %         DASAR_str='af'; %Use '*' to process all DASARs
        %
        date_str{1}='20070829';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20070829';
        
        %date_str{2}='20070830';  %End date
        DASAR_str='abcdfg'; %Use '*' to process all DASARs
    case 'Site 4 airgun',
        Site_vector=4;
        %         date_str{1}='20070913';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %         date_str{2}='20070916';  %End date
        %         DASAR_str='af'; %Use '*' to process all DASARs
        %
        date_str{1}='20070915';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20070920';  %End date
        DASAR_str='*'; %Use '*' to process all DASARs
    case 'Site 3 short',
        Site_vector=3;
        date_str{1}='20070823';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='200701008';  %End date
        
        date_str{1}='20070918';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        date_str{2}='20070918';  %End date
        
        DASAR_str='h'; %Use '*' to process all DASARs
    case 'Site 4 long',
        Site_vector=4;
        date_str{1}='20070823';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        %date_str{1}='20070920';  %Start date-NOTE HAS TO BE EXACT! OTHERWISE PROGRAM WILL NOT RUN
        
        date_str{2}='20071010';  %End date
        DASAR_str='abcdfg'; %Use '*' to process all DASARs
        
end




%Pull out algorithm type:  %'contour','morph','both','none':


%%Pull out detection stages.
Icount=1;
detection_stage=[];
if findstr(lower(keyword.stage),'raw')>0,
    detection_stage{Icount}='raw';Icount=Icount+1;
end
if findstr(lower(keyword.stage),'classified')>0,
    detection_stage{Icount}='classified';Icount=Icount+1;
end
if findstr(lower(keyword.stage),'crosschecked')>0,
    detection_stage{Icount}='crosschecked';Icount=Icount+1;
end
if findstr(lower(keyword.stage),'random')>0,
    detection_stage{Icount}='random';Icount=Icount+1;
end


%[detection_stage,dates_list,DASAR_list,autoCase,DASAR_list_total]
%%Master reference list of DASAR data...
