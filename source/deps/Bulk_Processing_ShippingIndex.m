%%%%Bulk_Processing_ShippingIndex.m%%
%% Aaron Thode
%% May 21, 2013
%
%%  Load raw acoustic data and compute shipping metrics using shipdet_metrics.m
% Input parameters:
%
%   interval: time interval to skip for generating metric samples
%   param.tsample:  data window in seconds over which to process data for a single metric sample
%   param.threshold=5; dB value to feed into shipdet_metrics
%   param.df_search=[10 20].  Bandwidth over which to seach for the threshold drop...
%   param.f_min,f_max:  frequency range to search for tones/bands...
%   param.ovlap,Nfft:  spectrogram generation parameters.
%   chc: {'med','avg'}: med: take median of spectrogram window rows; 'avg', take average

% Output Data structure:
%   metric{IIsites}{IIDASAR}
%   The Sites index has value of 1-6 (6=Site '0')
%   The DASAR index has 'A'=1, 'B'=2, etc.
%      has following fields:
%       tabs:  datenumbers for all dates deployed..
%       entropy (from median spectrum)
%       kurtosis (from median spectrum)
%       .isi., adp., isi or adaptive peak algorithm
%           .avg,.med: whether median or averaging used to obtain spectrum
%               .peaks{Npeaks}: different bandwidths uses on peak detectors.
%                   Fpeak, pwr: band frequencies (Hz) detected; total power of industrial activity in dB re
%                   1uPa
%
%  The last DASAR index contains all metrics from all DASARS combined into one series (Site-wide metrics).
%  Not used in 2012 analysis--removed too many calls

clear
rawdatadir='/Volumes/Shell2012_GSI_Data';
Sites=[ 3 4];
DASAR_str='*';
Icase='Shell12';
date_str{1}={'20120729'};
date_str{2}={'20121129'};

%date_str{1}={'20120929'};
%date_str{2}={'20120929'};

%%Ship metric parameters...
interval=0.5;  %minutes
tabs_interval=datenum(0,0,0,0,interval,0);
param.tsample=30; %seconds
param.threshold=5; %dB
param.df_search=[10 20]; %+-/Hz to examine
param.f_min	=	25;
param.f_max	=	450;
param.ovlap=0.75;
param.Nfft=512;  %FFT size for spectrogram
param.Idebug=0;  %If 0, no plotting of output...
chc={'med','avg'};

%%%%%%%%%%%%%%%%No adjustments below this line%%%%%%%%%%%%%%%%%
Npeaks=length(param.df_search);  %Number of paramater sets to run on peak picker...
datelist=expand_date_range(date_str);

for Isites=1:length(Sites)
    if Sites(Isites)==0
        IIsites=6;
    else
        IIsites=Sites(Isites);
    end
    
    for Idate=1:size(datelist,1)
        
        [goodDASAR,goodFile,goodNames,DASAR_str_out]=find_DASAR_dates(datelist(Idate,:),Sites(Isites),DASAR_str,rawdatadir,'Shell12');
        for Ifile=1:length(goodFile)
            fprintf('Starting %s\n',goodFile{Ifile});
            
            IIDASAR=double(DASAR_str_out(Ifile))-64; %'A' is put in index 1, etc...
            
            %Create a time vector
            head=readgsif_header(goodFile{Ifile});
            tabs=head.tabs_start:tabs_interval:head.tabs_end;
            
            %Initialize metric cells if they do not exist yet..
            initialize=0;
            try
                metric{IIsites}{IIDASAR};
            catch
                initialize=1;
            end
            
            if initialize==1||isempty(metric{IIsites}{IIDASAR})
                fprintf('Initializing Site %i, DASAR %s\n',IIsites,DASAR_str_out(Ifile));
                maxsize=ceil(90*24*60/interval); %90 day deployment, 24 hours a day, 60 minutes an hour, divided by interval in minutes
                metric{IIsites}{IIDASAR}.Icount=1;
                metric{IIsites}{IIDASAR}.tabs=zeros(1,maxsize);
                metric{IIsites}{IIDASAR}.entropy=zeros(1,maxsize);
                metric{IIsites}{IIDASAR}.kurtosis=zeros(1,maxsize);
                for III=1:Npeaks
                    metric{IIsites}{IIDASAR}.isi{III}.pwr=zeros(1,maxsize);
                    metric{IIsites}{IIDASAR}.isi{III}.Fpeaks=zeros(10,maxsize,'single');
                    metric{IIsites}{IIDASAR}.adp{III}.pwr=zeros(1,maxsize);
                    metric{IIsites}{IIDASAR}.adp{III}.Fpeaks=zeros(10,maxsize,'single');
                end
            end
            
            %%Initialize metric structure
            %	Jit: read whole file before processing
            Tmin	=	head.ctbc;
            Tmax	=	head.ctec;
            Tlen	=	Tmax - Tmin;
            it_win	=	floor(head.Fs*param.tsample);
            it_inc	=	floor(head.Fs*interval*60);
            %	Read the whole file
            fprintf('Loading %s\n',goodFile{Ifile});
            tic
            [X, Tx, ~] =	readgsi(goodFile{Ifile},0,Tlen);
            disp('Done.');
            toc
            X	=	X(1,:);		%	keepy only sound data
            Ntx		=	length(Tx);
            tic
            disp('Calibrating');
            X=calibrate_GSI_signal(X,'DASARC');
            disp('Done.');
            toc
            
            for I=1:(length(tabs)-1)
                it_start	=	1 + it_inc*(I-1);
                
                if rem(I,30)==0,fprintf('Loading %i s at %s, %6.2f percent done\n',param.tsample,datestr(tabs(I)),100*I/(length(tabs)-1)); end
                try
                    %[x,t,head]=readgsi(goodFile{Ifile},tabs(I),param.tsample,'datenum');
                    %[S,FF,TT,B] = spectrogram(x(1,:),hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,head.Fs);
                    
                    %Fetch current segment...
                    it_end		=	it_start + it_win - 1;
                    if it_end > Ntx
                        it_end	=	Ntx;
                    end
                    x		=	X(it_start:it_end);
                    [S,FF,TT,B] = spectrogram(x,hanning(param.Nfft),round(param.ovlap*param.Nfft),param.Nfft,head.Fs);
                    
                catch
                    disp('Crash! Recover and continue...');
                    metric{IIsites}{IIDASAR}.Icount=metric{IIsites}{IIDASAR}.Icount+1;
                    continue
                end
                data.PSD	=	B;
                data.T		=	TT;
                data.F		=	FF;
                %data.Fs		=	head.Fs;
                %data.Nfft	=	param.Nfft;
                
                %Compute metrics and add to growing structure.
                %                 if I==16
                %                     figure(1);
                %                     imagesc(TT,FF/1000,10*log10(B))
                %                     axis('xy');
                %
                %                     [single_metric,hf]	=	shipdet_metrics(data, param,1,[]);
                %                 else
                [single_metric,hf]	=	shipdet_metrics(data, param);
                
                %   end
                Icount=metric{IIsites}{IIDASAR}.Icount;
                metric{IIsites}{IIDASAR}.tabs(Icount)=tabs(I);
                metric{IIsites}{IIDASAR}.entropy(Icount)=single_metric.entropy.med;  %Median of one-minute interval;
                metric{IIsites}{IIDASAR}.kurtosis(Icount)=single_metric.kurtosis.med;  %Median of one-minute interval;
                
                for Ic=1:2
                    for III=1:Npeaks
                        
                        value1=single_metric.(chc{Ic}).peaks{III}.isi.F;
                        if ~isempty(value1)
                            metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).Fpeaks(1:length(value1),Icount)=value1';
                            metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).pwr(Icount)=single_metric.(chc{Ic}).peaks{III}.isi.PdBall;
                        end
                        
                        value2=single_metric.(chc{Ic}).peaks{III}.adp.F;
                        if ~isempty(value2)
                            metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).pwr(Icount)=single_metric.(chc{Ic}).peaks{III}.adp.pdBinttotal;
                            metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).Fpeaks(1:length(value2),Icount)=value2';
                        end
                    end
                end
                
                metric{IIsites}{IIDASAR}.Icount=metric{IIsites}{IIDASAR}.Icount+1;  %Advance counter.
                
            end %I (time sample)
        end %Ifile (DASAR...
    end %Idate
    disp('begin cleanup and union of frequencies');
    
    
    NDASAR=length(metric{IIsites});
    first_DASAR=0;
    for Ic=1:2
        for III=1:Npeaks
            
            for IIDASAR=1:NDASAR
                
                if ~isempty(metric{IIsites}{IIDASAR})
                    fprintf('DASAR: %i\n',IIDASAR);
                    if first_DASAR==0
                        first_DASAR=1;
                    else
                        first_DASAR=-1;
                    end
                    Igood=find(metric{IIsites}{IIDASAR}.tabs>0);  %Identify unused reserved memory
                    
                    metric{IIsites}{IIDASAR}.tabs=metric{IIsites}{IIDASAR}.tabs(Igood);
                    metric{IIsites}{IIDASAR}.entropy=metric{IIsites}{IIDASAR}.entropy(Igood);
                    metric{IIsites}{IIDASAR}.kurtosis=metric{IIsites}{IIDASAR}.kurtosis(Igood);
                    metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).pwr=metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).pwr(Igood);
                    metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).Fpeaks=metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).Fpeaks(:,Igood);
                    metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).pwr=metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).pwr(Igood);
                    metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).Fpeaks=metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).Fpeaks(:,Igood);
                    
                    %Merge frequecy detections across all DASARs
                    if first_DASAR==1
                        disp('First non-empty DASAR slot');
                        tabs_all=metric{IIsites}{IIDASAR}.tabs;
                        Fpeaks_isi_all=metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).Fpeaks;
                        Nf_isi_all=size(Fpeaks_isi_all,1);
                        
                        Fpeaks_adp_all=metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).Fpeaks;
                        Nf_adp_all=size(Fpeaks_adp_all,1);
                        
                        
                        %=
                        %metric{IIsites}{NDASAR+1}.isi{III}.(chc{Ic}).Fpeaks=;
                        %metric{IIsites}{NDASAR+1}.adp{III}.(chc{Ic}).Fpeaks=
                    else
                        
                        tabs=metric{IIsites}{IIDASAR}.tabs;
                        Fpeaks_isi=metric{IIsites}{IIDASAR}.isi{III}.(chc{Ic}).Fpeaks;
                        Nf_isi=size(Fpeaks_isi,1);
                        
                        Fpeaks_adp=metric{IIsites}{IIDASAR}.adp{III}.(chc{Ic}).Fpeaks;
                        Nf_adp=size(Fpeaks_adp,1);
                        
                        [tabs_all,Isort]=sort([tabs_all tabs]);
                        
                        if Nf_isi>Nf_isi_all %If more peaks available
                            Fpeaks_isi_all=[Fpeaks_isi_all; zeros(Nf_isi-Nf_isi_all,size(Fpeaks_isi_all,2))];
                            Nf_isi_all=Nf_isi;
                            fprintf('line 59: Nf_isi_all: %i\n',Nf_isi_all);
                        elseif Nf_isi<Nf_isi_all
                            Fpeaks_isi=[Fpeaks_isi; zeros(Nf_isi_all-Nf_isi,size(Fpeaks_isi,2))];
                        end
                        Fpeaks_isi_all=[Fpeaks_isi_all Fpeaks_isi];
                        Fpeaks_isi_all=Fpeaks_isi_all(:,Isort);
                        
                        if Nf_adp>Nf_adp_all %If more peaks available
                            Fpeaks_adp_all=[Fpeaks_adp_all; zeros(Nf_adp-Nf_adp_all,size(Fpeaks_adp_all,2))];
                            Nf_adp_all=Nf_adp;
                            fprintf('line 69: Nf_adp_all: %i\n',Nf_adp_all);
                        elseif Nf_adp<Nf_adp_all
                            Fpeaks_adp=[Fpeaks_adp; zeros(Nf_adp_all-Nf_adp,size(Fpeaks_adp,2))];
                        end
                        
                        
                        Fpeaks_adp_all=[Fpeaks_adp_all Fpeaks_adp];
                        Fpeaks_adp_all=Fpeaks_adp_all(:,Isort);
                        
                        for Itabs=1:(length(tabs_all)-1)
                            if rem(Itabs,20000)==0,disp(Itabs/length(tabs_all));end
                            if tabs_all(Itabs+1)==tabs_all(Itabs)
                                tmp=union(Fpeaks_isi_all(:,Itabs), Fpeaks_isi_all(:,Itabs+1));
                                tmp=tmp(2:end); %Remove 0
                                if length(tmp)>Nf_isi_all
                                    Fpeaks_isi_all=[Fpeaks_isi_all; zeros(length(tmp)-Nf_isi_all,size(Fpeaks_isi_all,2))];
                                    Nf_isi_all=length(tmp);
                                    fprintf('line 85: Nf_isi_all: %i\n',Nf_isi_all);
                                end
                                Fpeaks_isi_all(1:length(tmp),Itabs+1)=tmp;
                                Fpeaks_isi_all((1+length(tmp)):end,Itabs+1)=0;
                                
                                tmp=union(Fpeaks_adp_all(:,Itabs), Fpeaks_adp_all(:,Itabs+1));
                                tmp=tmp(2:end);  %Remove 0
                                if length(tmp)>Nf_adp_all
                                    Fpeaks_adp_all=[Fpeaks_adp_all; zeros(length(tmp)-Nf_adp_all,size(Fpeaks_adp_all,2))];
                                    Nf_adp_all=length(tmp);
                                    fprintf('line 95: Nf_adp_all: %i\n',Nf_adp_all);
                                end
                                Fpeaks_adp_all(1:length(tmp),Itabs+1)=tmp;
                                Fpeaks_adp_all((1+length(tmp)):end,Itabs+1)=0;
                                
                                tabs_all(Itabs)=-1;  %Marker for later removal
                            end
                            
                        end %Itabs
                        
                        Ipass=find(tabs_all>0);
                        tabs_all=tabs_all(Ipass);
                        Fpeaks_adp_all=Fpeaks_adp_all(:,Ipass);
                        Fpeaks_isi_all=Fpeaks_isi_all(:,Ipass);
                        fprintf('Final number of times: %i\n',length(Ipass));
                        
                        
                    end %first_DASAR
                end %isempty
            end %IIDASAR
            metric{IIsites}{NDASAR+1}.tabs=tabs_all;
            metric{IIsites}{NDASAR+1}.isi{III}.(chc{Ic}).Fpeaks=Fpeaks_isi_all;
            metric{IIsites}{NDASAR+1}.adp{III}.(chc{Ic}).Fpeaks=Fpeaks_adp_all;
        end %III
    end  %Ic
    disp('cleanup finished');
    save(sprintf('Shell2012_shipping_metrics_combined_Site%i',IIsites),'-v7.3','metric','param');
    
    
end  %Isite

