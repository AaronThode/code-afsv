%%%%%%%%Energy_detector_review.m%%%%%%%%%%%%%%%
%function Energy_detector_review(param,fname,tstart,twant,tlen,climm,hfigs, tview)
%
%  Aaron Thode
%  December 9. 2008
%  Run a short segment of time through Energy detector for insight and
%  adjustment.
%
%  Input:
%       param: structure of parameters for energy detector; must contain
%           field param.energy.
%       fname: complete pathname to raw data file
%       tstart: datenumber of start time of file
%       twant:  datenumber of desired start time energy detector processing
%       tlen:   number of seconds to process
%       climm:  intensity bounds for images
%       hfig:  handles to axes to plot final detections on
%       tview: two element vector giving bounds for viewing output (sec)
%            This permits viewing of a short transient window while
%            permitting a long "start up" time for the equalization.

%%Load parameters
%path(path,'../CommonScripts.dir');
%path(path,'../EvaluateLongRuns.dir');

function data_all=Energy_detector_review(param,fname,tstart,twant,tlen,climm,hfigs,tview)
%function data_all=Energy_detector_review(param,fname,tfs,twant,tsec,tview)
!rm *.sel *detsum *snips *.txt

%params_chc='Shell08_initial';
%param=TOC_params(params_chc);  %Unpack detection parameters from keyword
%param.energy.threshold=5;

%%Define target file
%fname='/Volumes/macmussel1/Arctic_2008/Site_05/S508A0/S508A0T20080821T000000.gsi';
%head=readgsif_header(fname);

%Define desired time
%twant=datenum(2008,8,21,0,0,0);
%twant=datenum(2008,8,21,0,50,0);

mydir=pwd;
if ~exist('climm')
    climm=[80 100];
end
if ~exist('hfigs')
    hfigs=[];
end
    
if ~exist('tview')
    xlimm=[0 tlen];
else
    xlimm=tview;
end

%%%%%%%%%%%%%%%%%%%%%%%
param.energy.Fs=param.Fs;
if twant>0&tstart>0
    tmp=datevec(twant-tstart);
    tmp=tmp(:,6)+60*tmp(:,5)+3600*tmp(:,4)+24*3600*tmp(:,3);
else
    tmp=0;
end
param.energy.nstart=floor(param.Fs*tmp);
param.energy.nsamples=floor(tlen*param.Fs);


%param.energy=alter_parameters(param.energy);
%%Write energy detector scripts
if ~isempty(findstr(lower(computer),'mac')),
    slashstr='/';
else
    slashstr='\';
end

Islash=max(findstr(fname,slashstr));
param.energy.exten=fname((Islash+1):end);
param.energy.file_dir=fname(1:(Islash-1));

hh=figure;
hfigs=[hfigs hh];

for J=1:3
    param.energy.debug=J;  %Equalized SNR
    
    scriptname=writeEnergyDetectorCshell(param.energy);
    disp('Running JAVA script');
    %!rm *snips *detsum
    eval(sprintf('! ./%s > outtt.txt',scriptname));
    cd(param.energy.dir_out)
    selnames=dir('*.sel');
    for I=1:length(selnames)
        [SEL,Tsec{I},Tabs{I},data]=read_Java_SEL(selnames(I).name,0,Inf);
        if I==1,
            SEL_all=zeros(length(selnames),length(SEL));
            
        end
        try
            SEL_all(I,1:length(SEL))=SEL';
        catch
            disp(selnames(I).name);
            disp(sprintf('size of matrix: %i, size of incoming: %i',size(SEL_all,2),length(SEL)));
        end
        fmin(I)=data.fmin;
        fmax(I)=data.fmax;
        
    end
    cd(mydir)
    
    [fmin,Isort]=sort(fmin);
    fmax=fmax(Isort);
    SEL_all=SEL_all(Isort,:);
    
    subplot(3,1,J);
    
    imagesc(Tsec{1},0.5*(fmin+fmax),10*log10(SEL_all));axis('xy')
    set(gca,'fontweight','bold','fontsize',14);
    xlabel('Time (seconds)');ylabel('Frequency (Hz)');
    xlim(xlimm)
    ylim([0 param.energy.f_high])
    switch(J)
        case 1
            title(sprintf('%s, Raw SEL (dB re 1 uPa^2-s) measured over %6.2f Hz segments, starting at %s', ...
                fname, param.energy.bandwidth,datestr(Tabs{1}(1))));
            caxis(climm);colorbar('EastOutside');
        case 2
            title(sprintf('Equalized background estimate (dB re 1 uPa^2-s), equalization timescale=%6.2f s',param.energy.eq_time));
            caxis(climm);colorbar('EastOutside');
        case 3
            title(sprintf('Detection function with lower threshold %6.2f dB, max %i dB',param.energy.threshold,30));
            caxis([param.energy.threshold 30]);
            colorbar('EastOutside');
    end
    %colorbar
end

cd(param.energy.dir_out)

detsum=dir([ '*detsum']);
[data_all,head]=readEnergySummary(detsum.name, Inf);

%cd(mydir);
if isempty(data_all.ctime)
    errordlg('No events detected!');
    return
end
tmin=datestr(datenum(1970,1,1,0,0,data_all.ctime(1)));
tmax=datestr(datenum(1970,1,1,0,0,data_all.ctime(end)));
fprintf('Start time: %s End time:%s\n',datestr(tmin),datestr(tmax));

%%Feature names
%     'min_freq'
%     'min_duration'
%     'max_freq'
%     'max_duration'
%     'peak_freq'
%     'peak_duration'
%     'peak'
%     'ctime_peak'
%     'total_duration'

for J=1:length(hfigs)
    if length(hfigs)==2&&J==1
        axes(hfigs(1));
        fscale=1000;
    else
        figure(hfigs(J));
        fscale=1;
    end
    for I=1:length(data_all.nstart)
        fmin=data_all.features(1,I)/fscale;
        tmin=data_all.features(2,I)/param.Fs; %duration of minimum frequency bin
        fmax=data_all.features(3,I)/fscale;
        tmax=data_all.features(4,I)/param.Fs;  %duration of maximum frequency bin
        fpeak=data_all.features(5,I)/fscale;
        tpeak=data_all.features(6,I)/param.Fs;  %duration of peak frequency bin
        ctime_peak=data_all.features(8,I);
        %Start time of peak detection in seconds.
            %For AFSV window, subtract the burn in time (represented by
            %tview(1).
        ttstart=(data_all.nstart(I)/param.Fs)-(2-J)*tview(1);  
        tduration=data_all.npt(I)/param.Fs;  %duration across all bands
        
        %Vertical line of frequency extent
        
        
        hold on; hh=line([1 1]*ttstart,[fmin fmax]); set(hh,'Color','w','linewidth',3);
        
        %hh=line([0 data_all.npt(I)/param.Fs]+data_all.nstart(I)/param.Fs,0.5*(fmin+fmax)*[1 1]); set(hh,'Color',[1 1 1],'linewidth',3);
        
        hh=line([0 tmin]+ttstart,[fmin-5/fscale fmin-5/fscale]); set(hh,'Color','g','linewidth',3);  %Min frequency
        hh=line([0 tmax]+ttstart,[fmax+5/fscale fmax+5/fscale]); set(hh,'Color','g','linewidth',3);   %Max frequency band
        %hh=line([0 tpeak]+tstart,[fpeak fpeak]); set(hh,'Color','r','linewidth',3);  %Peak frequency
        hh=line([0 tduration]+ttstart,[fpeak fpeak],'color','k','linewidth',5,'linestyle','-');
        plot(ttstart+(ctime_peak-data_all.ctime(I)),fpeak,'rs','linewidth',3);
    end
    hold off;
end

disp('Hit return to see alternate view of bottom subplot:')

pause;
caxis([0 20]);
hold off


end

%%%%%%%alter_parameters.m%%%%
function [param,change_flag]=alter_parameters(param)

names=fieldnames(param);
names{length(names)+1}='Exit';
Ichc=menu('Select a parameter:',names);
change_flag=0;
double_level=0;
if ~strcmp(names{Ichc},'Exit')&&isstruct(param.(names{Ichc}))
    double_level=1;
    names2=fieldnames(param.(names{Ichc}));
    Ichc2=menu('Select a parameter:',names2);
else
    double_level=0;
    
end


while Ichc~=length(names)
    change_flag=1;
    
    if double_level==0,
        current_val=param.(names{Ichc});
    else
        current_val=param.(names{Ichc}).(names2{Ichc2});
    end
    
    for I=1:length(current_val),
        
        if double_level==0,
            val=input(sprintf('Change %s = %6.2f into what? (Return to leave unchanged)',names{Ichc},current_val(I)));
            if ~isempty(val)
                param.(names{Ichc})(I)=val;
            end
            
        else
            val=input(sprintf('Change %s.%s = %6.2f into what? (Return to leave unchanged)',names{Ichc},names2{Ichc2},current_val(I)));
            if ~isempty(val)
                param.(names{Ichc}).(names2{Ichc2})(I)=val;
            end
            
            
        end
    end
    Ichc=menu('Select a parameter:',names);
    
    double_level=0;
    if ~strcmp(names{Ichc},'Exit')&&isstruct(param.(names{Ichc}))
        double_level=1;
        names2=fieldnames(param.(names{Ichc}))
        Ichc2=menu('Select a parameter:',names2);
    else
        double_level=0;
        
    end
end

end

