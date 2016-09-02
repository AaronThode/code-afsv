%%%%%review_snips_fast.m%%%%%
%  Demonstration of how to use readEnergySummary and readEnergySnips
%  Checks to see that 'persistent' assignment to fid allows for faster file
%  reading..
%  Start from within the Processed/Site_0* folder...
%   Variables you can change within the script:
%   ndays: number of days to display at a time in an interval plot
%   Idir_start:  Index of starting directory to process. To process entire
%       site, set the value to '1'.
%   Ifile_start: Index of starting file within a directory.  Usually set to
%       '1' unless you are trying to review a specific file.
%   yes_interactive: If set to one, program will ask you to click on the
%       figure and return a time and inteval value.  Useful as input into
%       GSI_specgram_viewer.m
%   Plots are stored in FinalResults directory.

clear all;close all
%interval_remove=[19.75 22];
ndays=5; %Number of days to display at a time...
Idir_start=1;
Ifile_start=1;
yes_interactive=0;

Fs=1000;
% mydir=pwd;
mydir='/Volumes/macmussel1/Arctic_2007/Processed/Site_04'
cd /Volumes/macmussel1/Arctic_2007/Processed/Site_04
!mkdir FinalResults
dirname=dir('D07*');
for Idir=Idir_start:length(dirname),
    cd(mydir);
    cd(dirname(Idir).name);
    fnames=dir('*.snips');
    cum_time=[];
    cum_interval=[];
    cum_npt=[];
    for Ifile=Ifile_start:length(fnames),

        %Ifile=input('Enter index for snips file:');

        Idot=max(findstr(fnames(Ifile).name,'.'))-1;
        fname=fnames(Ifile).name(1:Idot);

        indexx=Inf;
        %[data,head,Igood]=readEnergySummary([fname '.detsum'], indexx,interval_remove);
        [data,head,Igood]=readEnergySummary([fname '.detsum'], indexx);
        cum_time=[cum_time datenum(1970,1,1,0,0,data.ctime)];
        cum_interval=[cum_interval data.interval/head.Fs];
        cum_npt=[cum_npt data.npt(1,:)/head.Fs];
        figure(1);
        subplot(2,1,1);
        hh=plot(cum_time,cum_interval,'g.');ylabel('samples between snips in sec');
        set(hh,'markersize',2);
        datetick('x',6);ylim([0 30]);grid on;
        title(sprintf('Results of %i files',Ifile));
        subplot(2,1,2);
        hh=plot(cum_time,cum_npt,'g.');ylabel('length of snip in sec');
        set(hh,'markersize',2);
        datetick('x',6);grid on;
        keyboard;
        if rem(Ifile,ndays)==0,
            if yes_interactive==1,
                disp('Select a time');
                tmp=ginput(1);
                disp(sprintf('At time %s interval is %6.2f',datestr(tmp(1)),tmp(2)));
                pause;
            end

            orient tall
            subplot(2,1,1);
            title(sprintf('Detection intervals between %s and %s',datestr(min(cum_time),1),datestr(max(cum_time),1)));
            pause(1);
            currentdir=pwd;
            cd(mydir);
            cd('FinalResults');
            Idot=max(findstr(fnames(Ifile).name,'_tstart'))-1;
            fname1=fnames(Ifile).name(1:Idot);
            print('-djpeg',[fname1 '_intervals.jpg']);
            cd(currentdir);
            cum_interval=[];
            cum_npt=[];
            cum_time=[];


        end
    end

end