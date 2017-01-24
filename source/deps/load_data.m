%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function load_data.m%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,t,Fs,tmin,tmax,head]	=	load_data(filetype,tdate_start,tlen,Ichan,handles)

%%% tmin,tmax, t are datenumbers
%   x rows are samples, columns are channels
%   head.geom.rd:  If multiple channels are present, this gives spacing
persistent  Fs_keep keyword space

mydir=handles.mydir;
myfile=handles.myfile;
teager=get(handles.checkbox_teager,'Value');
set(handles.edit_normal_rotation,'Vis','Off');  %Set off for the moment
set(handles.text_normal_rotation,'Vis','Off');  %Set off for the moment

Fs=[];
x=[];
t=[];
tmin=[];
tmax=[];
head=[];
head.multichannel=false;
head.linked=false;

filetype	=	upper(filetype);

switch filetype
    case 'PSD'
        [x,F,t,Tabs,params]=read_Java_PSD(fullfile(mydir,myfile),tdate_start,tlen);
        if ~isempty(t)&&length(t)>1
            Fs=1./(t(2)-t(1));
            t=t-t(1);
        else
            Fs=1;
        end
        head=params;
        
        tmin=head.tstart_file;
        tmax=head.tend_file;
        t=linspace(tmin,tmax,length(x));
    case 'MAT'
        
        
        %%Look at THAW data
        if ~isempty(strfind(myfile,'rcv'))
            data=load(fullfile(mydir,myfile));
            Fs=data.sample_rate;
            x=data.rcvx;
            
            day=str2num(myfile(5:7));
            hour=str2num(myfile(8:9));
            minn=str2num(myfile(10:11));
            secc=str2num(myfile(12:13));
            tmin=datenum(2012,0,day,hour,minn,secc);
            tmax=tmin+datenum(0,0,0,0,0,length(x)/Fs);
            head.Nchan=1;
            head.Fs=Fs;
            
            if tdate_start==-1
                tdate_start=tmin;
                
            end
            dt=datevec(tdate_start-tmin);
            
            dt=dt(:,6)+60*dt(:,5)+3600*dt(:,4)+24*3600*dt(:,3);
            Istart=1+round(dt*Fs);
            Iend=Istart+(round(tlen*Fs));
            x=x(Istart:Iend);
            
            %keyboard
            return
        end
        simulated=load(fullfile(mydir,myfile));
        Fs=simulated.fs;
        
        %%Uncomment to conduct Bering Sea call selection
        %%% size(x)=[ 5(zplot)           4 (rplot)          8       16033]
        name='Bering Sea Simulations';
        numlines=1;
        if ndims(simulated.x)>3
            
            prompt={sprintf('Depths: %s',mat2str(simulated.zplot)), ...
                sprintf('Ranges: %s',mat2str(simulated.rplot)), ...
                sprintf('Modes: %s', mat2str(1:size(simulated.x,3)))};
            defaultanswer={'22','22000',''};
            
        else
            prompt={sprintf('Depths: %s',mat2str(simulated.zplot)), ...
                sprintf('Ranges: %s',mat2str(simulated.rplot))};
            defaultanswer={'22','22000'};
        end
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        [~,Iz]=min(abs(eval(answer{1})-simulated.zplot));
        [~,Ir]=min(abs(eval(answer{2})-simulated.rplot));
        if isempty(deblank(answer{3}))
            Imode=[];
        else
            Imode=eval(answer{3});
        end
        if ndims(simulated.x)>3
            if isempty(Imode)
                simulated.x=squeeze(sum(simulated.x,3));
               
            else
                Imode=min([Imode size(simulated.x,3)]);
                simulated.x=squeeze(simulated.x(:,:,Imode,:));
            end
            %figure;spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)
        end
         x=squeeze(simulated.x(Iz,Ir,:));
        
       %%%%%%%%%%%%%
        %x=simulated.x_sweep';
        if strcmp(Ichan,'all')
            Ichan=1:size(x,2);
        end
        if max(Ichan)>max(size(x,2)),
            disp(['Channel too high, restricting to ' int2str(max(size(x,2)))]);
            Ichan=max(size(x,2));
        end
        
        x=x(:,Ichan);
        
        
        if tdate_start==-1
            tmin=datenum(0);
        elseif tdate_start>0 %We are trimming beginning time
            tmin=tdate_start;
            tmp=datevec(tdate_start);
            dn=1+floor(Fs*tmp(end));
            x=x(dn:end,:);
            t=t(dn:end);
        else
            tmin=tdate_start;
        end
        
        if tlen*Fs<size(x,1)
            x=x(1:floor(tlen*Fs),:);
        end
        t=(1:length(x))/Fs;
        
        tmax=tmin+datenum(0,0,0,0,0,max(t));
        head.geom.rd=simulated.rd;
        head.Nchan=size(x,2);
        x=x';
    case 'GSI'
        
        button_chc=get(handles.togglebutton_ChannelBeam,'String');
        beamform_data=0;
        get_geometry=0;
        if strcmpi(button_chc,'angle')&&~strcmpi(Ichan,'all')
            beamform_data=1;
            get_geometry=1;
        elseif strcmpi(Ichan,'all')
            Ichan=1:3;
            beamform_data=0;
            get_geometry=1;
        end
        
        if beamform_data==1
            thta=Ichan;
            Ichan=1:3;
        end
        
        if beamform_data==0
            if max(Ichan)>3
                disp('Channel too high, restricting to 3');
                Ichan=3;
            end
            
            if min(Ichan)<1
                disp('Channel too low, restricting to 1');
                Ichan=1;
            end
        end
        
        [x,t,head]=readGSIfile([mydir '/' myfile],tdate_start,tlen,Ichan,'datenum','nocalibrate');
        head.multichannel=true;
        head.linked=true;
        head.Nchan=3;
        
        if isempty(keyword)
            prompt = {'Enter a keyword for GSI calibration [DASARC]:'};
            dlg_title = 'DASAR calibration';
            num_lines = 1;
            def = {'DASARC'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            keyword=answer{1};
        end
        % keyword=input('Enter a keyword for GSI calibration [DASARC]:','s');
        %             if isempty(keyword)
        %                 keyword='DASARC';
        %             end
        %         end
        x=calibrate_GSI_signal(x, keyword);
        
        Fs=head.Fs;
        
        tmin=datenum(1970,1,1,0,0,head.ctbc);
        tmax=tmin+datenum(0,0,1,0,0,0);
        %head.Nchan=length(Ichan);
        if beamform_data==1
            %y=x(:,1)-fudge_factor_velocity*(cos(angles(I)*pi/180)*x(:,3)+sin(angles(I)*pi/180)*x(:,2)); %switch x and y to get compass bearing
            
            thta=thta-head.brefa;
            x=x(:,1)+sin(thta*pi/180)*x(:,2)+cos(thta*pi/180)*x(:,3);
            x=x/2; %Turn into equivalent of one channel.
            
        end
        
        
    case 'MT'
        %[x,t,Fs]=load_mt_mult(handles.mydir,tdate_start,tlen);
        %head=read_mt_header([mydir filesep myfile]);
        head=read_mt_header(fullfile(mydir, myfile));
        
        tmin=head.tstart;
        tmax=head.tend;
        Fs=head.Fs;
        
        if tdate_start>0
            tdate_vec=datevec(tdate_start-tmin);
            nsec=tdate_vec(6)+60*tdate_vec(5)+3600*tdate_vec(4);
        else
            %uiwait(msgbox('No start time requested'));
            head.multichannel=false;
            head.linked=false;
            
            return
        end
        %Load accelerometer data as desired
        %% Template: CB_Accel_X_2014_02100958.mt
        if strfind(head.abbrev,'Accel')
            mystr='XYZ';
            Ispace=strfind(myfile,'_');
            Ispace=Ispace(2)+1;
            x=[];
            for Iacc=1:3
                myfile_temp=myfile;
                myfile_temp(Ispace)=mystr(Iacc);
                try
                    [x0,t]=load_mt([mydir filesep myfile_temp],nsec,tlen);
                    if isempty(x)
                        x=zeros(length(x0),3);
                    end
                    x(:,Iacc)=x0;
                catch
                    fprintf('Channel %s of %s not available...\n',mystr(Iacc),myfile_temp);
                    
                end
                
            end
            head.Nchan=3;
            
            %%Select channel requested...or rotate channels...
            % Axis conventions assume 2011 definitons: X is long axis, Y
            % and Z are normal.
            
            set(handles.edit_normal_rotation,'Vis','On');
            set(handles.text_normal_rotation,'Vis','On');
            normal_angle=(pi/180)*str2num(get(handles.edit_normal_rotation,'String'));
            
            if normal_angle~=0
                xrot=x;
                xrot(:,2)=x(:,2)*cos(normal_angle)+x(:,3)*sin(normal_angle);
                xrot(:,3)=-x(:,2)*sin(normal_angle)+x(:,3)*cos(normal_angle);
                x=xrot;
                head.normal_rotation=(180/pi)*normal_angle;
                head.transformation.description{1}=' Normal Rotation: %6.2f deg';
                head.transformation.value{1}=head.normal_rotation;
            else
                head.transformation.description{1}=' No rotation';
                head.transformation.value{1}=[];
            end
            
            x=x(:,Ichan);
            %myfile(Ispace)=mystr(Ichan);
            
            head.myfile=myfile;
            
            % guidata(handles.myfile);  %Saves new file name...
            %guidata(hObject, handles);  %Saves new file name...
            %guidata(hObject, handles);
            %
        else
            if Ichan>1
                disp('WARNING: load_data: MT can only have one channel');
                Ichan=1;
            end
            if tdate_start>0
                [x,t]=load_mt(fullfile(mydir , myfile),nsec,tlen);
            end
            head.Nchan=1;
        end
    case 'SIO'
        %x,t,Fs,tmin,tmax,head
        %load_data(filetype,tstart_min,tdate_start,tlen,Ichan,handles)gff
        %set(handles.text_channel,'String','Channel [-angle (deg)]');
        sio_chc=get(handles.togglebutton_ChannelBeam,'String');
        [~, head] = sioread(fullfile(mydir,myfile));
        head.multichannel=true;
        
        %Extract start date and time
        %if ~exist('tstart_min') || isempty(tstart_min) || tstart_min < 0
        Idot=strfind(myfile,'.');
        
        if isfield(head,'date')
            tmin=head.date;
        else
            success=0;
            for II=1:(length(Idot)-1)
                if success || Idot(II+1)-Idot(II)~= 12
                    continue
                end
                template=myfile((Idot(II)+1):(Idot(II+1)-1));
                try
                    tmin=datenum(2000+str2num(template(1:2)),0,str2num(template(3:5)), ...
                        str2num(template(6:7)),str2num(template(8:9)),str2num(template(10:11)));
                    %tstart_min=tmin;
                    success=1;
                catch
                    tmin= now;
                end
            end
        end
        %end
        
        if isfield(head,'Fs')
            Fs=head.Fs;
        else
            Fs_default=25000;
            if isempty(Fs_keep)
                Fs=input('Enter a sampling frequency in Hz per channel [25000]:','s');
                if isempty(Fs),Fs=Fs_default;end
                head.Fs=Fs;
                Fs_keep=Fs;
            else
                Fs=Fs_keep;
            end
        end
        
        %%	Data parameters needed
        head.np		=	head.PperChan;
        tmax =tmin+datenum(0,0,0,0,0,head.np/Fs);
        head.Nchan	=	head.N_Chan;
        
        if ~exist('tdate_start') || isempty(tdate_start) || tdate_start < 0
            x=[];
            return
        end
        
        if tdate_start>0
            np_start=1+round((tdate_start-tmin)*24*3600*Fs);
        else
            
            np_start=1;
        end
        if np_start>head.np
            errordlg('load_data:  SIO start time exceeds points available!');
            return
        end
        
        
        %%%If beamforming is desired for a look at a given direction...
        beamform_data=0;
        get_geometry=0;
        if strcmpi(sio_chc,'angle')&&~strcmpi(Ichan,'all')
            beamform_data=1;
            get_geometry=1;
        elseif strcmpi(Ichan,'all')
            Ichan=1:head.Nchan;
            beamform_data=0;
            get_geometry=1;
        end
        
        if beamform_data==1
            thta=-Ichan;
            Ichan=1:head.Nchan;
        end
        
        %If loading single channel, check that request is reasonable...
        if beamform_data==0
            if max(Ichan)>head.Nchan
                errordlg(sprintf('load_data:  SIO channel request: %i, max channels: %i',max(Ichan),head.Nchan));
                return
            end
            
            if max(Ichan)<1
                errordlg(sprintf('load_data:  SIO channel request: %i is less than 1',max(Ichan)));
                return
            end
            
        end
        
        npi=round(tlen*Fs);
        if np_start+npi>head.np
            errordlg('load_data:  SIO end time exceeds points available!');
            return
        end
        
        
        
        if get_geometry==1
            if isempty(space)
                prompt = {'Enter spacing[m] between elements for SIO file:'};
                dlg_title = 'SIO file spacing';
                num_lines = 1;
                def = {'0.1'};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
                space=eval(answer{1});
                fprintf('Half-wavelength frequency: %6.2f Hz\n',1500/(2*space));
            end
            head.geom.rd=(0:(head.Nchan-1))*space;
            
            
        end
        
        
        [x,~]=sioread(fullfile(mydir,myfile),np_start,npi,Ichan);
        %Data arranged so that time are rows, columns are channels
        
        %Flip data ....
        
        if size(x,2)>1
            x=fliplr(x);
        end
        
        if beamform_data==1
            try
                space=head.geom.rd(2)-head.geom.rd(1);
                xtot=delaynsum(x,thta,space,Fs,Ichan);
                x=xtot;
                head.thta=thta;
            catch
                disp('load_data: sioread beamform failure');
                keyboard
            end
        end
        
        %Set time
        t=(1:size(x,1))/Fs;
        
        
    case 'DAT'
        
        %           [x,head]=read_synchronized_mdat_files(fullfile(mydir,myfile),tdate_start,tlen);
        %         Fs=head.fs;
        %         sio_chc=get(handles.togglebutton_ChannelBeam,'String');
        %         head.multichannel=true;
        %         t=(1:length(x))/Fs;
        %
        %         tmin=head.tfs;
        %         tmax=head.tfe;
        %         head.Nchan=size(x,2);
        %         beamform_data=0;
        %         get_geometry=1;
        
        dat_chc=get(handles.togglebutton_ChannelBeam,'String');
        [x,tmin,tmax,Fs,head]=read_synchronized_dat_file([mydir '/' myfile],tdate_start,tlen,0); %output in uPa
        %[x,tmin,tmax]=read_dat_file([mydir '/' myfile],Fs,tdate_start,tlen,1);  %Voltage output
        t=(1:size(x,1))/Fs;
       
       
        head.Nchan=size(x,2);
        
        if head.Nchan>1
             head.multichannel=true;
        end
         %%%If beamforming is desired for a look at a given direction...
        beamform_data=0;
        get_geometry=0;
        if strcmpi(dat_chc,'angle')&&~strcmpi(Ichan,'all')
            beamform_data=1;
            get_geometry=1;
        elseif strcmpi(Ichan,'all')
            Ichan=1:head.Nchan;
            beamform_data=0;
            get_geometry=1;
        end
        
        if beamform_data==1
            thta=-Ichan;
            Ichan=1:head.Nchan;
            try
                space=head.geom.rd(2)-head.geom.rd(1);
                xtot=delaynsum(x,thta,space,Fs,Ichan);
                x=xtot;
                head.thta=thta;
            catch
                disp('load_data: dat beamform failure');
                keyboard
            end
        else
            if ~strcmp(Ichan,'all')
                x		=	x(:,Ichan);
            end
        end
        
        
        
%         switch Fs
%             case 50000
%                 head.calcurv=[
%                     1.179288464673746e+06
%                     -3.417289147406752e+06
%                     3.972100408634462e+06
%                     -2.459193259685826e+06
%                     8.904700994689314e+05
%                     -1.924134277822444e+05
%                     2.476608484423531e+04
%                     -2.235739303825218e+03
%                     2.904887584919255e+02
%                     -5.381149759460806e+00
%                     7.841554559708414e-03
%                     ];
%             case 6250
%                 head.calcurv=[
%                     -9.864342626384007e+06
%                     2.675183405132254e+07
%                     -3.072255757018830e+07
%                     1.946983114345214e+07
%                     -7.445224085881455e+06
%                     1.766054734429601e+06
%                     -2.570588847834060e+05
%                     2.188411119767746e+04
%                     -9.803725367146685e+02
%                     1.959124505642275e+01
%                     -2.811936435415921e-01];
%             case 12500
%                 head.calcurv=[
%                     9.262441626302190e+07
%                     -2.151487191990283e+08
%                     2.069375942078056e+08
%                     -1.063702102525421e+08
%                     3.153159716612202e+07
%                     -5.458152352141772e+06
%                     5.382152297627985e+05
%                     -2.765563363629215e+04
%                     6.113088605208859e+02
%                     -1.301582987521525e+00
%                     -1.634557871607174e-01];
%         end
%         
    case 'ADI'
        [x,tmin,tmax,fs]=read_adi_file(mydir,myfile,[],0,tlen,0);
        
        if isempty(fs)
            fs=input('Enter sampling rate in Hz:');
        end
        Fs=fs;
        [x,tmin,tmax]=read_adi_file(mydir,myfile,Fs,tdate_start,tlen,0); %Output in uPa
        %[x,tmin,tmax]=read_adi_file(mydir,myfile,Fs,tdate_start,tlen,1);  %Output in voltage
        t=(1:length(x))/fs;
        head.Nchan=1;
        
    case 'MDAT'
        
        
        [x,head]=read_synchronized_mdat_files(fullfile(mydir,myfile),tdate_start,tlen);
        Fs=head.fs;
        sio_chc=get(handles.togglebutton_ChannelBeam,'String');
        head.multichannel=true;
        t=(1:length(x))/Fs;
        
        tmin=head.tfs;
        tmax=head.tfe;
        head.Nchan=size(x,2);
        beamform_data=0;
        get_geometry=1;
        
        if strcmpi(sio_chc,'angle')&&~strcmpi(Ichan,'all')
            beamform_data=1;
            get_geometry=1;
        elseif strcmp(Ichan,'all')
            Ichan=1:head.Nchan;
            beamform_data=0;
            get_geometry=1;
        end
        
        if beamform_data==1
            thta=-Ichan;
            Ichan=1:head.Nchan;
        end
        
        if max(Ichan)>head.Nchan
            disp(sprintf('Channel too high, max channel %i',head.Nchan));
            Ichan=head.Nchan;
        end
        
        try
            x=x(:,Ichan);
            for II=Ichan
                fprintf('Logged depth of channel %i is %6.2f m\n',II,head.geom.rd(II));
            end
        end
        
        x=x';
        
        if beamform_data==1
            try
                space=head.geom.rd(2)-head.geom.rd(1);
                xtot=delaynsum(x',thta,space,Fs,Ichan);
                x=xtot;
                head.thta=thta;
            catch
                disp('load_data: MDAT beamform failure');
                keyboard
            end
        end
        
        %Set time
        t=(1:length(x))/Fs;
        head.calcurv=[
            -3.826740371096994e+07
            8.955964717194067e+07
            -9.041409824956748e+07
            5.195683344725087e+07
            -1.873830725900567e+07
            4.336813690670102e+06
            -6.251710636164501e+05
            5.339062841084976e+04
            -2.452316303554637e+03
            5.178423592184026e+01
            -1.229223450332553e+00];
        
    case {'M4V','MP4'}
        info	=	audioinfo(fullfile(mydir,myfile));
        info_vid = mmfileinfo(fullfile(mydir,myfile));
        %[~,Fs]		=	audioread(fullfile(mydir,myfile),1,'native');
        Nsamples	=	info.TotalSamples;
        handles.Fs	=	info.SampleRate;
        Fs=info.SampleRate;
        
        %Attempt to get date
        tmin=get_start_time_MP4(mydir,myfile,info);
        %tmin=datenum([1970 1 1 0 0 0]);
        tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
        %tlen=Nsamples/handles.Fs;
        
        %a 2 V pk-pk signal yields 0.9 pk-pak amp at 3 kHz or channel 2
        %(audio)
        sens=.9*2/0.9;  %Convert to volts
        sens_hydro=172; %dB re uPa/V
        sens=sens*10.^(sens_hydro/20);
        %sens=1;
        
        
        
        tdate_vec	=	datevec(tdate_start - tmin);
        nsec		=	tdate_vec(6) + 60*tdate_vec(5) + 3600*tdate_vec(4);  %Ignores differences in days
        N1			=	1 + round(nsec*handles.Fs);
        N2			=	N1 + round(tlen*handles.Fs);
        
        try
            [x,Fs]		=	audioread(fullfile(mydir,myfile),[N1 N2],'native');
            %Quick adjust for 6 dB rolloff of ultrasonic (channel 1)
            
        catch
            x=[];
            
            t=[];
            head.Nchan=0;
            return
        end
        
        head.Nchan	=	size(x,2);
        if head.Nchan>1
            head.multichannel=true;
        end
        
        %%Option for "adaptive processing" to remove interference...
        % Set channel equal to -1 to give difference between signal
        
        if ~strcmp(Ichan,'all')&Ichan>0
            try
                x(:,1)=[0; diff(x(:,1))];
            end
            x		=	x(:,Ichan);
            
        elseif Ichan<0  %for audio channel will be correct response
            x=x(:,2)-x(:,1);
        end
        
        t	=	(1:length(x))/Fs;
        
        x			=	double(x)*sens;
    case 'WAV' %% Includes SUDAR data
        
        %Check version of MATLAB
        
        if verLessThan('matlab', '8.2.0.29')
            [~,Fs]=wavread(fullfile(mydir,myfile),3);
            Nsamples=max(wavread(fullfile(mydir,myfile),'size') );
            handles.Fs=Fs;
            
            % Put code to run under MATLAB older than MATLAB 7.0.1 here
        else
            % Put code to run under MATLAB 7.0.1 and newer here
            info	=	audioinfo(fullfile(mydir,myfile));
            %[~,Fs]		=	audioread(fullfile(mydir,myfile),1,'native');
            Nsamples	=	info.TotalSamples;
            handles.Fs	=	info.SampleRate;
            Fs=info.SampleRate;
        end
        
        [head.cable_factor,sens]=get_ADAT24_cable_factor;
        
        try
            done=false;
            [SUDAR_true,tmin,tmax,FsSUDAR]=get_SUDAR_time(mydir,myfile); %Check whether a SUDAR file exists
            
            if SUDAR_true
                sens=(10^(186/20))/(2^15);
                Fs=FsSUDAR;
                done=true;
                
            end
            
            if ~done
                try
                    [Berchok_true,tmin,tmax]=get_Berchok_time(mydir,myfile,Nsamples,Fs); %Check whether a Catherine Berchok file exists
                    if Berchok_true
                        sens=1;
                        
                        % kludge Code for testing beta=-3;
%                         disp('WARNING!  DOWNLOADING TEST FILE');
%                         data=load('/Users/thode/Projects/RightWhaleDetection/Publications/Refractive_Invariant/DataExample.dir/Animal_faint/AveragedCallBeta1Range15kmModes1n2/StackedResult.mat');
%                       
%                         x=data.xtot;
%                         %x=data.xref;
%                         Fs=data.Fs;
%                         t=1:length(x);
%                         tmin=datenum([1970 1 1 0 0 0]);
%                         tmax=tmin+datenum(0,0,0,0,0,length(x)/Fs);
%                         head.Nchan=1;
%                         
%                         return
                        
                        %%End of kludge
                        done=true;
                    end
                catch
                    %keyboard
                end
            end
            
            
            if ~isempty(strfind(myfile,'LL017_Set4_3min.wav.x.wav'))
                sens=3162;
            end
            
            
            if ~done
                tmin	=	convert_date(myfile,'_');
                if isempty(tmin)
                    tmin=datenum([1970 1 1 0 0 0]);
                end
                tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
            end
        catch
            disp([myfile ': convert_date failure']);
            try
                tmin=datenum(get(handles.text_mintime,'String'));
            catch
                minn	=	input('Enter start date in format [yr mo day hr min sec] or hit return: ');
                if isempty(minn)
                    minn=[1970 1 1 0 0 0];
                end
                tmin	=	datenum(minn);
            end
            tmax	=	tmin + datenum(0,0,0,0,0,Nsamples/Fs);
        end
        
        
        tdate_vec	=	datevec(tdate_start - tmin);
        nsec		=	tdate_vec(6) + 60*tdate_vec(5) + 3600*tdate_vec(4);  %Ignores differences in days
        N1			=	1 + round(nsec*handles.Fs);
        N2			=	N1 + round(tlen*handles.Fs);
        
        try
            if verLessThan('matlab', '8.2.0.29')
                [x,Fs]		=	wavread(fullfile(mydir,myfile),[N1 N2],'native');
            else
                [x,Fs]		=	audioread(fullfile(mydir,myfile),[N1 N2],'native');
            end
        catch
            x=[];
            
            t=[];
            head.Nchan=0;
            return
        end
        
        head.Nchan	=	size(x,2);
        if head.Nchan>1
            head.multichannel=true;
        end
        if ~strcmp(Ichan,'all')
            x		=	x(:,Ichan);
        end
        
        t	=	(1:length(x))/Fs;
        
        x			=	double(x)*sens;
        
end

if isempty(tmin) || isempty(tmax)
    disp('load_data: Warning, tmin and tmax should never be empty when exiting..');
end


%%Store whether multichannel data, regardless of original file format.
if min(size(x))>1
    head.multichannel=true;
end

if ~isfield(head,'linked')
    head.linked=false;
end

if ~isfield(head,'multichannel')
    head.multichannel=false;
end

%%%Optional Teager-Kaiser filtering...
if teager
    %%Assume that x is in form [ channel time]
    x=x(:,2:end-1).^2-x(:,1:end-2).*x(:,3:end);
end


end  %function load_data

