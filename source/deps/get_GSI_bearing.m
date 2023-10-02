function [thet, kappa, tsec, xf] = get_GSI_bearing(app, ~, ~, handles, tmp)
% Freq or temp in kiloHz
%tsec: seconds into time series that data are selected...
% tmp is previous ginput
thet=-1;  %Start with failed result
kappa=-1;
tsec=-1;
tdate_start=handles.tdate_start;
tlen=handles.tlen;
%yes_wav=get(handles.togglebutton_getwav,'value');

mydir=pwd;

%Ichan='all';  %Hardwire first channel
%Ichan=str2double(get(handles.edit_chan,'String'));
azigram_flag=handles.radiobutton_azimuth.Value;
if azigram_flag
    disp('Can''t get bearing from azigram')
    return
end
DIFAR_flag=strcmpi(handles.filetype,'wav')&~isempty(strfind(handles.myfile,'DIFAR'));
GSI_flag=strcmpi(handles.filetype,'gsi')&~azigram_flag;

if GSI_flag %GSI file displayed in spectrogram view, use traditional bearing method
    [x,t,Fs,tstart,tend,head]=load_data(handles.filetype,tdate_start,tlen,'all',handles,app);
elseif DIFAR_flag
    %channel=str2num(handles.edit_chan.String);
    disp('Cannot use this button for DIFAR files:  use Select Points instead');
    set(handles.text_filename,'String','Cannot use this button for DIFAR files:  use Select Points instead');

    return
end
disp('Click on two extreme corners, click below axis twice to reject:');
if ~exist('tmp','var')
    tmp=ginput(2);
end
if (isempty(tmp)||any(tmp(:,2)<0))
    return
end
tsec=min(tmp(:,1));
n=round(Fs*sort(tmp(:,1)));

freq=1000*sort(tmp(:,2));
contents=get(handles.popupmenu_Nfft,'String');
Nfft=str2double(contents{get(handles.popupmenu_Nfft,'Value')});

if GSI_flag
    params.f_transition=300;
    x=x(n(1):n(2),:);
    [thet0,kappa,sd,xf]=extract_bearings(app, x,0.05,Nfft,Fs,freq(1),freq(2),200,params.f_transition);

    if ~isempty(strfind('T2007',handles.myfile))
        cal07flag=1;
        handles.calibration_DASAR2007_dir='/Users/thode/Projects/Greeneridge_bowhead_detection/Macmussel_Mirror/RawData';

        brefa_table=calibrate_bearing_Shell2007(handles.calibration_DASAR2007_dir,handles.myfile);
    else
        cal07flag=0;
    end

    if cal07flag==0
        thet=bnorm(thet0+head.brefa);
        if isnan(thet)
            keyboard
        end
    else
        [~,Icol]=calibrate_bearing_Shell2007(handles.calibration_DASAR2007_dir,handles.myfile,1);
        thet= interp1(0:360,brefa_table(:,Icol),bnorm(thet0));
    end
end

set(handles.text_filename,'String',sprintf('%s/%s %6.2f degrees... ',handles.mydir,handles.myfile,thet));
%keyboard;
end
