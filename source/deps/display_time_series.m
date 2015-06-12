function display_time_series(handles,x,Fs,hdr)
    
    %%Check that we are looking at acoustic data and are in spectrogram
    %%mode
    if isempty(strfind(handles.myfile,'Press'))
        % msgbox('Enter min and max frequency, or hit return to skip filtering:','modal');
        
        if strcmp(handles.old_display_view,'Spectrogram')
            ButtonName = questdlg('Filter? (If yes, click on two frequency values in spectrogram)');
        else
            ButtonName='No';
        end
        
        %if ~isempty(tmp)&&size(tmp,1)==2
        %%Check that are on current spectrogram view
        if strcmp(ButtonName,'Yes')
            tmp=ginput(2);
            
            freq=sort(tmp(:,2))*1000;
            minfreq=freq(1);maxfreq=freq(2);
            %y=quick_filter(x(:,1),Fs,freq(1),freq(2))
            frange=[0.8*minfreq minfreq maxfreq maxfreq+0.2*minfreq];
            [N,Fo,Ao,W] = firpmord(frange,[0 1 0],[0.05 0.01 0.1],Fs);
            B = firpm(N,Fo,Ao,W);
            
            y=filter(B,1,x(:,1)-mean(x(:,1)));
        else
            %y=x(:,1)-mean(x(:,1));
            y=x(:,1);
        end
    else
        y=x(:,1)/1000; %Pressure conversion
    end
    
    t=(1:length(x(:,1)))/Fs;
    xlabel('Time (sec)');
    if isfield(hdr,'calunits')&&~isempty(strfind(hdr.calunits,'mPa'))
        plot(handles.axes1,t,1000*y);grid on;
        ylabel('uPa');
    else
        plot(handles.axes1,t,y);grid on;
        try
            ylabel(hdr.calunits);
        catch
        ylabel('Amplitude');
    end
    end
    end