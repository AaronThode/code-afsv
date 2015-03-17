function [modes,Nmode,choice] = modeSelect(x_ok,Fe,N,r,c,soundsamp)

% Warping
[s_w, Fe_w]=warp_temp(x_ok,Fe,r,c);    % s_w: warped signal, Fe_w: new warping frequency
clear x_ok
M=length(s_w);
t_w=(0:M-1)/Fe_w;                           % Warped time
f_w=(0:M-1)*Fe_w/M;                         % Warped frequencies
Nww=1+2*ceil(ceil(M/4)/2);                               % Number of points for the tfr of the warped signal
rtf=tfrstft(s_w,1:M,M,hamming(Nww));

RTF=abs(rtf);
RTF(floor(M/2):end,:)=[];                   % Delete negative frequencies
RTF=RTF/max(max(RTF));

fig=imagescFun(t_w,f_w,RTF,[],'Warped signal');
axis ij

choice=str2double(input(['\nAre you satisfied with the source deconvolution ? \n' ...
'1: Back to signal selection\n'...
'2: Back to source deconvolution\n'...
'3: Back to time selection (for warping)\n'...
'4: Save and continue\n'...
'5: Continue\n'...
'6: Exit\n'],'s'));

if choice==4
    save([soundsamp '_warped_modes'])
end

if choice>=4 && choice <6

    disp('Beginning warping')
    disp('Choose the lowest and highest frequencies to zoom on the signal')
    title('Frequency zoom (select 2 points)')
    [~,y_lim]=ginput(2);
    [fw_min,fw_max]=deal(max(1,min(y_lim)),min(M/2,max(y_lim)));
    close(fig)

    % Mode number
    Nmode=nan;
    while isnan(Nmode) % if the user make a misclic in the mode selection, press esc to do it again
        fig=imagescFun(t_w,f_w,RTF,[fw_min,fw_max],'Enter a mode number');
        axis ij
        Nmode=str2double(input('Number of modes ? ','s'));
        close(fig)
    end

    modes=zeros(Nmode,N);

    % Filtering (Mask)
    for mm=1:Nmode
        BW=[];
        while isempty(BW) % if the user make a misclic in the mode selection, press esc to do it again
            disp(strcat('Filtering mode ',int2str(mm)))
            fig=imagescFun(t_w,f_w,RTF,[fw_min,fw_max],['Draw a polygon to select the mode ' num2str(mm) ' (Press esc to do it again)']);
            axis ij
            BW = roipoly ;
            close(fig)
        end
        masque=zeros(size(rtf));
        masque(1:length(RTF(:,1)),:)=BW;
        mode_rtf_warp=masque.*rtf;
        % We don't forget the window normalization factor to calculate the time...
        % signal from the tf response, as well as the *2 factor because the
        % selection is done only on the positive frequencies:
        [mode_temp_warp]=real(sum(mode_rtf_warp,1))/M/max(hamming(Nww)/norm(hamming(Nww)))*2;

        modes(mm,:)=iwarp_temp(mode_temp_warp,Fe_w,r,c,Fe,N);
        disp(strcat('Mode ',int2str(mm),' filtered'))
        clear masque mode_rtf_warp mode_temp_warp
    end

    else
        Nmode=0;
        modes=[];
end
disp('End of warping')
end