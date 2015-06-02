function [modes,Nmode,choice,filt] = modeSelect(x_ok,Fe,N,Nww,r,c,clims,filt)

% Warping
[s_w, Fe_w]=warp_temp_exa(x_ok,Fe,r,c);    % s_w: warped signal, Fe_w: new warping frequency
M=length(s_w);
t_w=(0:M-1)/Fe_w;                           % Warped time
f_w=(0:M-1)*Fe_w/M;                         % Warped frequencies
rtf=tfrstft(s_w,1:M,M,hamming(Nww));
rtf(floor(M/2):end,:)=[];                   % Delete negative frequencies
RTF=abs(rtf);

workspace;
choice=str2double(input(['\nAre you satisfied with the source deconvolution ? \n' ...
'1: Back to source deconvolution\n'...
'2: Back to time selection (for warping)\n'...
'3: Continue\n'],'s'));

if choice==3

    disp('Beginning warping')

    % Mode number
    Nmode=nan;
    while isnan(Nmode) % if the user make a misclic in the mode selection, press esc to do it again
        fig=figure;
        imagescFun(t_w,f_w,RTF,'ij')
        ylim(filt.fwlims)
        caxis(clims)
        title('Enter a mode number')
        Nmode=str2double(input('Number of modes ? ','s'));
        close(fig)
    end
    
    filt.('mask')=zeros(Nmode,length(rtf(:,1)),M);
    modes=zeros(Nmode,N);

    % Filtering (Mask)
    for mm=1:Nmode
        BW=[];
        while isempty(BW) % if the user make a misclic in the mode selection, press esc to do it again
            disp(strcat('Filtering mode ',int2str(mm)))
            fig=figure('units','normalized','outerposition',[0 0 1 1]);
            imagescFun(t_w,f_w,RTF,'ij')
            caxis(clims)
            ylim(filt.fwlims)
            title(['Draw a polygon to select the mode ' num2str(mm) ' (Press esc to do it again)'])
            BW = roipoly ;
            close(fig)
        end
        masque=zeros(size(rtf));
        masque(1:length(RTF(:,1)),:)=BW;
        mode_rtf_warp=masque.*rtf;
        % We don't forget the window normalization factor to calculate the time...
        % signal from the tf response, as well as the *2 factor because the
        % selection is done only on the positive frequencies:
        [mode_temp_warp]=real(sum(mode_rtf_warp,1))*2/M/max(hamming(Nww)/norm(hamming(Nww)));

        modes(mm,:)=iwarp_temp_exa(mode_temp_warp,Fe_w,r,c,Fe,N);
        disp(strcat('Mode ',int2str(mm),' filtered'))
        filt.('mask')(mm,:,:)=masque;       
        clear masque mode_rtf_warp mode_temp_warp
    end

    else
        Nmode=0;
        modes=[];
end
disp('End of warping')
end