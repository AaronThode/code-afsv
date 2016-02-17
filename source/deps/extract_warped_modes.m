function [modes,Nmode,choice,filt] = extract_warped_modes(x_ok,Fe,N,Nww,r,c,clims,filt)

% Warping
[s_w, Fe_w]=warp_temp_exa(x_ok,Fe,r,c);    % s_w: warped signal, Fe_w: new warping frequency
M=length(s_w);
t_w=(0:M-1)/Fe_w;                           % Warped time
f_w=(0:M-1)*Fe_w/M;                         % Warped frequencies
rtf=tfrstft(s_w,1:M,M,hamming(Nww));
rtf(floor(M/2):end,:)=[];                   % Delete negative frequencies
M1=M/2;
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

    modes=zeros(Nmode,N);
    filt.('mask')=zeros(Nmode,length(rtf(:,1)),M);
    
    figure,
    imagescFun(t_w,f_w,RTF,'ij');
    ylim(filt.fwlims)
    title('Select the mode peaks (if you can''t see it press a key)')
    hold on,
    
%% Filtering (Mask)

    % The user selects approximately the peak of mode m
    [peaks_x,peaks_y,button]=ginput(Nmode);
    peaks=round([peaks_x*Fe_w peaks_y*M1/Fe_w]);
    
    % Choose directions around the peak to draw a polygon at max-thresh dB
    angle1=[0 0.005 0.01 0.03];                   % First quadran
    angle2=[angle1 0.25 flip(0.5-angle1(2:end))]; % 2 Firsts quadrans
    vect=exp(1i*2*pi*[angle2 angle2(1:end)+0.5]);
    nb_dir=length(vect);
    thresh=5; % [dB]

    points=zeros(Nmode,2,nb_dir);
   
    for mm=1:Nmode
        % Search for the real peak in a box around the selected points
        i_lims=200; 
        j_lims=15;
        maxi=max(max(RTF(max(1,peaks(mm,2)-j_lims):min(length(RTF(:,1)),peaks(mm,2)+j_lims),max(1,peaks(mm,1)-i_lims):min(length(RTF(1,:)),peaks(mm,1)+i_lims))));
        [mi,mj]=ind2sub(size(RTF),find(RTF==maxi));
        plot(mj/Fe_w,mi*Fe_w/M1,'kx')
        
        % Search in each direction the point where the threshold is reached
        for k=1:nb_dir
            v=[imag(vect(k)) real(vect(k))];
            bool=1;
            dv=1;
            while bool && mi+round(dv*v(1))>0 && mi+round(dv*v(1))<M && mj+round(dv*v(2)) && mj+round(dv*v(2))<M
                if 10*log10(RTF(mi+round(dv*v(1)),mj+round(dv*v(2))))<=10*log10(maxi)-thresh
                    bool=0;
                end
                dv=dv+1;
            end
            [points(mm,1,k),points(mm,2,k)]=deal(mi+round(dv*v(1)),mj+round(dv*v(2)));
        end
    end

    h=hggroup;
    for mm=1:Nmode
        if button(mm)==1 || button(mm)==2 || button(mm)==3 % if the user clicked on the peak
            % Create a draggable polygon
            for k=1:nb_dir-1
                if k*mm==1
                    tmp=imline(h,[points(mm,2,k)/Fe_w points(mm,1,k)/M1*Fe_w;points(mm,2,k+1)/Fe_w points(mm,1,k+1)/M1*Fe_w]);
                    line=repmat(tmp,nb_dir,Nmode);
                else
                    line(k,mm)=imline(h,[points(mm,2,k)/Fe_w points(mm,1,k)/M1*Fe_w;points(mm,2,k+1)/Fe_w points(mm,1,k+1)/M1*Fe_w]);
                end
            end
                line(nb_dir,mm)=imline(h,[points(mm,2,nb_dir)/Fe_w points(mm,1,nb_dir)/M1*Fe_w;points(mm,2,1)/Fe_w points(mm,1,1)/M1*Fe_w]);

            % Binds 2 consecutives lines
            for k=1:nb_dir-1
                line(k,mm).addNewPositionCallback(@(pos) link1(line,pos,k+1,mm));
                line(k+1,mm).addNewPositionCallback(@(pos) link2(line,pos,k,mm));
            end
            line(nb_dir,mm).addNewPositionCallback(@(pos) link1(line,pos,1,mm));
        end
    end
        % Let the user adjust his selection before going on
    msgbox('Press any key to continue')
    pause()
        
    for mm=1:Nmode
        if button(mm)==1 || button(mm)==2 || button(mm)==3 % If the user clicked
            % Define the inside of the polygone
            [area_j,area_i]=meshgrid(1:M,1:length(rtf(:,1)));

            for k=1:nb_dir-1
                tmp=line(k,mm).getPosition;
                tmp=round([tmp(:,1)*Fe_w tmp(:,2)*M1/Fe_w]);
                points(mm,:,k)=flip(tmp(1,:));
                points(mm,:,k+1)=flip(tmp(2,:));
                u_dir=points(mm,:,k+1)-points(mm,:,k);  % direction vector : points(mm,:,k+1)-points(mm,:,k)
                n_vec=[-u_dir(2) u_dir(1)];             % n_vec=(-dir_y,dir_x)

                proj=(area_i-points(mm,1,k))*n_vec(1)+(area_j-points(mm,2,k))*n_vec(2);

                ind=find(proj<0);
                [area_i,area_j]=deal(area_i(ind),area_j(ind));
            end
            % Save and apply the mask
            filt.('mask')(mm,unique(area_i),unique(area_j))=1;
            
        else % The user pressed a key
            title('Draw a polygon')
            BW = roipoly ;
            filt.('mask')(mm,:,:)=BW;
        end

        mode_rtf_warp=squeeze(filt.('mask')(mm,:,:)).*rtf;
        
        % We don't forget the window normalization factor to calculate the time...
        % signal from the tf response, as well as the *2 factor because the
        % selection is done only on the positive frequencies:
        [mode_temp_warp]=real(sum(mode_rtf_warp,1))/M1/max(hamming(Nww)/norm(hamming(Nww)));

        modes(mm,:)=iwarp_temp_exa(mode_temp_warp,Fe_w,r,c,Fe,N);
        disp(strcat('Mode ',int2str(mm),' filtered'))
        clear mode_rtf_warp mode_temp_warp
    end

    else
        Nmode=0;
        modes=[];
end
disp('End of warping')

end

% Link fonctions
function link1(h,pos,nn,mm)
    pos1=h(nn,mm).getPosition();
    h(nn,mm).setPosition([pos(2,:);pos1(2,:)]);
end
function link2(h,pos,nn,mm)
    pos1=h(nn,mm).getPosition();
    h(nn,mm).setPosition([pos1(1,:);pos(1,:)]);
end