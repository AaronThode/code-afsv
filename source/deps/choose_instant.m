function [fw_min,fw_max,xmax] = choose_instant(x_deconv,jj,Fe,r,c,delay,flims)

                if nargin==6
                    x_ok=x_deconv(delay:end);
                    [s_w, Fe_w]=warp_temp_exa(x_ok,Fe,r,c);
                    xmax=(length(s_w)-1)/Fe_w;
                else
                    xmax=delay;
                end
                
                x_min=max(1,jj);               
                x_ok=x_deconv(x_min:end);

        %% Warping
                [s_w, Fe_w]=warp_temp_exa(x_ok,Fe,r,c);     % s_w: warped signal, Fe_w: new warping frequency
                clear x_ok
                M=length(s_w);
                t_w=(0:M-1)/Fe_w;                           % Warped time
                f_w=(0:M-1)*Fe_w/M;                         % Warped frequencies
                Nww=1+2*ceil(M/8);                          % Number of points for the tfr of the warped signal
                rtf=tfrstft(s_w,1:M,M,hamming(Nww));
                
                RTF=abs(rtf);
                RTF(floor(M/2):end,:)=[];                   % Delete negative frequencies
                
                % Adjust the size of the window and keep the same for each
                % time instant
                if nargin==6  
                    fw_min=1;
                    fw_max=floor(M/2)-1;
                    while(max(RTF(fw_max,:))/max(RTF(:))<=1e-3)
                        fw_max=fw_max-1;
                    end
                    fw_min=fw_min*Fe_w/M;
                    fw_max=fw_max*Fe_w/M; 
                else
                    fw_min=flims(1);
                    fw_max=flims(2);
                end

                imagesc(t_w,f_w,RTF);
                
                % Pekeris cutoff frequencies
                [c1,c2]=deal(1490,1707);
                D=55;
                pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
                df=10;
                
                imagesc(t_w,f_w,RTF)
                hold on
                for mm=1:3
                    [x_rect,y_rect]=modes_rect(RTF,6,round(pek_cutoff*(mm-0.5)*M/Fe_w),df*M/Fe_w);
                    plot(x_rect/Fe_w,y_rect*Fe_w/M,'r-')
                end
                
                xlim([0 xmax])
                ylim([fw_min fw_max])
                axis ij
                colorbar
                title(['delay : ' num2str(jj)])
end

