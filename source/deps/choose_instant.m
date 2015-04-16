function [fwlims,clims,xmax,values] = choose_instant(x_deconv,jj,Fe,r,c,Nmode,delay,flims,clims)
                args_max=9;
                values=zeros(2,Nmode);
                
                if nargin==args_max-1
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
                M1=floor(M/2);
                rtf(M1:end,:)=[];                   % Delete negative frequencies
                
                RTF=abs(rtf).^2;
                % Adjust the size of the window and keep the same for each
                % time instant
                if nargin==args_max-1
                    fw_min=0;
                    fw_max=M1-1;
                    while(max(RTF(fw_max,:))/max(RTF(:))<=5e-2)
                        fw_max=fw_max-1;
                    end
                    fw_min=fw_min*Fe_w/M1;
                    fw_max=fw_max*Fe_w/M1;
                    clims=[0 max(RTF(:))];
                else
                    fw_min=flims(1);
                    fw_max=flims(2);
                end

                imagescFun(t_w,f_w,RTF,'ij')
                ylim([fw_min,fw_max])
                xlim([0 xmax])
                hold on,

                % Pekeris cutoff frequencies
                %[c1,c2]=deal(1490,1707);
                [c1,c2]=deal(1439,1673);
                D=55;
                pek_cutoff=c1*c2/2/D/sqrt(c2^2-c1^2);
                df=pek_cutoff*0.4;
                threshold=[1 3]; % thresholds to determine the height and the width of the modes (dB)
                colorMode=['r','w','c','y','g','m'];
                hold on
                
                for mm=1:Nmode
                    [x_rect,y_rect]=modes_rect(RTF,threshold,round(pek_cutoff*(mm-0.5)*M1/Fe_w),ceil(df*M1/Fe_w));
                    plot(x_rect/Fe_w,y_rect*Fe_w/M1,'Color','k')
                    plot(t_w(ceil(M/2):end),linspace(pek_cutoff*(mm-0.5),pek_cutoff*(mm-0.5),M-ceil(M/2)+1),'Color',colorMode(1+mod(mm-1,length(colorMode))))
                    values(:,mm)=evalWarp(RTF,x_rect,y_rect);
                end
                fwlims=[fw_min,fw_max];
                caxis(clims)
                    
end

