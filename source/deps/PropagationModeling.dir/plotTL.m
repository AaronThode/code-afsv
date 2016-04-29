%%%%compute_and_plot_TL%%%%
repeat=1;
while repeat==1
    
    if flag_suppress_io==0
        %yess=menu('TL slice:','range-frequency','range_depth','power-law fit');
        yess=2;
    else
        yess=2;
    end
    figure;
    if yess==1  %%Range-frequency plot
        disp(zplot);
        if flag_suppress_io==0
            zwant=input('Enter desired depth:');
            if isempty(zwant)
                zwant=max(zplot);
            end
        else
            zwant=max(zplot);
        end
        [junk,Iwant]=min(abs(zwant-zplot));
        TLdB=squeeze(TLdB0(Iwant,:,:));
        %imagesc(rplot/1000,freq,TLdB');axis('xy')
        if run_options.incoherent_addition==1
            subplot(2,1,1)
            [cc,hh]=contourf(rplot/1000,freq,TLdB',20:5:140);
            
            clabel(cc,hh,'fontsize',14,'color','k','rotation',0)
            set(gca,'yscale','linear');
            %set(gca,'ytick',[50 100 200 400 800 1600 3200])
            grid on
            xlabel('range(km)','fontweight','bold','fontsize',14);ylabel('frequency (Hz)','fontweight','bold','fontsize',14);
            set(gca,'fontweight','bold','fontsize',14)
            plotname=sprintf('Source depth: %6.2f m, receiver depth: %6.2f, water depth: %6.2f, bottom type: %s, ', ...
                sd,zplot(Iwant),D,case_bottom');
            title(plotname);
            
            mfig=gcf;
            prop_law=compute_simple_TL_law(rplot,TLdB,freq,nfreq);
            figure(mfig);
            subplot(2,1,2)
            hist(prop_law/10);
            xlabel('Modeled alpha');
            title(sprintf('Distribution of alpha between %6.2f and %6.2f Hz',min(freq),max(freq)));
            grid on
            
            save TL_result rplot freq TLdB sd zplot Iwant D case_bottom nfreq
        else
            subplot(2,1,1)
            imagesc(rplot/1000,[],TLdB');axis('xy')
            set(gca,'ytick',1:length(freq),'yticklabel',round(freq));
            %colormap(flipud(jet));
            set(gca,'fontweight','bold');
            %title(num2str(rd));
            %caxis([40 120]);
            caxis([20 100]);
            colorbar;
            
            grid on
            xlabel('range(km)','fontweight','bold','fontsize',14);ylabel('frequency (Hz)','fontweight','bold','fontsize',14);
            set(gca,'fontweight','bold','fontsize',14)
            plotname=sprintf('%s \n%s, Source depth: %6.2f m, receiver depth: %6.2f, water depth: %6.2f, bottom type: %s ', ...
                savename, code_chc,sd,zplot(Iwant),D,case_bottom');
            title(plotname,'interp','none');
            
            mfig=gcf;
            [prop_law, prop_law_mean]=compute_simple_TL_law(rplot,TLdB,freq,nfreq);
            figure(mfig);
            subplot(2,1,2)
            hist(prop_law/10);
            xlabel('Modeled alpha');
            title(sprintf('Distribution of alpha between %6.2f and %6.2f Hz, arithmetic value %6.2f',min(freq),max(freq),prop_law_mean));
            grid on
            
        end
        
        
        %plotname=sprintf('%s_source%5.2fmRcvr%5.2fmDepth%5.2fmBottom%s_Incoh%i.jpg',case_output,sd,zplot(Iwant),D,case_bottom,run_options.incoherent_addition);
        figure(1);
        print('-djpeg',sprintf('%s_%6.2f.jpg',savename,pwavatten(2)))
        
    elseif yess==2  %Range-depth plot
        disp(freq);
        if length(freq)==1
            Iwant=1;
        elseif flag_suppress_io==0
            fwant=input('Enter desired frequency:');
            [junk,Iwant]=min(abs(freq-fwant));
        else
            fwant=max(freq);
            Iwant=length(freq);
        end
        TLdB=squeeze(TLdB0(:,:,Iwant));
        
        %imagesc(rplot/1000,zplot,TLdB);
        if size(TLdB,1)==1
            plot(rplot/1000,TLdB,'k');
            set(gca,'fontweight','bold','fontsize',14)
            xlabel('range(km)','fontweight','bold','fontsize',14);
            ylabel('TL (dB)','fontweight','bold','fontsize',14);axis('ij');grid on
            
        else
            if run_options.incoherent_addition==1
                [cc,hh]=contourf(rplot/1000,zplot,TLdB,10:5:140);axis('ij')
                clabel(cc,hh,'fontsize',14,'color','k','rotation',0);
            else
                imagesc(rplot/1000,zplot,TLdB);axis('ij');caxis([10 140]);colorbar;grid on
            end
            colormap(flipud(jet))
            xlabel('range(km)','fontweight','bold','fontsize',14);ylabel('depth (m)','fontweight','bold','fontsize',14);
            set(gca,'fontweight','bold','fontsize',14)
            hold on;
            if ~isempty(bath)
                if size(bath,1)==1
                    bath=[bath;max(rplot/1000) bath(1,2)];
                end
                line(bath(:,1),bath(:,2));
            end
            
            
        end
        %title(sprintf('%s, Source depth: %6.2f m, frequency: %6.2f Hz, bottom type/depth: %s/%i m, azimuth %i model:%s', ...
        %     case_output,sd,freq(Iwant),case_bottom',D, azi(Iazi),code_chc));
        title(sprintf('%s, Source depth: %6.2f m, D: %6.2f m,%6.2f Hz,%s', ...
            case_output,sd,D, freq(Iwant),code_chc));
        
        plotno=input('Enter number for plot:');
        if ~isempty(plotno)
            orient landscape
            print('-djpeg',sprintf('TLplot_%i',plotno));
        end
        
    else %TL fit plot
       % disp(zplot);
        prop_law=zeros(length(zplot),nfreq);
        for Iwant=1:length(zplot)
            TLdB=squeeze(TLdB0(Iwant,:,:));
            %imagesc(rplot/1000,freq,TLdB');axis('xy')
            %disp('Max range modeled is 
            [prop_law(Iwant,:), prop_law_mean]=compute_simple_TL_law(rplot,TLdB,freq,nfreq);
            
            if 1==0
                subplot(2,1,1)
                imagesc(rplot/1000,[],TLdB');axis('xy')
                set(gca,'ytick',1:length(freq),'yticklabel',round(freq));
                %colormap(flipud(jet));
                set(gca,'fontweight','bold')
                %title(num2str(rd));
                %caxis([40 120]);
                caxis([20 100]);
                colorbar;
                
                grid on
                xlabel('range(km)','fontweight','bold','fontsize',14);ylabel('frequency (Hz)','fontweight','bold','fontsize',14);
                set(gca,'fontweight','bold','fontsize',14)
                plotname=sprintf('%s \n%s, Source depth: %6.2f m, receiver depth: %6.2f, water depth: %6.2f, bottom type: %s ', ...
                    savename, code_chc,sd,zplot(Iwant),D,case_bottom');
                title(plotname,'interp','none');
                
                mfig=gcf;
                figure(mfig);
                subplot(2,1,2)
                hist(prop_law(Iwant,:));
                xlabel('Modeled alpha');
                title(sprintf('Distribution of alpha between %6.2f and %6.2f Hz, depth: %6.2f, arithmetic value %6.2f',min(freq),max(freq),zplot(Iwant),prop_law_mean));
                grid on
                pause(0.25);
            end
            
        end
        
        %plotname=sprintf('%s_source%5.2fmRcvr%5.2fmDepth%5.2fmBottom%s_Incoh%i.jpg',case_output,sd,zplot(Iwant),D,case_bottom,run_options.incoherent_addition);
        figure(2);
       % pcolor(freq,zplot,prop_law);shading flat; caxis([15 25]);colormap(jet);colorbar
        [cc,hh]=contourf(freq,zplot,prop_law,8:1:20);
        clabel(cc,hh);colorbar
        axis('ij')
        title(sprintf('Power-law coefficient fit to coherent transmission loss: %s',savename));
        xlabel('Frequency(Hz)')
        ylabel('Source depth (m)');
        orient landscape
        text(100,50,sprintf('Bottom absorption %6.2f dB/lambda',pwavatten(2)));
        print('-djpeg',sprintf('%s_%3.2fatten.jpg',savename,pwavatten(2)));
        TLdB0=uint16(-20*log10(abs(TL(:,:,:,Iazi))+eps));
            
        save(sprintf('%s_%3.2fatten.mat',savename,pwavatten(2)),'TLdB0','freq','zplot','prop_law','ro','pwavatten','savename');
        
    end  %if yes==1
    
    %colormap(flipud(jet));
    set(gca,'fontweight','bold')
    %title(num2str(rd));
    %caxis([40 120]);
    %caxis([20 90]);
    colorbar;
    
    %title(sprintf('Average of %i frequencies between %i and %i Hz',length(freq),freq(1),freq(end)));
    if flag_suppress_io==0
        repeat=menu('Repeat?','Yes','No');
    else
        repeat=2;
    end
end
