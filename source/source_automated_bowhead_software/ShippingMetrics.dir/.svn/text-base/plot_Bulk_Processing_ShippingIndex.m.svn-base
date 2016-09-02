%%plot_Bulk_Processing_ShippingIndex.m%%

clear
close all
!rm *jpg
Site_vec=[5];

%load Shell2012_shipping_metrics_combined_Site0_KatMacPro
%Npeaks=length(param.df_search);
Npeaks=1;


chc='avg';
xlimm=[datenum(2012,9,16,0,0,0) datenum(2012,9,19,0,0,0)];


%strr='abcdefghijklmnopqrstuvwxyz';


for I=Site_vec
    close all
    if I==0
        DASAR_str='xyz';
    elseif I==1
        DASAR_str='def';
    else
        DASAR_str='aceg';
    end
    %DASAR_str='aceg';
    
    %%%%%%%Convert desired DASAR string reference into a number..
    for II=1:length(DASAR_str)
        if strcmp(DASAR_str(II),'*')
            DASAR_vec=1:27;
        else
            DASAR_vec(II)=double(lower(DASAR_str(II)))-double('a')+1;
        end
    end
    load(sprintf('ShippingMetrics.dir/Shell2012_shipping_metrics_basic_Site%i.mat',I));  %loads metric; param
    
    if I==0
        I=6;
    end
    for J=1:length(metric{I})
        if all(J~=DASAR_vec)
            continue
        end
        DASAR_str=upper(char(J+double('a')-1));
        if ~isempty(metric{I}{J})
            Igood=find(metric{I}{J}.tabs>xlimm(1)&metric{I}{J}.tabs<xlimm(2)); 
            Igood2=find(Igood<=length(metric{I}{J}.adp{1}.(chc).pwr));
            Igood=Igood(Igood2);
            
            metric{I}{J}.tabs=metric{I}{J}.tabs(Igood);
            try
                metric{I}{J}.entropy=metric{I}{J}.entropy(Igood);
                metric{I}{J}.kurtosis=metric{I}{J}.kurtosis(Igood);
            end
            for III=1:Npeaks
                %metric{I}{J}.isi{III}.pwr=metric{I}{J}.isi{III}.(chc).pwr(Igood);
                %metric{I}{J}.isi{III}.Fpeaks=metric{I}{J}.isi{III}.(chc).Fpeaks(:,Igood);
                metric{I}{J}.adp{III}.pwr=metric{I}{J}.adp{III}.(chc).pwr(Igood);
                metric{I}{J}.adp{III}.Fpeaks=metric{I}{J}.adp{III}.(chc).Fpeaks(:,Igood);
            end
        end
        
        
        %%Plots.
        tabs=metric{I}{J}.tabs;
        try
            entropy=metric{I}{J}.entropy;
            kurtosis=metric{I}{J}.kurtosis;
            yes_entropy_flag=1;
        catch
            disp('entropy does not exist');
            yes_entropy_flag=0;
        end
        for III=1:Npeaks
            
            %Erase previous values...
            adp_Nf{III}=[];
            
            if yes_entropy_flag==1
                isi_Nf{III}=[];
                isi_pwr{III}= metric{I}{J}.isi{III}.pwr;
                isi_F{III}= metric{I}{J}.isi{III}.Fpeaks;
                isi_Nf{III}=sum(isi_F{III}>0);
            end
            %for K=1:length(isi_F{III})
            %isi_Nf{III}(K)=sum(isi_F{III}(:,K)>0);
            %end
            
            adp_pwr{III}= metric{I}{J}.adp{III}.pwr;
            adp_F{III}= metric{I}{J}.adp{III}.Fpeaks;
            adp_Nf{III}=sum(adp_F{III}>0);
            
            %             for K=1:length(adp_F{III})
            %                 adp_Nf{III}(K)=sum(adp_F{III}(:,K)>0);
            %             end
        end
        
        
        
        %Plotting
        
        %%Determine tick label
        tspan=datevec(diff(xlimm));
        if tspan(3)>1
            datetickk=6;
        else
            datetickk=14;
        end
        
        if yes_entropy_flag==1
            figure;set(gcf,'pos',[ 69         647        1515         437]);
            subplot(2,1,1);
            plot(tabs,entropy);
            if ~isempty(xlimm)
                xlim(xlimm);
            end
            datetick('x',datetickk,'keeplimits');
            grid on;
            xlabel('Date');ylabel('Relative entropy of linear PSD');
            title(sprintf('%i%s,%s-%s',I,DASAR_str,datestr(xlimm(1)),datestr(xlimm(2))));
            
            subplot(2,1,2);set(gcf,'pos',[ 69         647        1515         437]);
            plot(tabs,kurtosis);
            if ~isempty(xlimm)
                xlim(xlimm);
            end
            datetick('x',datetickk,'keeplimits');grid on;
            xlabel('Date');ylabel('Kurtosis of linear PSD');
            orient landscape
            print('-djpeg',sprintf('KurtosisEntropy_%i_%s_%s_%s_%s.jpg',I,DASAR_str,chc,datestr(xlimm(1)),datestr(xlimm(2))));
        end
        
        %%Plot peak detections...
        for III=1:Npeaks
            figure;set(gcf,'pos',[ 69         647        1515         437]);
            
            if yes_entropy_flag==1
                subplot(2,1,1)
                plot(tabs,isi_pwr{III},tabs,adp_pwr{III});
                if ~isempty(xlimm)
                    xlim(xlimm);
                end
                datetick('x',datetickk,'keeplimits');
                grid on;
                legend('isi','adp');
                xlabel('Date');ylabel('Band power (dB re 1uPa)');
                title(sprintf('%i%s,%s,%s-%s',I,DASAR_str,chc,datestr(xlimm(1)),datestr(xlimm(2))));
                
                subplot(2,1,2);
                plot(tabs,isi_Nf{III},tabs,adp_Nf{III});
            else
                subplot(2,1,1);
                plot(tabs,adp_pwr{III},'.');
                if ~isempty(xlimm)
                    xlim(xlimm);
                end
                datetick('x',datetickk,'keeplimits');
                grid on;
                legend('adp');
                xlabel('Date');ylabel('Band power (dB re 1uPa)');
                title(sprintf('%i%s,%s,%s-%s',I,DASAR_str,chc,datestr(xlimm(1)),datestr(xlimm(2))));
                
                subplot(2,1,2);
                plot(tabs,adp_Nf{III},'.');
                
            end
            
            if ~isempty(xlimm)
                xlim(xlimm);
            end
            datetick('x',datetickk,'keeplimits');
            grid on;
            xlabel('Date');ylabel('Number of bands');
            title(sprintf('param.df_search: %6.2f',param.df_search(III)));
            
            orient landscape
            print('-djpeg',sprintf('BandSummary_%i_%s_%s_%i_%s_%s.jpg',I,DASAR_str,chc,param.df_search(III),datestr(xlimm(1)),datestr(xlimm(2))));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%plot detected frequencies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure;
            set(gcf,'pos',[ 69         647        1515         437]);
            
            if yes_entropy_flag==1
                subplot(2,1,1)
                plot(tabs,isi_F{III},'.');
                if ~isempty(xlimm)
                    xlim(xlimm);
                end
                datetick('x',datetickk,'keeplimits');
                grid on;
                xlabel('Date');ylabel('isi Frequency (Hz)');
                title(sprintf('%i%s,%s,%s-%s',I,DASAR_str,chc,datestr(xlimm(1)),datestr(xlimm(2))));
                
                subplot(2,1,2)
            end
            
            plot(tabs,adp_F{III},'.');
            if ~isempty(xlimm)
                xlim(xlimm);
            end
            datetick('x',datetickk,'keeplimits');
            grid on;
            xlabel('Date');ylabel('adp Frequency (Hz)');
            title(sprintf('param.df_search: %6.2f',param.df_search(III)));
            orient landscape
            print('-djpeg',sprintf('BandFreq_%i_%s_%s_%i_%s_%s.jpg',I,DASAR_str,chc,param.df_search(III),datestr(xlimm(1)),datestr(xlimm(2))));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Histogram entire deployment...
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure;
            set(gcf,'pos',[ 69         647        1515         437]);
            if yes_entropy_flag==1
                subplot(2,1,1)
                df_size=10;
                [NN,FF]=hist(isi_F{III}(:),0:df_size:450);
                barh(FF(2:end),NN(2:end)/length(tabs));grid on
                grid on;
                xlabel('Fraction of Time');ylabel('isi Frequency (Hz)');
                title(sprintf('%i%s,%s,%s-%s',I,DASAR_str,chc,datestr(xlimm(1)),datestr(xlimm(2))));
                
                subplot(2,1,2)
            end
            [NN,FF]=hist(adp_F{III}(:),0:10:450);
            barh(FF(2:end),NN(2:end)/length(tabs));grid on
            grid on;
            xlabel('Fraction of Time');ylabel('adp Frequency (Hz)');
            title(sprintf('%i%s,%s,%s-%s',I,DASAR_str,chc,datestr(xlimm(1)),datestr(xlimm(2))));
            
            
            orient landscape
            print('-djpeg',sprintf('HistFreq_%i_%s_%s_%i_%s_%s.jpg',I,DASAR_str,chc,param.df_search(III),datestr(xlimm(1)),datestr(xlimm(2))));
            
            
        end
        
    end %J, DASAR review...
    
end  %I, Isite