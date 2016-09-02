
clear ;close all
strr='abcdefg';
for Isite=1:5
    
    plotlabel=sprintf('Shell07_Site%i_NeuralNetupdated',Isite);
    datadir=sprintf('/Users/khkim/Thode/2007_Arctic_Analysis/Processed2007/Site_0%i/%s/morph',Isite,plotlabel);
    searchstr='Filter';
    
    if Isite==1
        Nstation=1:7;
    else
        Nstation=[1:7 ];
    end
    %Nstation=7;
    
    [N{Isite},Nplot{Isite},Ndiff{Isite},tdate{Isite}]=plot_automated_stage_statistics(datadir,searchstr,plotlabel,Nstation);
    
    save FigureAutomatedStatistics_results
    
    if 1==0
        subplot(5,1,Isite)
        hh=area(tdate{Isite},Ndiff{Isite}'/1000);
         datetick('x',6);
       
        set(gca,'fontweight','bold','fontsize',14);
        plot_letter_label([strr(Isite) ')']);
        grid on;
        set(gca,'layer','top');   %Put grid lines on top
        title(plotlabel);
        if Isite==5
            legend(fliplr(hh),{'interval removal stage','image processing stage','neural net stage','linking stage','localization stage','final result'})
            xlabel('Date')
        
        end
    end
end

diary off

close all
for Isite=1:5
    
    
    subplot(5,1,Isite)
    hh=area(tdate{Isite},Ndiff{Isite}'/1000);
    
    if Isite==1
        xlimm=xlim;
        xtick=get(gca,'xtick');
    else
       xlim(xlimm)
       %set(gca,'xtick',xtick);
    end
      datetick('x',6,'keeplimits');
  
    set(gca,'fontweight','bold');
    plot_letter_label([strr(Isite) ')']);
    grid on;
    set(gca,'layer','top');   %Put grid lines on top
    title(plotlabel);
    ylabel('Detections (thousands)')
    if Isite==5
        %legend(fliplr(hh),{'interval removal stage','image processing stage','neural net stage','linking stage','localization stage','final result'})
        xlabel('Date')
        
    end
    
end

