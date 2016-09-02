for Isite=0:5
    clear metric param
    load(sprintf('Shell2012_shipping_metrics_combined_Site%i',Isite));
    disp(sprintf('Shell2012_shipping_metrics_combined_Site%i',Isite));
    
    
    for I=1:length(metric{end})
        if isempty(metric{end}{I})
            continue
        end
        try
            metric{end}{I}=rmfield(metric{end}{I},'isi');
            metric{end}{I}=rmfield(metric{end}{I},'entropy');
            metric{end}{I}=rmfield(metric{end}{I},'kurtosis');
            
            metric{end}{I}.adp=metric{end}{I}.adp(1);
            metric{end}{I}.adp{1}=rmfield(metric{end}{I}.adp{1},'med');
            metric{end}{I}.adp{1}=rmfield(metric{end}{I}.adp{1},'pwr');
            metric{end}{I}.adp{1}=rmfield(metric{end}{I}.adp{1},'Fpeaks');
        catch
            fprintf('Isite: %i, DASAR: %i,  isi, entropy, or kurtosis do not exist\n',Isite,I);
        end
    end
    
    %15 MB vs 132 MB at this point%
    
    %This only saves 3 MB
    %for I=1:length(metric{end})
    %   metric{end}{I}.adp{1}.avg.Fpeaks=single(metric{end}{I}.adp{1}.avg.Fpeaks);
    %end
    
    
    save('-v7.3',sprintf('Shell2012_shipping_metrics_basic_Site%i',Isite), 'metric', 'param');
end
