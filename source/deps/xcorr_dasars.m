function [tdoa_mat]=xcorr_dasars(xsample,tsec,DASAR_coords)
% cross-correlates signals between all DASAR pairs

% default number of DASARs
Nd=size(DASAR_coords,1);

tdoa_mat=nan(Nd,Nd);
for I=1:Nd
    for II=1:Nd
        
        if isempty(xsample{I})==1 || isempty(xsample{II})==1
            tdoa_mat(I,II)=nan; % set to NaN if no data 
        elseif I==II 
            tdoa_mat(I,II)=0;   % set to 0 when cross-correlating same DASARs
        else        
            [cc,lags] = xcov(xsample{I}(:,1),xsample{II}(:,1));
            [~,Ilag]=max(abs(hilbert(cc)));
            tau=lags(Ilag)/1000;
            tau_corrected=tsec(I)-tsec(II)+tau; % adjust according to start time of boxes
            
            % populates matrix 
            tdoa_mat(I,II)=tau_corrected;
        end
        
        
    end
end