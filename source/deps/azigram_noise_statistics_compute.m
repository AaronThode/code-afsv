function output=azigram_noise_statistics_compute(x,Fs,Nfft,ovlap,params)
%%%Inputs:x:  cell matrix with x{Idasar} being a 3 by N column matrix from
%%%         DASAR Idasar
%%%%Parameters
%params.tabs_start:  %%Start time of x in datenumber.
%params.threshold=5*4;  %Hz
%params.f_transition;  %Frequency at which DASAR p and v go out of phase.
%params.vertical_line=round(params.vertical_line/dF);
%params.min_TBW=round(params.min_TBW/(dT*dF));
%params.brefa=vector the same size as x;
%params.frange
%params.debug
%params.azi_grid, da
%params.threshold
%params.time_chunk=10*60;  %How much data to read into memory at a
%           time--solves clock drift problem too
%
% Outputs:
%  (Ipp,It,If,var):  dimensions of output histograms
%  PdB_percentiles(It,If,Ipp)

Npercentiles=length(params.percentiles)-1;

%%%Determine how big a chunk to import

params.time_chunk=min([params.time_chunk size(x,2)/Fs]);
%Nt=floor(params.time_chunk/dT);
%Npt=floor(1+((size(x,2)/Nfft)-1)/(1-ovlap));
Iwindow=unique([1:(params.time_chunk*Fs):size(x,2) size(x,2)]);  %Index of chunks

if params.time_chunk==size(x,2)/Fs  %%%Chunk is whole window
    Iwindow=[1 size(x,2)+1];
end

dT=Nfft*(1-ovlap);
dF=Fs/Nfft;

psd_fac=2./(Fs.*sum(hanning(Nfft).^2));  %converts FFT into one-sided power spectral density
Ntotal_window=(length(Iwindow)-1);

for It=1:Ntotal_window
    %fprintf('Processing chunk %i of %i\n',It,length(Iwindow)-1);
    Iindex=Iwindow(It):(Iwindow(It+1)-1);
    
    
    [TT,FF,azi,PdB,~,Ix,Iy]=compute_directional_metrics ...
        (x(:,Iindex),{'Directionality','Directionality','ItoERatio','ItoERatio', ...
        'KEtoPERatio','KEtoPERatio','IntensityPhase','IntensityPhase'},Fs,Nfft,ovlap,params, ...
        'gsi',[false true false true false true false true]);
    [~,Icut]=min(abs(FF-params.f_transition));
    
    %     [~,~,temp,~,~,Ixr,Iyr]=compute_directional_metrics ...
    %         (x,'Directionality',Fs,Nfft,ovlap,params_temp, ...
    %         'gsi',true);
    
    Iff=find(FF>=params.frange(1)&FF<=params.frange(2));
    
    %%%Correct for reactive/active intensity flip around 300 Hz...
    azi{1}((Icut+1):end,:)=azi{2}((Icut+1):end,:);
    Ix((Icut+1):end,:)=1i*Ix((Icut+1):end,:);  %Shift 90 degrees
    Iy((Icut+1):end,:)=1i*Iy((Icut+1):end,:);  %Shift 90 degrees
    
    ItoE=azi{3};  %Transport ratio
    ItoE((Icut+1):end,:)=azi{4}((Icut+1):end,:);  %%Correct for reactive error
    
    KEtoPE=azi{5};  %Transport ratio
    KEtoPE((Icut+1):end,:)=azi{6}((Icut+1):end,:);  %%Correct for reactive error
    
    IntensityPhase=azi{7};
    IntensityPhase((Icut+1):end,:)=azi{8}((Icut+1):end,:);  %%Correct for reactive error
    
    azi=azi{1};
    azi=azi(Iff,:);
    
    ItoE=ItoE(Iff,:);
    KEtoPE=KEtoPE(Iff,:);
    Ix=Ix(Iff,:);
    Iy=Iy(Iff,:);
    PdB=PdB(Iff,:);
    IntensityPhase=IntensityPhase(Iff,:);
    
    
    FF=FF(Iff);  %Restricted spectrum
    
    if It==1  %Initialize ouputs, now that we have FF
        output.PdB_percentiles=zeros(Ntotal_window,length(FF), Npercentiles+1);
        output.PdB_hist=zeros(Npercentiles,Ntotal_window,length(FF), length(params.PdB_grid)-1);
        output.azi_hist=zeros(Npercentiles,Ntotal_window,length(FF), length(params.azi_grid)-1);
        output.ItoE_hist=zeros(Npercentiles,Ntotal_window,length(FF), length(params.ItoE_grid)-1);
        output.KEtoPE_hist=zeros(Npercentiles,Ntotal_window,length(FF), length(params.KEtoPE_grid)-1);
        output.IntensityPhase_hist=zeros(Npercentiles,Ntotal_window,length(FF), length(params.IntensityPhase_grid)-1);
        
    end
    
    if params.debug.image
        figure(5);
        plot_azigram_overview;
        keyboard;
        close(5);
    end
    
    %%%%Compute azigram difference
    %%%Convert parameter values into pixels
    dF=FF(2)-FF(1);
    dT=TT(2)-TT(1);
    Tabs=params.tabs_start+datenum(0,0,0,0,0,Iwindow(It)/Fs+TT);
    
    Nffts_per_avg=size(azi,2);  %%%Number of FFT samples within a time window
    Iper=floor(params.percentiles*Nffts_per_avg);
    
    if Iper(1)==0
        Iper(1)=1;
    end
    if Iper(end)>Nffts_per_avg
        Iper(end)=Nffts_per_avg;
    end
    
    
    [PdB_sort,Isort]=sort(PdB,2);
    
    
    output.PdB_percentiles(It,:,:)=PdB_sort(:,Iper);
    %output.results{It}.Tabs=Tabs;
    output.tabs(It)=Tabs(1);
    
    %plot(output{1}.PdB_percentiles,FF)
    
    
    for If=1:length(FF)
        for Ipp=1:Npercentiles
            Iwant=Isort(If,Iper(Ipp):Iper(Ipp+1));
            output.PdB_hist(Ipp,It,If,:)=histcounts(PdB(If,Iwant),params.PdB_grid);
            output.azi_hist(Ipp,It,If,:)=histcounts(azi(If,Iwant),params.azi_grid);
            output.ItoE_hist(Ipp,It,If,:)=histcounts(ItoE(If,Iwant),params.ItoE_grid);
            output.KEtoPE_hist(Ipp,It,If,:)=histcounts(log10(KEtoPE(If,Iwant)),params.KEtoPE_grid);
            output.IntensityPhase_hist(Ipp,It,If,:)=histcounts(IntensityPhase(If,Iwant),params.IntensityPhase_grid);
        end
    end
    
    %subplot(1,3,1);imagesc(params.azi_grid,FF,output.results{It}.azi_hist);
    %subplot(1,3,2);imagesc(params.ItoE_grid,FF,output.results{It}.ItoE_hist);
    %subplot(1,3,3);imagesc(params.KEtoPE_grid,FF,output.results{It}.KEtoPE);
    
end  %It  window


output.params=params;
output.FF=FF;



