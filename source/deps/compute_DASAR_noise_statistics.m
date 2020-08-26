function output=compute_DASAR_noise_statistics(x,Fs,Nfft,ovlap,params)
%%%Inputs:x:  cell matrix with x being a 3 row by Nt column matrix from
%%%         DASAR
%%%%Parameters
%params.tabs_start:  %%Start time of x in datenumber.
%params.brefa=vector the same size as x;
%params.time_chunk=60;  %How much data to read into memory at a
%           time--solves clock drift problem too
%
% Outputs:
%  (var_grid,Ilevel,If,It):  dimensions of output histograms, with var_grid
%           being the grid of the dependent variable. Ilevel=intensity or
%           ItoE ratio, If, frequency, It=time index
%  tabs:  absolute times associated with data

if size(x,1)>3&&size(x,2)==3
    disp('Flipping matrix');
    x=x';
end

%%%Determine how big a chunk to import

time_chunk=min([params.time_chunk size(x,2)/Fs]);
%Nt=floor(params.time_chunk/dT);
%Npt=floor(1+((size(x,2)/Nfft)-1)/(1-ovlap));
Iwindow=unique([1:round(time_chunk*Fs):size(x,2) size(x,2)]);  %Index of chunks

if time_chunk==size(x,2)/Fs  %%%Chunk is entire data
    Iwindow=[1 size(x,2)+1];
end

Ntotal_window=(length(Iwindow)-1);

J_azi=1;
J_ItoE=2;
J_KEtoPE=3;
J_phase=4;

%tabs=zeros(1,Ntotal_window);
tabs=params.tabs_start+datenum(0,0,0,0,0,0.5*time_chunk+Iwindow(1:(end-1))/Fs);
for It=1:Ntotal_window
    fprintf('Processing chunk %i of %i\n',It,length(Iwindow)-1);
    
    Iindex=Iwindow(It):(Iwindow(It+1)-1);
    
    
    [TT,FF,azi,PdB,~,Ix,Iy]=compute_directional_metrics ...
        (x(:,Iindex),{'Directionality','ItoERatio', ...
        'KEtoPERatio','IntensityPhase',},Fs,Nfft,ovlap,params, ...
        'gsi',[false  false  false  false ]);
    
    Iff=find(FF>=params.frange(1)&FF<=params.frange(2));
    
    for In=1:4
        azi{In}=azi{In}(Iff,:);
    end
    Ix=Ix(Iff,:);
    Iy=Iy(Iff,:);
    PdB=PdB(Iff,:);
    FF=FF(Iff);  %Restricted spectrum
    
    
    %  (var_grid1,var_grid2,If,It):  dimensions of output histograms, with
    %  var_grid1 and var_grid2
    %           being the grid of the dependent variable. Ilevel=intensity or
    %           ItoE ratio, If, frequency, It=time index
    
    td.permitted_labels={'Azimuth','ItoE', 'KEtoPE','IntensityPhase','PSD'};
    
    
    precision='single';
    if It==1  %Initialize ouputs, now that we have FF
        
        %output.AziVsPdB=zeros(Nazi,length(params.grid.PdB)-1,length(FF), Ntotal_window,precision);
        %output.AziVsItoE=zeros(Nazi,length(params.grid.ItoE)-1,length(FF), Ntotal_window,precision );
        %output.AziVsKEtoPE=zeros(Nazi,length(params.grid.KEtoPE)-1,length(FF), Ntotal_window,precision );
        %output.AziVsIntensityPhase=zeros(Nazi,length(params.grid.IntensityPhase)-1,length(FF), Ntotal_window,precision );
        %output.PdBVsItoE=zeros(length(params.grid.PdB)-1,length(params.grid.ItoE)-1,length(FF), Ntotal_window,precision );
        mid_grid.azi=0.5*(params.grid.azi(1:(end-1))+params.grid.azi(2:end));
        mid_grid.PdB=0.5*(params.grid.PdB(1:(end-1))+params.grid.PdB(2:end));
        mid_grid.ItoE=0.5*(params.grid.ItoE(1:(end-1))+params.grid.ItoE(2:end));
        mid_grid.KEtoPE=0.5*(params.grid.KEtoPE(1:(end-1))+params.grid.KEtoPE(2:end));
        mid_grid.IntensityPhase=0.5*(params.grid.IntensityPhase(1:(end-1))+params.grid.IntensityPhase(2:end));
        
        output.AziVsPdB=MatrixND([],{mid_grid.azi,mid_grid.PdB,FF,tabs},{'Azimuth','PSD','Frequency','Time'});
        output.AziVsItoE=MatrixND([],{mid_grid.azi,mid_grid.ItoE,FF,tabs},{'Azimuth','ItoE','Frequency','Time'});
        output.AziVsKEtoPE=MatrixND([],{mid_grid.azi,mid_grid.KEtoPE,FF,tabs},{'Azimuth','KEtoPE','Frequency','Time'});
        output.AziVsIntensityPhase=MatrixND([],{mid_grid.azi,mid_grid.IntensityPhase,FF,tabs},{'Azimuth','IntensityPhase','Frequency','Time'});
        output.PdBVsItoE=MatrixND([],{mid_grid.PdB,mid_grid.ItoE,FF,tabs},{'PSD','ItoE','Frequency','Time'});
      
    end
    
    if params.debug.image
        figure(5);
        plot_azigram_overview;
        keyboard;
        close(5);
    end
    
    %%%%Compute azigram difference
    %%%Convert parameter values into pixels
    
    %(var_grid1,var_grid2,If,It)
    
    %ItoE=azi{2};  %Transport ratio
    %KEtoPE=azi{3};  %Transport ratio
    %IntensityPhase=azi{4};
    %azi=azi{1};
    %     J_azi=1;
    % J_ItoE=2;
    % J_KEtoPE=3;
    % J_phase=4;
    
    for If=1:length(FF)
        %Iwant=Isort(If,Iper(Ipp):Iper(Ipp+1));
        
        output.AziVsPdB.N(:,:,If,It)=single(histcounts2(azi{J_azi}(If,:),PdB(If,:),params.grid.azi,params.grid.PdB));
        output.AziVsItoE.N(:,:,If,It)=single(histcounts2(azi{J_azi}(If,:),azi{J_ItoE}(If,:),params.grid.azi,params.grid.ItoE));
        output.AziVsKEtoPE.N(:,:,If,It)=single(histcounts2(azi{J_azi}(If,:),azi{J_KEtoPE}(If,:),params.grid.azi,params.grid.KEtoPE));
        output.AziVsIntensityPhase.N(:,:,If,It)=single(histcounts2(azi{J_azi}(If,:),azi{J_phase}(If,:),params.grid.azi,params.grid.IntensityPhase));
        output.PdBVsItoE.N(:,:,If,It)=single(histcounts2(PdB(If,:),azi{J_ItoE}(If,:),params.grid.PdB,params.grid.ItoE));
        
    end
    
    
end  %It  window

%%%Add extra information to output so it is self contained
output.params=params;
output.params.grid.FF=FF;
output.params.grid.tabs=tabs;


