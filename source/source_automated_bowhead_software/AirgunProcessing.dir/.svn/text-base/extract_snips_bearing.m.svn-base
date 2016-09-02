%%%%%%extract_snips_bearing.m%%%%%%%%
% function thet=extract_snips_bearing(data_all,run_options,param,short_fname)
%
%  Extract a time series from a JAVA CFAR snips file and compute its
%    bearing.
%  Input:
%       data_all: structure containing 'features', a [Nfeatures,Ndetection]
%       array, output of readEnergySummary.
function thet=extract_snips_bearing(data_all,run_options,param,short_fname)
fclose('all');
Igood=1:length(data_all.ctime);
thet=-1*ones(size(Igood));
Iwant=[];
nnames={'min_freq'  'max_freq'};

%%Determine what indicies contain minimum and maximum frequencies of each
%%  detection
for JJ=1:length(nnames),
    for KK=1:length(data_all.names),
        if strcmp(data_all.names{KK},nnames{JJ})
            Iwant=[Iwant KK];
        end
    end
end

%%%If 2007 data, preload nonlinear calibration information.
 if ~isempty(findstr('T2007',short_fname))
        cal07flag=1;
        brefa_table=calibrate_bearing_Shell2007(param.calibration_dir,short_fname);
    else
        cal07flag=0;
 end
 
  
    
%%To prevent loading all data into RAM, load only Ncalls_to_sample at a
%%time
for I=1:ceil(length(Igood)/run_options.Ncalls_to_sample)  
    Iabs=run_options.Ncalls_to_sample*(I-1)+(1:run_options.Ncalls_to_sample);
    Iabs=Iabs(Iabs<=length(Igood));  %Ensure we don't run over end of file
    Iref=Igood(Iabs);  %Indicies to extract from snips file
    snips_name=dir([param.energy.dir_out '/' short_fname '*.snips']);
    
    [x,nstarts_snips,npts_snips,feq,eq,Ireturn,head]=readEnergySnips([param.energy.dir_out '/' snips_name.name], Iref,'double','cell','keep_open');
   
    %Find min and max frequencies of each snip detection
    fmin=data_all.features(Iwant(1),Iref);
    fmax=data_all.features(Iwant(2),Iref);
    
   
    
    for II=1:length(Iref)
        [thet0,kappa,sd]=extract_bearings(x{II}',param.energy.bufferTime,param.Nfft,param.Fs,fmin(II),fmax(II),run_options.bearing_alg,2);
        
        %%Convert acoustic azimuth to geographic azimuth using data loaded
        %%from the header...
        if cal07flag==0
            thet(Iref(II))=bnorm(thet0+head.brefa);
        else
            [junk,Icol]=calibrate_bearing_Shell2007(param.calibration_dir,short_fname,1);
            
            thet(Iref(II))= bnorm(interp1(0:360,brefa_table(:,Icol),bnorm(thet0)));
        end
        
        
        
         
    end
    
    
end

[x,nstarts_snips,npts_snips,feq,eq,Ireturn,head]=readEnergySnips([param.energy.dir_out '/' snips_name.name], 1,'double','cell');
   
end


