%%%%compute_replicas.m%%%%%%%%%%%
%[TL,rplot,zplot,sd,kr,Umode]=compute_replicas(code_chc,case_bottom,freq,nfreq,sspopt,lthick,cmedtop,cmedbot, rho, pwavatten, ...
%    svp,bath,ro,rd,sd,dr_factor,dz_factor,fu,rms_roughness,run_options)%% Output:
%%  TL: transmission loss, dimensions of  [length(rd),length(ro)];
%
% Input:
%       run_options:
%           nofield: do not compute field for KRAKEN
%           incoherent_addition: compute TL incoherently
%
%%          case_output specifies type of acoustic output.  Categoreis include
%                TL_C:  coherent transmission loss
%                TL_I:  incoherent transmission loss
%                TL_S:  semicoherent transmission loss (Lloyd mirror source pattern)-BELLHOP
%                time_series:  time series
%                group_velocity:  group velocity curves
%                ray:  ray tracing.


function [TL,rplot,zplot,sd,Umode]=compute_replicas(code_chc,case_bottom,freq,sspopt,lthick,cmedtop,cmedbot, rho, pwavatten, ...
    svp,bath,ro,rd,sd,dr_factor,dz_factor,fu,rms_roughness,run_options)

if ~isempty(dir('inp_file*'))
   !rm inp_file.*
end
TL=[];rplot=[];
zplot=[];Umode=[];
%%%**  START FREQUENCY LOOP ************

incoherent_addition=run_options.incoherent_addition;
nofield=run_options.nofield;
case_output=run_options.case_output;

compute_own_field=1;  %If one, compute field in MATLAB using modes...
Umode=[];
to=clock;
kr=[];
rmax=max(ro);
zmax=max(rd);  %Maximum receiver depth
%dr=0.1*1500/fu;
%dz=0.05*dr;
%disp('Warning, using coarse RAM step');
%%Get ranges and depths of grid sizes...
dr=dr_factor*1500/fu;
dz=dz_factor*dr;
Nmesh=[round(lthick/dz)];

if strcmp(code_chc,'RAM'),
    nr=length(ro);
    nz=length(rd);
    
    %covert bath variable ranges to units of meters
    try
        bath(:,1)=bath(:,1)*1000;
    end
    [ro_grid,rd_grid,ndr,ndz]=writeRAM(sd,rd(1),dz,dr,nr,nz,zmax, ...
        rmax,case_bottom,freq(1),lthick, cmedtop, cmedbot, rho, pwavatten, svp,bath);
    nr=length(ro_grid);
    nz=length(rd_grid);
    if size(rd_grid,2)>1,
        rd_grid=rd_grid';
    end
    % keyboard;
end


%%Interpolate sound speed profile to estimate small horizontal range
%%offsets...assumes phones are always in the water column


if strcmpi(code_chc,'RAM')
    
    disp(sprintf('RAM frequency: %6.2f Hz',freq(1)));
    [ro_grid,rd_grid,ndr,ndz]=writeRAM(sd,rd(1),dz,dr,nr,nz,zmax, ...
        rmax,case_bottom,freq,lthick, cmedtop, cmedbot, rho, pwavatten, svp,bath);
    %% Most recent version of RAM
    %!./ramgeo_exec_ibm
    % What was working in Feb 04
    mypath=which('compute_replicas.m');
    Islash=findstr(mypath,'/');
    mypath=mypath(1:max(Islash));
    eval(sprintf('!%s/ramgeo1.5_aaron',mypath));
    %disp('finished execution');
    if length(rd)>1
        [p,TLL,rplot,zplot,nz,dz,nzplt,zmplt,rmax,dr,ndr]=readRAMgrd;
        if size(rplot,1)>1,rplot=rplot';end
        if size(zplot,2)>1,zplot=zplot';end
        
    else
        [p,rplot]=readRAM;
        zplot=rd(1);
        
    end
    TL=p.*(ones(size(zplot))*exp(i*2*pi*freq*rplot./1500));
    
    %keyboard;
elseif ~isempty(findstr(lower(code_chc),lower('KRAKEN')))
    
    if isempty(bath)
        %bath=[0 D;
        %    max(r0/1000) D];
        
        bath=[0 D];
    elseif size(bath,1)>1
        disp('WARNING, range-dependent bathymetry not reliable in KRAKEN');
        
    end
    for Id=1:size(bath,1)
        if ~isempty(dir('inp_file*'))
            !rm inp_file.*
        end
        D=bath(Id,2);
        lthick(1)=D;
        if ~isempty(svp)
            Igood=find(svp(:,1)<=D);
            svp=svp(Igood,:);
            svp(end,1)=bath(Id,2);
        end
        
        Nmesh=ceil(50*lthick./(1500/max(freq)));
        Nmesh=max([Nmesh;100*ones(1,length(Nmesh))]);  %Ensure some minimum number of points...
        writeKRAKEN(sd,rd,max(ro/1000),'inp_file',freq, Nmesh, lthick, cmedtop, ...
            cmedbot, rho, pwavatten,svp,[],[],rms_roughness)
        if strcmpi(code_chc,'KRAKEN'),
            !kraken.exe inp_file
        else
            !krakenc.exe inp_file
        end
        
        try
            [kr,ztab,phi]=readm('inp_file');
            
            %Check normalization
            %dz=ztab(2)-ztab(1);
           % dz*trapz(phi.^2)
           
            Nkr=min([run_options.max_kr length(kr)]);
            Umode{Id}.z=ztab;
            Umode{Id}.phi=phi;
            Umode{Id}.kr=kr(1:Nkr);
            %disp(sprintf('Number modes: %i',length(kr)));
            kraken_fail=false;
            
        catch
            kraken_fail=true;
            ztab=[];
            Umode{Id}.phi=[];
            Umode{Id}.kr=[];
            Umode{Id}.z=ztab;
        end
        if kraken_fail||nofield==1 || size(bath,1)>1  %If rd bathymetry
            TL=[];zplot=ztab;rplot=[];
            continue
        end
        
        if compute_own_field==0
            writeflp(0,ro/1000,sd,rd,'inp_file',zeros(size(rd)),9999,incoherent_addition)
            !field inp_file > field_outout.txt
            %[TL,rplot,sd,zplot,totime]=readshd_corrected('inp_file', 1);
            try
                [title,freq,nsd,nrd,nrr,sd,zplot,rplot,TL]=read_shd_bin( 'inp_file.shd' );
            catch
                TL=[];zplot=[];rplot=[];
            end
            if size(rplot,1)>1,rplot=rplot';end
            
        else
            for Iz=1:length(sd)
                [~,Is]=min(abs(Umode{Id}.z-sd(Iz)));
                %a = squeeze(model.U( isd, Iwant(If), : ));
                %phi = squeeze(model.U(:,Iwant(If),:)) * diag( a, 0 );	% scale modes by a
                
                a=phi(Is,:);
                phis=phi*diag(a,0);
                
                % ******************************************************
                % form pressure field
                % ******************************************************
                
                phase = diag( 1.0 ./ sqrt( Umode{Id}.kr ) ) * exp( 1i * Umode{Id}.kr * ro ) * diag( sqrt( 2 * pi ./ ro ) );
                
                Ibad = find(isnan(real(phase)));
                phase(Ibad)=0;
                
                TL = 1i*exp(-1i*pi/4)*phis * phase;  % [rd Nm] x [Nm ranges] = [rd ranges]
                zplot=rd;
                rplot=ro;
                
            end
            
            %             figure
            %             subplot(2,1,1)
            %             imagesc(20*log10(abs(TL)))
            %             subplot(2,1,2)
            %             imagesc(20*log10(abs(TL2)))
            %
            %              figure
            %             subplot(2,1,1)
            %             imagesc(angle(TL))
            %             subplot(2,1,2)
            %             imagesc(angle(TL2))
            %keyboard
        end %compute_own_fiel
    end  %bathymetry step
elseif strcmpi(code_chc,'scooter')
    writeKRAKEN(sd,rd,max(ro/1000),'inp_file',freq, Nmesh, lthick, cmedtop, ...
        cmedbot, rho, pwavatten,svp)
    !scooter inp_file
    writeflps(ro/1000,rd,'inp_file')
    !fields inp_file
    [TL,rplot,sd,zplot,totime]=readshd('inp_file', 1);
    if size(rplot,1)>1,rplot=rplot';end
    
elseif strcmpi(code_chc,'bellhop')
    %         da=sqrt(1500/(6*max(ro)*freq));
    %        	bellhop.nbeams=floor(3+(178*pi/180)/da);
    %
    
    %%Attempt to limit angular spread
    %crit_angle=(min(round(acos(cmedtop(1)/cmedtop(2))*180/pi)+5,89));
    crit_angle=89;
    da=sqrt(1500/(6*max(ro)*freq));
    
    %Compute number of beams
    bellhop.nbeams=floor(2*(pi/180)*crit_angle/da);
    %bellhop.nbeams=0;
    angle_vec=[crit_angle*[-1 1] bellhop.nbeams];
    
    
    %%Set ray option correctly, if this is output desired
    if strcmpi(case_output,'ray')
        case_output='TL_RG';
        prompt1={'ray angles desired (examples: [0 1 5], -90:90,linspace(0,10,5)):'};
        def1={'[-89:89]'};
        dlgTitle1='Ray trace parameters';
        answer=inputdlg(prompt1,dlgTitle1,1,def1);
        
        %crit_angle=str2num(answer{1});
        %bellhop.nbeams=str2num(answer{2});
        %angle_vec=[crit_angle*[-1 1] bellhop.nbeams];
        angle_vec=eval(answer{1});
        
    end
    
    Iscore=findstr(case_output,'_')+1;
    run_type=case_output(Iscore:end);
    if isempty(run_type)
        run_type='CB';
    elseif length(run_type)==1
        run_type=[run_type 'B']; %B stands for Gaussian beams
    end
    
    %        run_type='aB' ;    % amplitude-delay file with Gaussian beams
    
    
    
    if ~isempty(bath)
        if size(bath,1)==1
            disp('bathymetry only specified at one point; assume flat bathymetry and add one more bathymetry point');
            bath=[bath; max(ro)/1000 bath(1,2)];
        end
        writeBELLHOP(sd,rd,ro/1000,'inp_file',freq, sspopt, lthick, cmedtop, svp,run_type,angle_vec,cmedbot, rho,pwavatten,bath);
    else
        writeBELLHOP(sd,rd,ro/1000,'inp_file',freq, sspopt, lthick, cmedtop, svp,run_type,angle_vec,cmedbot, rho,pwavatten);
    end
    !bellhop.exe inp_file
    
    if strcmpi(run_type,'rg') %if a ray trace
        close all
        global units
        units='km';
        plotray('inp_file');
        %[zt,rt,pressure] = plotshd([fname '.shd']);  % For C options
        %[ Arr, Pos ] = read_arrivals_asc( 'inp_file.arr',10 );  %For A
        %    option
        pause;
        
        print('-djpeg','raytrace.jpg');
        return
        
    end
    
    if strcmp(run_type,'AG')|strcmp(run_type,'AB') %if a delay amplitude file
        close all
        [arr,pos]=plot_bellhop_arr('inp_file');
         %[ Arr, Pos ] = read_arrivals_asc( 'inp_file.arr',10 );  %For A
        
        pause;
        
        print('-djpeg','inp_file_bellhop_arr.jpg');
        return
        
    end
    [TL,rplot,sd,zplot,totime]=readshd('inp_file', 1);
    %[ title, freq, nsd, nrd, nrr, sd, rplot, rr, tlt ] = read_shd_bin( 'inp_file.shd' );
    %keyboard;
    if length(zplot)>1,
        Igood=find(diff(zplot)>0);
        Igood=[Igood ; length(zplot)];
        zplot=zplot(Igood);
        TL=TL(Igood,:);
    end
    if size(rplot,1)>1,rplot=rplot';end
elseif strcmpi(code_chc,'Lloyd'),
    rplot=ro;
    zplot=rd.';  if size(zplot,2)>1,zplot=zplot';end
    k=2*pi*freq/1500;
    rd=rd'*ones(1,length(ro));
    ro=ones(length(zplot),1)*rplot;
    R=sqrt(ro.^2+(rd-sd).^2);
    TL=-2*i.*exp(-i*k*R).*sin(k.*sd.*rd./R)./R;  %This is actually the direct pressure field, not TL
    
else
    error('propagation code not correct.')
end

if ~isempty(TL)
    nr=length(rplot);
    nz=length(zplot);
    TL=TL(1:nz,1:nr);
end


