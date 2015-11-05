function arr=plotray( rayfil )

% useage: plotray( rayfil )
% where rayfil is the ray file (without the extension)
% e.g. plotray( 'foofoo' )
%
% plots the RAYfil produced by Bellhop
% MBP July 1999

if ( strcmp( rayfil, 'RAYFIL' ) == 0 )
    rayfil = [ rayfil '.ray' ]; % append extension
end

% plots a BELLHOP ray file

zr = 90.0;	% use this to just plot eigenrays

% open the file

fid = fopen( rayfil, 'r' );
if ( fid == -1 )
    warndlg( 'No ray file exists; you must run BELLHOP first (with ray ouput selected)', 'Warning' );
end

% read header stuff

TITLE  = fgetl(  fid );
tmp=fgetl(fid);
FREQ   = sscanf( tmp, '%f', 1 );
NBEAMS = fscanf( fid, '%i', 1 );

DEPTHT = fscanf( fid, '%f', 1 );
tmp=fgetl(fid);
tmp=fgetl(fid);
DEPTHB = sscanf( tmp, '%f', 1 );

ii = findstr( TITLE(3:end), '''');   % find last quote
TITLE = deblank( TITLE(3:1:ii-1) );  % remove whitespace

% read rays

set( gca, 'YDir', 'Reverse','fontweight','bold','fontsize',14 );   % plot with depth-axis positive down

xlabel( 'Range (km)' )
ylabel( 'Depth (m)' )
hold on

% axis limits
rmin = +1e9;
rmax = -1e9;
zmin = +1e9;
zmax = -1e9;

for ibeam = 1:NBEAMS
    tmp=fgetl(fid);
    alpha0 = sscanf( tmp, '%f', 1 );
   % alpha0 = fscanf( fid, '%f', 1 );
    if isempty(alpha0),
        disp('end of line');
        break
    else
        fprintf('Angle: %6.2f\n',alpha0);
  
        arr.src_ang(ibeam)=alpha0;
        tmp=fgetl(fid);
        tmp2=sscanf(tmp,'%i',3);
        arr.nsteps(ibeam) = tmp2(1);
        arr.NumTopBnc(ibeam) = tmp2(2);
        arr.NumBotBnc(ibeam) = tmp2(3);
        if isempty( arr.nsteps(ibeam) ) break; end
        
        %%Problem, the ASCII text file often drops colons into the mix
        Ns=arr.nsteps(ibeam);
        r=zeros(1, Ns);
        z=zeros(1, Ns);
        for Istep=1:Ns
           tline=fgetl(fid);
           if isempty(tline),
               keyboard
           end
           Ispace=1+findstr(tline(2:end),' ');
           aaa{1}=tline(1:(Ispace-1));
           aaa{2}=tline(Ispace:end);
           for J=1:2
              Ipt=findstr(aaa{J},'.'); 
              Icol=min(findstr(aaa{J},':')); 
              Ie=findstr(aaa{J},'E');
              if isempty(Icol)
                  aa(J)=str2num(aaa{J});
              else
                 if Icol<Ipt
                     aa(J)=NaN;
                 else
                     try
                         aa(J)=str2num([aaa{J}(1:(Icol-1)) aaa{J}((Icol+1):end)]);
                     catch
                         aa(J)=NaN;
                     end
                 end
              end
           end
            r(Istep)=aa(1);
            z(Istep)=aa(2);
%             if aa(2)<2
%                 aaa
%             end
           
        end
        %ray = fscanf( fid, '%f', [2 arr.nsteps(ibeam)] );
        %r = ray( 1, : );
        %z = ray( 2, : );
        
        lincol = 'kbgrcmy';
        ii = arr.NumBotBnc(ibeam);
        ii = 1+min([ii length(lincol)-1]);
        plot( r/1000, z, [ lincol(ii) ]);
        %if arr.NumTopBnc > 1 & arr.NumBotBnc > 1
        %  plot( r, z, 'k--' )	% hits both boundaries
        %elseif arr.NumBotBnc > 1
        %  plot( r, z, 'r-' )	% hits bottom only
        %elseif arr.NumTopBnc > 1
        %  plot( r, z, 'b-' )	% hits surface only
        %elseif arr.NumTopBnc == 0 & NBotBnc == 0
        %  plot( r, z, 'b-' )
        %else
        %  plot( r, z, 'y-' )
        %end
        % update axis limits
        rmin = min( [ r/1000 rmin ] );
        rmax = max( [ r/1000 rmax ] );
        zmin = min( [ z zmin ] );
        zmax = max( [ z zmax ] );
        axis( [ rmin, rmax, zmin, zmax ] )
        if rem(ibeam,fix(NBEAMS/10))==0,
            drawnow
        end
        %keyboard
        %end
    end
    
end	% next beam
title(sprintf('Freq: %6.2f, source depth: %6.2f bottom depth: %6.2f',FREQ,z(1),DEPTHB));
grid on
fclose( fid );

hold off
zoom on
