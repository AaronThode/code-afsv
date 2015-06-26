function varargout = plotray( rayfil )

% Plot the RAYfil produced by Bellhop
% usage: plotray( rayfil )
% where rayfil is the ray file (extension is optional)
% e.g. plotray( 'foofoo' )
%
% MBP July 1999

global units jkpsflag

if ( strcmp( rayfil, 'RAYFIL' ) == 0 && isempty( findstr( rayfil, '.ray' ) ) )
    rayfil = [ rayfil '.ray' ]; % append extension
end

% plots a BELLHOP ray file

fid = fopen( rayfil, 'r' );   % open the file
if ( fid == -1 )
    disp( rayfil );
    errordlg( 'No ray file exists; you must run BELLHOP first (with ray ouput selected)', 'Error' );
end

% read header stuff

TITLE  = fgetl(  fid );
FREQ   = fscanf( fid, '%f', 1 );
NBEAMS = fscanf( fid, '%i', 1 );
DEPTHT = fscanf( fid, '%f', 1 );
DEPTHB = fscanf( fid, '%f', 1 );

% Extract letters between the quotes
nchars = strfind( TITLE, '''' );   % find quotes
TITLE  = [ TITLE( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 7 - ( nchars( 2 ) - nchars( 1 ) ) ) ];
TITLE  = deblank( TITLE );  % remove whitespace

% read rays

% set up axis lengths for publication
if ( jkpsflag )
    set( gcf, 'Units', 'centimeters' )
    set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
    
    set( gca, 'Position', [ 2 2 14.0  7.0 ] )
    %set( gcf, 'PaperPosition', [ 3 3 19.0 10.0 ] )
end

set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down

xlabel( 'Range (m)' )
if ( strcmp( units, 'km' ) )
    xlabel( 'Range (km)' )
end
ylabel( 'Depth (m)' )
title( TITLE )
hold on

% axis limits
rmin = +1e9;
rmax = -1e9;
zmin = +1e9;
zmax = -1e9;

for ibeam = 1:NBEAMS
    alpha0    = fscanf( fid, '%f', 1 );
    nsteps    = fscanf( fid, '%i', 1 );
    NumTopBnc = fscanf( fid, '%i', 1 );
    NumBotBnc = fscanf( fid, '%i', 1 );
    if isempty( nsteps ); break; end
    ray = fscanf( fid, '%f', [2 nsteps] );
    r = ray( 1, : );
    z = ray( 2, : );
    
    if ( strcmp( units, 'km' ) )
        r = r / 1000;   % convert to km
    end
    
    %lincol = 'kbgrcmy';
    %ii = NumBotBnc;
    %ii = mod( ii, 3 ) + 1;
    %plot( r, z, lincol(ii) );
    if NumTopBnc > 1 && NumBotBnc > 1
        plot( r, z, 'k' )    % hits both boundaries
    elseif NumBotBnc >= 1
        plot( r, z, 'b'  )	% hits bottom only
    elseif NumTopBnc >= 1
        plot( r, z, 'k'  )	% hits surface only
    else
        plot( r, z, 'r' )
    end
    % update axis limits
    rmin = min( [ r rmin ] );
    rmax = max( [ r rmax ] );
    zmin = min( [ z zmin ] );
    zmax = max( [ z zmax ] );
    if ( zmin == zmax ) % horizontal ray causes axis scaling problem
        zmax = zmin + 1;
    end
    axis( [ rmin, rmax, zmin, zmax ] )
    
    if rem( ibeam, fix( NBEAMS / 10 ) ) == 0,    % flush graphics buffer every 10th ray
        drawnow
    end;
    %end
end	% next beam

fclose( fid );

hold off
zoom on

if ( nargout == 1 )
    varargout{ 1 } = findobj( 'Type', 'Line' );   % return a handle to the lines in the figure
end

% fixed size for publications
if ( jkpsflag )
    set( gca, 'Units', 'centimeters' )
    set( gca, 'Position', [ 2 2 14.0  7.0 ] )
    set(gcf, 'PaperPositionMode', 'auto');

    %set( gcf, 'Units', 'centimeters' )
    %set( gcf, 'PaperPosition', [ 3 3 19.0 10.0 ] )
    set( gcf, 'Units', 'centimeters' )
    set( gcf, 'Position', [ 3 15 19.0 10.0 ] )
end