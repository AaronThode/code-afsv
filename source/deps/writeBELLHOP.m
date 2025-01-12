
%%%%%%%%%%%%%%%%%%%writeBELLHOP.m%%%%%%%%%%%%%%%%%%%%%
% function writeBELLHOP(sd,rd,r,envname,freq, sspopt, lthick, cmedtop, c,run_chc,angle_vec,cmedbot,rho,pwavatten,Darray,zbox,stepsize)
% Aaron Thode
% September, 2004 
% 
% Added option to write detailed ssp and rhofor two layers
% Added rigid bottom feature Sept 18, 1997
% Jan 2006: Bottom bathymetry added..
%Define list of input variables
% Can handle constant and linear sound-speed profiles
% sd array of source depths in meters
% rd array of receiver depths in meters
% r array of ranges in km
% envname -name of environment file desired.
% lthick-thickness of each layer. First element is water depth
% cmedtop: if greater than 10000 use a perfectly reflecting bottom
%   profile...
%  ssp- Specifies a more complex sound speed profile and density within a layer.
%	Option#1: A two column matrix with depth in col one, sound speed in col two
%	Option#2: A cell array with each cell representing a layer,
%		  col one is depth in m, sound speed col two, 
%			shear col 3, optional density (g/cm^3) col 4
%  run_chc-whether one uses "R", "E", "A", "a' "I" (incoherent) or "C"(coherent),"S" (Lloyd's Mirror)
% angle_vec: vector of angles: if [amin; amax; number_of_angles] --either row or column
%       vector..
% Darray: two-column bathymetry matrix of ranges in km and depths in
%       meters, first row is depth under source.  Needs to be at least two rows long...
% zbox: Specifies max depth to trace rays--useful for killing bottom
%           reflections.
% sspopt: three chars descripting ssp options e.g. 'SVN'

function writeBELLHOP(sd,rd,r,envname,freq, sspopt, lthick, cmedtop, c,run_chc,angle_vec,cmedbot,rho,pwavatten,Darray,zbox,stepsize)
%Define list of variables
Nmesh=1000;  %Dummy variables for Bellhop
Nmedia=length(Nmesh);  %Number of media layers.
Nmesh=round(Nmesh);
mediumdepth=[0 cumsum(lthick)];
options1=sspopt ;  % Added by ADS 2006.04.23
%options1='SVN' ;  %splice interpolation

if exist('Darray','var')&&~isempty(Darray)
    if size(Darray,2)>2,Darray=Darray';end
    fid=fopen([envname '.bty'],'w');
    fprintf(fid,'L \n');
    fprintf(fid,'%i \n',size(Darray,1));
    fprintf(fid,'%6.2f %6.2f\n',Darray');
    fclose(fid);
end

%Check for shear
if size(cmedtop,1)==1&size(cmedtop,2)>1
	cmedtop(2,:)=zeros(1,Nmedia+1);
	cmedbot(2,:)=zeros(1,Nmedia+1);
	pwavatten(2,:)=zeros(1,Nmedia+1);
    %elseif size(cmedbot,1)==1,
%	error('cmedbot or pwavatten needs to include shear too!');
end

rmax=max(r);
rmsrns=zeros(size(Nmesh));
%Write the environmental file!
fid=fopen([envname '.env'],'w');
fprintf(fid,'''%s'' \n', envname);
fprintf(fid,'%10.6f \n',freq);
fprintf(fid,'%i \n',Nmedia);
fprintf(fid,'''%s'' \n',options1);
for I=1:Nmedia
    fprintf(fid,'%i %i %10.4f \n',Nmesh(I),rmsrns(I),mediumdepth(I+1));
    %This section changed by Aaron Jan 2001 to allow multilayer profiles
    if (I==1 && iscell(c)==0&&size(c,1)>1) %I just have a complicated water profile
        fprintf(fid,'%16.12f %16.12f 0 1 0 0 \n',c');
    elseif (iscell(c)==1&&I<=length(c)),
        disp(sprintf('Writing profile to layer %i',I));
        co=c{I};
        if size(co,2)==3, %no density information, assume values of one
            co=[co ones(size(co,1),1) pwavatten(1,I)*ones(size(co,1),1)];
        elseif size(co,2)==4,
            co=[co pwavatten(1,I)*ones(size(co,1),1)];
        end
        fprintf(fid,'%16.12f %16.12f %16.12f %6.4f %6.4f 0\n',co');
    else
        fprintf(fid,'    %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f \n', ...
            mediumdepth(I),cmedtop(1,I),cmedtop(2,I), ...
            rho(I),pwavatten(1,I),pwavatten(2,I));
        fprintf(fid,'    %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f \n', ...
            mediumdepth(I+1),cmedbot(1,I),cmedbot(2,I), ...
            rho(I),pwavatten(1,I),pwavatten(2,I));
    end
end
if cmedtop(1,end)<10000 %i.e bottom fast
    if ~exist('Darray','var')||isempty(Darray)||all(diff(Darray(:,2))==0)  %If a flat bathymetry
        fprintf(fid,'''A'' , 0 \n');
        fprintf(fid,'    %10.4f %10.4f %10.4f %10.6f %10.6f 0 \n', ...
                mediumdepth(Nmedia+1),cmedtop(1,Nmedia+1), ...
                0,rho(Nmedia+1),pwavatten(1,Nmedia+1));
    else
        fprintf(fid,'''A*'' , 0 \n');
        fprintf(fid,'    %10.4f %10.4f %10.4f %10.6f %10.6f 0 \n', ...
                mediumdepth(Nmedia+1),cmedtop(1,Nmedia+1), ...
                0,rho(Nmedia+1),pwavatten(1,Nmedia+1));
    end
        
elseif ~exist('Darray','var')||isempty(Darray)
    fprintf(fid,'''R'',0 \n');
else 
    fprintf(fid,'''RB'',0 \n');
end

if length(sd)==1
    fprintf(fid,'%i \n %6.2f \n',1,sd);
else
    nsd=length(sd);
    %fprintf(fid,'%i %6.2f %6.2f / \n',nsd,sd(1),sd(nsd));
    fprintf(fid,'%i \n ',nsd);
    fprintf(fid,'%6.2f ',sd);
    fprintf(fid,'\n');
    
end

if length(rd)==1
    fprintf(fid,'%i \n %6.2f \n',1,rd);
elseif rd(1)<0, %Assume evenly spaced
    disp('evenly spaced');
    nrd=length(rd);
    fprintf(fid,'%i \n %6.2f %6.2f / \n',nrd,abs(rd(1)),abs(rd(nrd)));
else
    nrd=length(rd);
    fprintf(fid,'%i \n',nrd);
    fprintf(fid,'%6.2f ',rd);
    fprintf(fid,'\n');
end
if length(r)==1
    fprintf(fid,'%i \n %6.3f \n',1,r);
elseif r(1)<0 %Assume evenly spaced
    disp('evenly spaced');
    nr=length(r);
    fprintf(fid,'%i \n %6.3f %6.3f / \n',nr,abs(r(1)),abs(r(nr)));
else
    nr=length(r);
    fprintf(fid,'%i \n',nr);
    fprintf(fid,'%6.3f ',r);
    fprintf(fid,'\n');
end

fprintf(fid,'''%s''\n',run_chc);

if exist('angle_vec','var')&&~isempty(angle_vec)
    if length(angle_vec)==3  %Compute a range
        fprintf(fid,'%i \n %6.2f %6.2f / \n',angle_vec(3),min(angle_vec(1:2)),max(angle_vec(1:2)));
    elseif length(angle_vec)>3
        fprintf(fid,'%i \n',length(angle_vec));
        fprintf(fid,'%8.2f   ',angle_vec);
        fprintf(fid,'\n');
        
        % fprintf(fid,'%i \n %6.2f / \n',length(angle_vec),angle_vec);
    end
    
elseif ~exist('angle_vec','var')||isempty(angle_vec)
    fprintf(fid,'0 \n -89 89 /\n');  %For coherent computation let program compute automatically
else
    fprintf(fid,'%i \n',length(angle_vec));
    fprintf(fid,'%6.1f',angle_vec);
    fprintf(fid,'\n');
    %fprintf(fid,'%i \n %6.2f %6.2f / \n',angle_vec(3),min(angle_vec(1:2)),max(angle_vec(1:2)));
end

%%Now box box of output results
%%Numerical step size should be on order of sound speed profile step.  If
%%dense profile, explicitly request a step size from a hopefully interpolated sound speed...
if exist('zbox','var')&~isempty(zbox)
    Dmax=zbox;
elseif exist('Darray','var')&~isempty(Darray)
    Dmax=1.01*max(Darray(:,2));
else
    Dmax=1.01*lthick(1);
end
if ~exist('stepsize','var')
    try
        stepsize=c(2,1)-c(1,1);
    catch
        stepsize=1;
    end
end
if stepsize>15
    fprintf(fid,'0 %6.2f %6.2f',Dmax,1.01*max(rmax));
else
    fprintf('Fine-grained ssp: setting ray step size to %8.4f\n',0.5*stepsize);
    fprintf(fid,'%6.2f %6.2f %6.2f',0.5*stepsize,Dmax,1.01*max(rmax));
    
end
fprintf(fid,'\n');
fclose(fid);

