%%%%%%%%%%find_similar_elements_ovlap.m%%%%%%%%
%function [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match,Iy_match_diff,Iy_nomatch_diff]=find_similar_elements_ovlap(units,x,dx,y,dy,ovlap_tol,xfeature,yfeature);
%   x, y:  two vectors with times in either seconds or datenumber, generally x is considered the "automated" and y the "manual," or
%           "true"
%   dx,dy: uncertainty in x and y (e.g. duration of signal if x and y are times).  Same units as x an dy
%   overlap_tol: what fraction overlap x(I) and y(J) need to have to be considered a "match".
%   xfeature,yfeature: Optional features of x and y matricies, used to select
%         best x if multiple instances of x within tol of y. [2 x Nsamples].
%         Currently assume xfeature represents frequency; best x match is the one that has largest area.
%   If xfeature not available, best match x is one with longest duration (dx).
%   units:  'datenum' or 'sec'
%   tol:  time tolerance for automated result to be considered a candidate to a given manual detection, in seconds.  Used to speed up
%               processing:  if start time of a given x is within tol of y, and end time of x is withing tol of y, conduct further processing
% Output:
%  Ix_nomatch:  Indicies of x that have no counterpart in y
%  Ix_match:  Indicies of x that match y
%  Iy_match:  Indicies of y that match x, thus
%       x(Ix_match(k))=y(Iy_match(k));
%  Iy_match:  Indicies of y that have no counterpart in x
%  Iy_nomatch_diff:  how far this element of x was to nearest y element
%  Ix_redundant=other elements of x that match this element of y.
%
%  NOTE: the sum of lengths of Ix_match and Ix_nomatch might not equal
%  length of x, if more than one x element lies within tolerance of a given
%  y.

%  NOTE: be sure to compare x and y over same time periods.  For example,
%  if 'y' only covers midnight to noon, make sure the x vector does not
%  cover midnight to 24 hours...

function [Ix_nomatch,Iy_nomatch,Ix_match,Iy_match,Iy_nomatch_diff,Iy_match_diff,Ix_redundant]=find_similar_elements_ovlap(x,dx,y,dy,ovlap_tol,xfeature,yfeature,units,tol)
if ~exist('tol','var')
    tol=60; %Time tolerance window to pick candidates from--used to speed algorithm up
end

Ix_match=[];
Ix_nomatch=[];
Iy_match=[];
Iy_nomatch=[];
Iy_nomatch_diff=[];
Iy_match_diff=[];
Ix_redundant=[];
Nredundant=0;

Nx=[];
if exist('xfeature','var')&&~isempty(xfeature)
    Nfeatures=size(xfeature,1);
    feature_flag=1;
else
    feature_flag=0;
end

if isempty(x)
    Iy_nomatch=1:length(y);
    return
end

if isempty(y)
    Ix_nomatch=1:length(x);
    return;
end

%Make all vectors horizontal
if size(dx,1)>1
    dx=dx';
end

if size(dy,1)>1
    dy=dy';
end

if size(x,1)>1
    x=x';
end

if size(y,1)>1
    y=y';
end

%Convert time units if needed; assume default is second
if ~exist('units','var')||isempty(units)
    units='sec';
end

if ~isempty(findstr(units,'sec'))
    %Do nothing
elseif ~isempty(findstr(units,'date'))
    %%Refer all seconds to the same time origin
    tabs_min=min([x y]);
    x=datenum2sec(x,tabs_min);
    y=datenum2sec(y,tabs_min);
    dx=datenum2sec(dx,0);
    dy=datenum2sec(dy,0);
    
end


%Now all operations are in terms of seconds
Nx=length(x);
Ny=length(y);

%March through reference times
for I=1:Ny  %these are the "true" manual values, for example
    Ix_redundant{I}=[];
    
    fail_flag=0;
    %Before committing resources, check what values match criterial
    Ixp=find((abs(y(I)-x)<=2*tol)&(abs(y(I)+dy(I)-x-dx)<=2*tol));  %Automated signals that that lie within time tolerance
    if isempty(Ixp)
        fail_flag=1;
        ovlap=-1;
        Imin=[];
    else  %potential matches exist
        ovlap=-1*ones(size(Ixp));
        if feature_flag==1
            feature_ovlap=-1*ones(size(Ixp));
        end
        
        for J=1:length(Ixp)
            index=Ixp(J);
            tovlap=min([y(I)+dy(I) x(index)+dx(index)])-max([y(I) x(index)]);  %Time of overlap in sec
            ovlap(J)=(tovlap./min([dx(index) dy(I)]));  %Negative if no overlap
            if feature_flag==1 %If extra features provided
                fovlap=min([yfeature(2,I) xfeature(2,index)])-max([yfeature(1,I) xfeature(1,index)]);
                feature_ovlap(J)=(fovlap./min([diff(yfeature(:,I)) diff(xfeature(:,index))]));  %Negative if no overlap
                
            end
        end
        
        if feature_flag==1
            Imin=find((ovlap>=ovlap_tol)&(feature_ovlap>=ovlap_tol));
        else
            Imin=find(ovlap>=ovlap_tol);
        end
    end
    
    if isempty(Imin)||fail_flag==1  %No automated matches
        Iy_nomatch=[Iy_nomatch I];
        Iy_nomatch_diff=[Iy_nomatch_diff max(ovlap)];
        continue;
    end
    
    Iy_match=[Iy_match I];
    
    
    %%Are there multiple x matches for this y?
    if length(Imin)==1
        Ix_match=[Ix_match (Ixp(Imin))];
        Iy_match_diff=[Iy_match_diff ovlap(Imin)];
        
    else  %multiple x candidates, choose one with largest box area
        %%%Place other x candidates that qualify into Ix_redundant
        if feature_flag==1
            areatb=diff(xfeature(:,Ixp(Imin))).*(dx(Ixp(Imin)));  %dx is time, diff(xfeature) is bandwidth
            [junk,Ibest]=max(areatb);
        else  %Select best match based on longest time duration, best overlap may be a tiny fragment...
            [junk,Ibest]=max(dx(Ixp(Imin)));
        end
        Ix_redundant{I}=Ixp(Imin(setdiff(1:length(Imin),Ibest)));
        Ix_match=[Ix_match Ixp(Imin(Ibest))];
        Iy_match_diff=[Iy_match_diff ovlap(Imin(Ibest))];
        
        %Should we remove other candidates from "false" category, if we have feature information?
        Nredundant=Nredundant+length(Imin)-1;
    end
    
end

%%%Fixed bug Apr 4, 2010, Ix_match may not be monotonic in time...
[Ix_match,Isort]=sort(Ix_match);
Iy_match=Iy_match(Isort);
Iy_match_diff=Iy_match_diff(Isort);
if isempty(Ix_match)
    Ix_nomatch=1:length(x);
    return;
end

%Now look for false alarms in automated data...
Ix_nomatch=setdiff(1:length(x),Ix_match);  %These include true false matches and redundant matches

if Ny==1  %only one 'true' element to compre to.
    return
end
%Sanity check for false detections-do some actually qualify?
Nfunny=length(Ix_nomatch);

Ix_nomatch_ovlap=zeros(1,Nfunny);
for I=1:Nfunny  %For each false detection
    indexx=Ix_nomatch(I);
    Iyp=find((abs(y-x(indexx))<=2*tol)&abs(y+dy-x(indexx)-dx(indexx))<=2*tol);
    if isempty(Iyp)
        Ix_nomatch_ovlap(I)=-1;
    else
        
        for J=1:length(Iyp)
            indexy=Iyp(J);
            tovlap=min([y(indexy)+dy(indexy) x(indexx)+dx(indexx)])-max([y(indexy) x(indexx)]);
            Ix_nomatch_ovlap(I)=max([Ix_nomatch_ovlap(I) (tovlap./min([dx(indexx) dy(indexy)]))]);  %Negative if no overlap
            if Ix_nomatch_ovlap(I)>ovlap_tol
                if feature_flag==1  %If extra features provided
                    fovlap=min([yfeature(2,indexy) xfeature(2,indexx)])-max([yfeature(1,indexy) xfeature(1,indexx)]);
                    feature_ovlap=(fovlap./min([diff(yfeature(:,indexy)) diff(xfeature(:,indexx))]));  %Negative if no overlap
                    if feature_ovlap<ovlap_tol
                        Ix_nomatch_ovlap(I)=-1;
                    end
                    
                end
                
            end
            
        end
    end
end
Itruly_bad=find(Ix_nomatch_ovlap<ovlap_tol);  %These are indeed false alarms
Ix_nomatch=Ix_nomatch(Itruly_bad);  %Redundant matches removed...
Ix_nomatch_ovlap=Ix_nomatch_ovlap(Itruly_bad);
%fprintf('Out of %i original false matches, %i were redundant and %i were actually overlapped more than %6.2f  with a true value and not included in final false result \n', Nfunny,Nredundant, ...
%    Nfunny-length(Itruly_bad),ovlap_tol);


%disp(sprintf('Detected calls: %i  plus missed calls: %i = %i should be %i',length(Iy_match),length(Iy_nomatch),length(Iy_match)+length(Iy_nomatch),length(y)));
disp('');
% Iskip=find(diff(Iy_match)>1);
% for I=1:length(Iskip),
%     Istart=Iy_match(Iskip(I))+1;
%     Iend=Iy_match(Iskip(I)+1)-1;
%     Iy_nomatch=[Iy_nomatch Istart:Iend];
% end

end

function x=datenum2sec(x,tabs0)
tmp=datevec(x-tabs0);
x=tmp(:,6)+60*tmp(:,5)+3600*tmp(:,4)+24*3600*tmp(:,3);
end

