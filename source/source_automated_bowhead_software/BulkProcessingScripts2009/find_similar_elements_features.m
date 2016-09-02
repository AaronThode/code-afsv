%%%%%%%%%%find_similar_elements.m%%%%%%%%
%function [Iy_match,y_match,xy_diff]=find_similar_elements(x,y,tol,x_feature,y_feature,feature_tol,Idebug);

% x, y:  two vectors with times in same units
% tol: how close x(I) and y(J) have to be to be considered a "match".
% maxtol:  Sometimes, x or y start or end midway through the other vector.
% In that case, if x-max(y) or min(y)-x exceeds maxtol, remove those
%   x-points from consideration
% Output:
%  Ix_nomatch:  Indicies of x that have no counterpart in y
%  Ix_match:  Indicies of x that match y
%  Iy_match:  Indicies of y that match x, thus
%       x(Ix_match(k))=y(Iy_match(k));
%  Iy_match:  Indicies of y that have no counterpart in x
%  Ix_nomatch_diff:  how far this element of x was to nearest y element
%  Nx: number of x points that survive maxtol filter...

function [Iy_match,y_match,xy_diff]=find_similar_elements(x,y,tol,x_feature,y_feature,feature_tol,Idebug);
 


if isempty(x)||isempty(y)
    return;
end
Iy_match=zeros(1,length(x));
 xy_diff=Iy_match;
y_match=Iy_match;
Nx=length(x);
Ny=length(y);

for I=1:Nx,

    Iposs=find(abs(y-x(I))<=tol);

    
    %Check each possibility
    [compare,Ibest]=min(abs(x_feature(I)-y_feature(Iposs)));
    if compare<=feature_tol,
        Iy_match(I)=Iposs(Ibest);
        xy_diff(I)=y(Iy_match(I))-x(I);
        y_match(I)=y(Iy_match(I));
    else
        Iy_match(I)=-1;
        xy_diff(I)=NaN;
        y_match(I)=-1;
    end
    
    
    if Idebug.yes>0,
        tlen=4;
        figure(3);
        subplot(1+length(Iposs),1,1);
        file_dir=sprintf('%s/Site_0%i/%s',Idebug.rawdatadir,Idebug.Isite,Idebug.xdirname);
        fname=dir(sprintf('%s/*%s*.sio',file_dir,Idebug.current_date ));

        [xx,t,head]=readsiof([file_dir '/' fname(1).name],x(I)-tlen/2,tlen);

        spectrogram(xx,128,96,128,head.Fs,'yaxis');
        title(sprintf('To match: %ith element at %s',I,Idebug.xdirname));
        for II=1:length(Iposs),
            subplot(1+length(Iposs),1,II+1);
            file_dir=sprintf('%s/Site_0%i/%s',Idebug.rawdatadir,Idebug.Isite,Idebug.dirname);
            fname=dir(sprintf('%s/*%s*.sio',file_dir,Idebug.current_date ));

            [xx,t,head]=readsiof([file_dir '/' fname(1).name],y(Iposs(II))-tlen/2,tlen);

            spectrogram(xx,128,96,128,head.Fs,'yaxis');
            infostr=sprintf('%ith element at %s, Feature value %6.2f, wanted %6.2f, time diff %6.2f s', ...
                    Iposs(II),Idebug.dirname,y_feature(Iposs(II)),x_feature(I),y(Iposs(II))-x(I));
            if II==Ibest,
                title([infostr ', selected']);
            else
                title(infostr);
            end


        end
       keyboard;
    end

end



