%%%contour_2D_histogram.m%%
%function [xgrid,ygrid,count]=contour_2D_histogram(x,y,Nx,Ny),

function [xgrid,ygrid,count]=contour_2D_histogram(x,y,Nx,Ny)

if length(Nx)==1,
    xgrid=linspace(min(x),max(x),Nx);
    ygrid=linspace(min(y),max(y),Ny);
else
    xgrid=Nx;
    ygrid=Ny;
end


Nx=length(xgrid);
Ny=length(ygrid);
count=0.01*ones(Ny,Nx);

for I=1:(Nx-1),
    Isub=find(x>=xgrid(I)&x<=xgrid(I+1));
    Ny=histc(y(Isub),ygrid);
    count(:,I)=Ny;
end

%imagesc(xgrid,ygrid,log10(count));colorbar;
%set(gca,'fontweight','bold','fontsize',14);grid on;

%contour(xgrid,ygrid,count,1:2:12);
%keyboard;
