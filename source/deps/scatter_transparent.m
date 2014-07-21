%function scatter_transparent(x,y,sizze,collor,alpha)
% x: row or column vector of x plotting data
% y: row or column vector of y plotting data
% sizze: size of marker in terms of x and y units
% collor: color desired.  e.g. 'k' is black, 'b' is blue; see 'plot' for
%   choices
% alpha: transparency.  0 is invisible, 1 is opaque
%   Quick test:
%   x=[5 10];y=x;close all;scatter_transparent(x,y,.3,'k',0.05)

function scatter_transparent(x,y,sizze,collor,alpha)

Np=length(x);
if Np~=length(y)
    error('x and y need to be the same length');
    return
end

if size(x,1)>1
    x=x';
end

if size(y,1)>1
    y=y';
end

xplot=ones(4,1)*x+[zeros(1,Np);  -0.5*sizze*ones(1,Np); zeros(1,Np);0.5*sizze*ones(1,Np)];
yplot=ones(4,1)*y+[-0.5*sizze*ones(1,Np); zeros(1,Np);  0.5*sizze*ones(1,Np);zeros(1,Np) ];

p=patch(xplot,yplot,collor,'marker','o','facecolor',collor,'FaceAlpha',alpha);
set(p,'edgecolor','none','FaceAlpha',alpha)
