
%%%%%%%%%%%%%WNC.m%%%%%%%%%%%%%%
% March 19, 1996
%function [S,out]=WNC(U,D,d,Cnst,edB0),
% U is eigenvectors of CSDM,
% D is vertical vector of eigenvalues of CSDM
% d is replica vector, edB0 is initial guess for level of noise
% Cnst is the constraint in linear units, must be less than 10*log(size(K,2));
% edB0 should be set to 10*log10(mean(D))-60;
function S=WNC(U,D,d,Cnst,edB0)
global Yabs
% fprintf(fid,'\tKinv=eiginv(Ktot,10e-6); \n');
%         fprintf(fid,'\t[U,SIG,V]=svd(Ktot);        \n');
%         fprintf(fid,'\tSIG=(diag(SIG));SIGmax=mean(SIG) \n');
%         fprintf(fid,'\tedB0=10*log10(real(mean(SIG)))-60; \n');
%         %Cnst=input('Enter a negative dB value for WN constaint:');
%         %Cnst=-5;
%         fprintf(fid,'\tCnstdB=%8.4f;Cnst=10^(CnstdB/10); \n',Cnst);
%

if isnan(d),
    disp('isnan');
    S=NaN;edB=edB0;
    return;
end
Y=U'*d;
Yabs=abs(Y).^2;
out=WNCf(edB0);
if out>0,
    edB=edB0;
else
    edB=fzero('WNCf',edB0,1e-3);
end
e=10^(edB/10);
D0=diag((1./(SIG+e)));
D1=diag(D)*D0.^2;
S=(Y'*D1*Y)/(Y'*D0*Y).^2;

    function out=WNCf(edB00)
        %pause;
        e=10^(edB00/10);
        %e=edB;
        D0=D+e;
        D1=1./D0;
        %GdB=20*log10(Yabs'*D1)-10*log10(Yabs'*(D1.^2));
        GdB=(Yabs'*D1).^2/(Yabs'*(D1.^2));
        out=GdB-Cnst;
        
    end
end


