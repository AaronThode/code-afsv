function status = write_covmat(~, freqvec, Cov, depth, srcname, titlestr)
% write_covmat(farray,Ks,depth,srcname,titlestr)
%
% Write Covariance Matrices to file cov_dpss.in
%
% farry = vector of frequencies in Hertz
% Ks     = matrix of concatenated covariance matrices
%           this must be of size [Nel Nel nfreq]
% depth  = hydrophone depths.  First element is shallowest depth
% srcname =file name
% titlestr = title to write to top of input file.
%[n]= size(Cov,1);
%nf=size(Cov,3);
%if(nf ~= length(freqvec))
%   fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
%           size(Cov),length(freqvec));
%  error('sizes of Cov and freqvec are not compatible');
%end

[n,m,F]= size(Cov);
if(m ~= n)
    fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
        size(Cov),length(freqvec));
    error('sizes of Cov and freqvec are not compatible');
end
if(F ~=length(freqvec))
    fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
        size(Cov),length(freqvec));
    error('sizes of 3Cov and freqvec are not compatible');
end

if(n ~= length(depth))
    fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
        size(Cov),length(freqvec));
    error('sizes of Cov and depth are not compatible');
end


%depth=rd;
%disp('WARNING! CONJUGATING COV to be used by SAGA!');
%Cov=conj(Cov);

fidout =  fopen(sprintf('%s.in',srcname),'w');

for ifreq=1:length(freqvec)
    fprintf(fidout,' %s\n',titlestr);
    fprintf(fidout,' %f     0.000 dB\n',freqvec(ifreq));
    fprintf(fidout,' %d\n',n);
    fprintf(fidout,' %8.2f\n',depth);
    %   cols = [ (ifreq-1)*n+1:ifreq*n ];
    %   Cx = Cov(:,cols); % cut one cov-matrix out of the bunch
    Cx=squeeze(Cov(:,:,ifreq));
    for row=1:n
        for col=1:n
            fprintf(fidout,'%10d%10d (%E,%E) \n',row,col, ...
                real(Cx(row,col)),imag(Cx(row,col)));
        end
    end
end
fprintf(fidout,'!  \n');
status=fclose(fidout);
unix('ls -l *.in');
end
