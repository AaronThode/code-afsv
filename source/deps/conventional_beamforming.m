%function [B,wout]=conventional_beamforming(Ks,angles,freq,Lz,c,yesnorm)
function [B,wout]=conventional_beamforming(Ks,angles,freq,Lz,c,yesnorm)
    %Inputs:
    %  Ks:  cross-spectral density matrix [Nel, Nel, Nfreq];
    %  freq: frequencies in Hz
    %  Lz: positions in meters
    %  c: sound speed m/sec
% output: wout(element,freq,angle));

B=zeros(length(freq),length(angles));
if size(Lz,2)>1
    Lz=Lz.';
end

if nargout==2
    wout=zeros(length(Lz),length(freq),length(angles));
end
winn=hanning(length(Lz));
for If=1:length(freq)
    for Iang=1:length(angles)
        % lambda=1500/freq(If);
        w=exp((-1i*2*pi*Lz*freq(If)/c)*sind(angles(Iang)));
        w=w.*winn;
        
        w=w/norm(w);
        K=squeeze(Ks(:,:,If));
        if exist('yesnorm', 'var')
            K=K/norm(K);
        end
        B(If,Iang)=real(w'*K*w);
        if nargout==2
            wout(:,If,Iang)=w;
        end
    end
    
    
end

%plot(angles,10*log10(B(If,:)));
end
