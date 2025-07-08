%function [B]=conventional_beamforming(Ks,angles,freq,Lz,c,yesnorm)
function [B]=conventional_beamforming_fast(Ks,angles,freq,Lz,c,yesnorm)
%Inputs:
%  Ks:  cross-spectral density matrix [Nel, Nel, Ntimes, Nfreq];
%  angles:  angles in degrees
%  freq: frequencies in Hz
%  Lz: positions in meters
%  c: sound speed m/sec
% output: wout(element,freq,angle));


if ndims(Ks)==3
    B=zeros(length(freq),length(angles));
elseif ndims(Ks)==4
    B=zeros(length(freq),length(angles), size(Ks,3));
end

if size(Lz,2)>1
    Lz=Lz.';
end

%winn=hanning(length(Lz));

if ndims(Ks)==3
    for If=1:length(freq)

        w=exp((-1i*2*pi*Lz*freq(If)/c)*sind(angles))./sqrt(length(Lz));
        K=squeeze(Ks(:,:,If));
        if exist('yesnorm', 'var')
            K=K/norm(K);
        end

        Kw=K*w;
        B(If,:)=sum(Kw.*conj(w));
    end

    %%pagemtimes will multiply out whole matrix

elseif ndims(Ks)==4
    for If=1:length(freq)

        w=exp((-1i*2*pi*Lz*freq(If)/c)*sind(angles))./sqrt(length(Lz));
        K=squeeze(Ks(:,:,:,If));
        

        Kw=pagemtimes(K,w);

        
        B(If,:,:)=real(squeeze(sum(Kw.*conj(w))));
    end

end



%plot(angles,10*log10(B(If,:)));
end
