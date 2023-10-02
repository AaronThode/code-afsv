function B = MV_beamforming(~, Ks, angles, freq, Lz, c)
B=zeros(length(freq),length(angles));
if size(Lz,2)>1
    Lz=Lz.';
end
for If=1:length(freq)
    for Iang=1:length(angles)
        % lambda=1500/freq(If);
        w=exp((1i*2*pi*Lz*freq(If)/c)*sin(angles(Iang)*pi/180));

        %B(If,Iang)=real(w'*Ks{If}*w)/(norm(Ks{If})*norm(w).^2);
        w=w/norm(w);
        B(If,Iang)=1./real(w'*(squeeze(Ks(:,:,If))\w));

    end


end

%plot(angles,10*log10(B(If,:)));
end
