function [B,dt,t,wout]=conventional_beam_timedelay(Ks,angles,freq,Lz,c,Nfft,Fs,yesnorm)
%function [B,dt,t,wout]=conventional_beam_timedelay(Ks,angles,freq,Lz,c,Nfft,Fs,yesnorm)

% angles is two-element vector [look1 look2];
%   plotting (t,B) shows the time delay between beam (angle(1))=beam 1 and beam
%   (angle(2)) =beam2.  A positive value means that on element
%   one (first Ks column, first x column) beam 2 signal leads (arrives
%   before) beam 1.  A negative value means that beam 2 arrives after beam
%   1 on element 1.  So x-axis of plot(t,B) can be interprets as how beam 1
%   'lags' beam 2.
%
%   dt: time delay in msec. The sign of dt corresponds with convention
%   above

B=zeros(Nfft,1);
FF=linspace(0,Fs,Nfft);
[junk,Iff]=min(abs(FF-freq(1)));
Iff=Iff+(0:(length(freq)-1));


if size(Lz,2)>1
    Lz=Lz.';
end


winn=hanning(length(Lz));

B=zeros(Nfft,length(angles));
for If=1:length(freq)
    K=squeeze(Ks(:,:,If));
    if exist('yesnorm', 'var')
        K=K/norm(K);
    end
    % lambda=1500/freq(If);
    w1=exp((-1i*2*pi*Lz*freq(If)/c)*sin(angles(1)*pi/180));
    w1=w1.*winn;
    w1=w1/norm(w1);
    
    w1K=w1'*K;
    for Ia=1:length(angles)
        w2=exp((-1i*2*pi*Lz*freq(If)/c)*sin(angles(Ia)*pi/180));
        w2=w2.*winn;
        w2=w2/norm(w2);
        
        B(Iff(If),Ia)=(w1K*w2);
        
    end
end

t=((-Nfft/2):(Nfft/2-1))/Fs;
figure(get(gcf,'Number')+1)

for Ia=1:length(angles)
    B(:,Ia)=fftshift(real(ifft(B(:,Ia),Nfft)));
    subplot(length(angles),1,Ia);
    plot(t,B(:,Ia));grid on
    [junk,Ibest]=max(abs(B(:,Ia)));
    titlestr=sprintf(' Beam 1 (%6.2f deg) lags Beam %i: (%6.2f deg) by: %6.2f msec',angles(1),Ia, angles(Ia), t(Ibest)*1000);
    disp(titlestr);
    title(titlestr);
    xlabel('lag (msec)');
    dt(Ia)=t(Ibest)*1000;
end


end

%plot(angles,10*log10(B(If,:)));

