function [vx,vy,T,F,P]=demultiplex_DIFAR(x,Fs,Nfft,ovlap)


%Nfft=2^15;
%ovlap=0.75;
dn=(1-ovlap)*Nfft;

[S,F,T,P] = spectrogram(double(x),Nfft,round(ovlap*Nfft),Nfft,Fs);
P=10*log10(abs(P));

% imagesc(T,F,P);
% grid on;axis('xy');ylim([0 500]);
% %
% Bmean=median(P')';

%%%%Locate pilot signals
Bwidth=50;
[I75]=find(F>=7500-0.5*Bwidth&F<=7500+0.5*Bwidth);
[~,Ipilot]=max(P(I75,:));
Ipilot=mode(Ipilot);
Ipilot75=I75(Ipilot);

[I15]=find(F>=2*7500-0.5*Bwidth&F<=2*7500+0.5*Bwidth);
[~,Ipilot]=max(P(I15,:));
Ipilot=mode(Ipilot);
Ipilot15=I15(Ipilot);
pilot_phase=(exp(-1i*angle(S(Ipilot15,:))));

%%%Check relative phases of signals
%rel_phase1=angle(S(Ipilot15,:));
%rel_phase2=2*angle(S(Ipilot75,:));
%plot(T,unwrap(rel_phase1-rel_phase2))

P0=S(2:(Ipilot75-1),:);
P=P(2:(Ipilot75-1),:);
N=Ipilot75-2;
%Sm=flipud(S((Ipilot75+1):(Ipilot15-1),:));
Sm=flipud(S((Ipilot75+1)+(0:N-1),:));
Sp=S((Ipilot15+1)+(0:N-1),:);

[Nf,Nt]=size(P0);

F=F(1:Nf);
%%%%Multiplex by phase of 15 kHz pilot

Sm=Sm.*(ones(Nf,1)*pilot_phase);
Sp=Sp.*(ones(Nf,1)*pilot_phase);

%%%%%Build Mi and Mq

ReMi=real(Sp+Sm);
ReMq=imag(Sp+Sm);

ImMi=imag(Sp-Sm);
ImMq=-real(Sp-Sm);

Mi=ReMi+1i*ImMi;  %%NS
Mq=ReMq+1i*ImMq;  %EQ

%vx=squeeze(real((P.*conj(Mi))));
%vy=squeeze(real((P.*conj(Mq))));

vx=-squeeze(real((P0.*conj(Mi))));
vy=squeeze(real((P0.*conj(Mq))));

