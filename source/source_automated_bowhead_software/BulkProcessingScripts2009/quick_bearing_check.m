%%%%%quick_bearing_check.m

run_options.bearing_alg='sel_ratio'
[locations,tet]=compute_bearings(locations,station,goodFile,param,run_options.bearing_alg,0);
tet_time=tet;
run_options.bearing_alg='spectrogram';
[locations,tet]=compute_bearings(locations,station,goodFile,param,run_options.bearing_alg,0);
tet_specgram=tet;

t=1:(7*length(locations));

Igood=find(abs(tet_time.ang(:))>0);
subplot(2,1,1)
 x=tet_time.ang(:);
 y=tet_specgram.ang(:);
 plot(t(Igood),x(Igood),'kx',t(Igood),y(Igood),'go')

subplot(2,1,2)
x=tet_time.sd(:);
y=tet_specgram.sd(:);
plot(t(Igood),x(Igood),'kx',t(Igood),y(Igood),'go')


[Ntime,Xtime]=hist(x(Igood),0:.2:5);
[Ns,Xs]=hist(y(Igood),0:.2:5);
bar(Xs'*[1 1],[Ntime' Ns'],'grouped')
xlabel('Standard bearing error(deg');
ylabel('Samples');
legend('time','spectrogram')