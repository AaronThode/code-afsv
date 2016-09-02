 function [mt,b,bcsd,sigdb,stndb,dfreq,duration,bspread,mu,kappa,sd]= ...
   mtbrg4(idb,brefd,om,cosch,sinch,ts,DasarID,moreplots)
% function [mt,b,bcsd,sigdb,stndb,dfreq,duration,bspread,mu,kappa,sd,elscat,pewf,pnsf,pzchf]= ... 
%    mtbrg4(idb,brefd,om,cosch,sinch,zch,ts,DasarID,moreplots)
% cloned from mtbrg2 intended for use with dd using Chris Nations bootstrap 
% measures and returns time and bearing of a ping or other acoustic event
% whose absolute bearing is known , and has been displayed in spectrogram form
% on Dasar number idb
% by Matlab script dd.
% Modified to additionally process the DASAR B vertical Channel
% Outputs:
% In 2008 modified to remove vertical channels
% mt    is measured start time as Matlab time vector
% b     is measured bearing in degrees from scatter plot
% bcsd  is bearing measured by cross-spectral density
% sigdb is rms value of signal in selected band during selected interval in db//1 muPa
% stndb is signal to noise ratio in db in selected band
% dfreq(1:2) is lower and upper freq limits of call in Hz
% duration is length of call from start to end in sec
% bspread is angular spread of lobe in scatterplot
% mu is the average bearing by the bootstrap routine
% kappa is the kappa parameter from the bootstrap routine
% sd is the sd (~standard deviation) from ditto
% elstat is elevation angle in degrees as measured from scatterplot
% pewf,pnsf,pzchf are components of scatterplots
%
% inputs:
% refb is the grid bearing of the DASAR cosine axis
% om,cosch,sinch are the time series in volts
% ts is the Matlab time vector of start of screen
% DasarID is Dasar 2-char id.
%
% The operator must click at 3 points:
%  a)at a time at the start of the event. This becomes the event time.
%  b)at the start and lowest freq of the clearest part of the signal
%   (may or may not be the same point as a)
%  c)at the end time and frequency of the clearest part of the signal
% This allows separate timing and bearing measurement.
% mbrg then 
%   a)filters the input signal in the selected band
%     in the time range selected.
%   b)replays the sound
%   c)calculates rms energy in band, correcting for transducer cal
%   d)measures the direction
%   e)repeats a-c above for 1-sec sample starting 2-sec before start.
% mbrg does filtering with Butterworth Band Pass defined in window

%=====================================
rpdg = pi/180;
%
% define window in time and frequency
h=figure(20)
set(h,'Name','Measure 1 Dasar');
Fs=1000;
[wt,e] = ginput(3);
% Operator can change mind and not process by entering <cr>
if length(wt) < 3 
    mt=[0 0 0 0 0 0];
    b=NaN;
    bcsd=0;
    sigdb=0;
    stndb=0;
    dfreq=0;
    duration=0;
    bspread=0;
    mu=0;
    kappa=0;
    sd=0;
    elscat = 0;
    pewf(1:10)=0;
    pnsf(1:10)=0;
    pzchf(1:10)=0;
    return
end
% wt(3) and e(3) must be sorted into increasing order of wt
if wt(3)<wt(2)
    x=wt(3);
    wt(3)=wt(2);
    wt(2)=x;
    x=e(3);
    e(3)=e(2);
    e(2)=x;
end
if wt(2)<wt(1)
    x=wt(2);
    wt(2)=wt(1);
    wt(1)=x;
    x=e(2);
    e(2)=e(1);
    e(1)=x;
end
if wt(3)<wt(2)
    x=wt(3);
    wt(3)=wt(2);
    wt(2)=x;
    x=e(3);
    e(3)=e(2);
    e(2)=x;
end
%wt
%e
%
% If first time is less than 2 sec, noise cannot be measured.
% display error message directing earlier start.
if wt(1)<2
    figure(16)
    clf
    h = axes('Position',[0 0 1 1],'Visible','off');
    text(.05,.9,'First time is too close to start of Window','Fontsize',20)
    text(.05,.8,'Shift window 30 sec earlier and re-measure','Fontsize',20)
    text(.05,.7,'Hit any key to proceed','Fontsize',20)
    beep
    pause
    close(16)
    mt=[0 0 0 0 0 0];
    b=NaN;
    bcsd=0;
    sigdb=0;
    stndb=0;
    dfreq=0;
    duration=0;
    bspread=0;
    mu=0;
    kappa=0;
    sd=0;
    elscat = 0;
    pewf(1:10)=0;
    pnsf(1:10)=0;
    pzchf(1:10)=0;
    return
end
%
% save start and end frequencies
fstart=e(2);
fend=e(3);
% put freqs in ascending order
if e(2) > e(3)
   s=e(2);
   e(2)=e(3);
   e(3)=s;
end
dfreq(1:2)=e(2:3);
%
 % compute event time as start of window
 
i1=round(wt(1)*Fs);  %Start index of event time
i2=round(wt(2)*Fs);  %start index of event+noise time series
i3=round(wt(3)*Fs);  %end index of event+noise time series
i4=i1-2*Fs;         %start index of noise sample
i5=i1-Fs;           %end index of noise sample
i6=i3+Fs;           % 1 sec beyond end of event sample
duration=wt(3)-wt(1);
ize=[i1 i2 i3 i4 i5 i6];
% rotate directional channels into grid alignment
sinbr=sin(.0174533*brefd);
cosbr=cos(.0174533*brefd);
css=cosch(i2:i3);
sss=sinch(i2:i3);
ns=css*cosbr -sss*sinbr;
ew=(sss*cosbr +css*sinbr);
%filter out spurious frequencies
omf=om(i2:i3); %volts
nsf=ns;
ewf=ew;
%zchf=zch(i2:i3);
%length(nsf)
%length(omf)
% measure rms value of omni channel during sample
vrms=std(omf);
dbvrms=20*log10(vrms);
% design and apply filter to pass only between e1 and e2
w1=2*e(2)/Fs;
w2=2*e(3)/Fs;
[B,A]=butter(2,[w1 w2]);
omf=filter(B,A,omf);  %volts
nsf=filter(B,A,nsf);
ewf=filter(B,A,ewf);
%zchf=filter(B,A,zchf);
% listen to filtered selection
omstr=resample(omf,8000,1000);
soundsc(omstr,8000);
%disp('Sound Fs is')
%Fs
time5=ts;
time5(6)= time5(6)+wt(1);
time5=normtv(time5);
mt=time5;
title5 = [DasarID ' at ' datestrb(time5)] ;
%---------------------------------------
% compute acoustic levels in the band and S/N ratio
Ts = 1/Fs;
t = Ts * (0:i6-i4)';
oml=10^(134/20)*om(i4:i6);
oml=filter(B,A,oml);    %muPa
%oml = filter(B,A,om(i4:i6,idb));  %volts
%oml=oml*10^(134/20); %convert volts to muPa(valid at 100 Hz)
%IIR filter to simulate flat hydrophone
% +2 db at 100 Hz (ideally 0 db)
% -8 db at 250 Hz  (ideally 20*log10(100/250)=-8
b = [0.3984  0.3984];
a = [1  -.9920];
oml = oml-mean(oml); %muPa
oml = filter(b,a,oml);  %muPa
% Make quick & dity plot of selected, filtered time series
%figure (11);
%plot(t,oml),grid;
%pt =[DasarID,Blanks(1),'Oml through inverse hyd filter',... 
%        Blanks(1),datestrb(ts)];
%title(pt);
rmsn=std(oml(1:i5-i4+1));
rmsspn = std(oml(i2-i4+1:i3-i4+1));
sigsq = max([(rmsspn^2 - rmsn^2) rmsn^2]);  %Don't permit s/n<0
sigdb=10*log10(sigsq);
stndb = sigdb -20*log10(rmsn);
% Make plot of selected, filtered time series
if moreplots
	figure (8);
	plot(t,oml),grid;
	pt =[DasarID Blanks(1) 'Selected, Filtered Time Series Starting:' ... 
            Blanks(1) datestrb(ts)];
	title(pt);
	ylabel ('MicroPascals');
	xlabel('Sec after Start-2');
end
%----------------------------------------
% Do psd on selected times
figure(5)
N = 512;
df=500/N;
psdsig=psd(om(i2:i3),N,Fs,hanning(N),N/2); %plots PSD
psd(om(i2:i3),N,Fs,hanning(N),N/2); %plots PSD
title(title5)
%get freq of max amplitude in selected band
ifstart=round(e(2)/500*length(psdsig));
ifend=  round(e(3)/500*length(psdsig));
[mxval,ifmax]=max(psdsig(ifstart:ifend));
fmax=(ifmax+ifstart)/length(psdsig)*500;
if ~moreplots
    close(5)
end
% plot omni versus vertical channel
% figure( 25)
% plot(omf,zchf)
% axis equal
% xlabel('Omni channel, Volts')
% ylabel('Vert channel, Volts')
%title('Filtered Vert Chan vs Omni')
%-------------------------------------
tic    %start stopwatch
% plot cosine vs sine signals
figure(21)
set(gca,'Fontsize',14);
pewf=ewf.*omf;
pnsf=nsf.*omf;
%pzchf=zchf.*omf;
% Plot either as 
%   a)point-to-point, preserving time relations, or
%   b)Origin-to-point-to-origin (looks cooler,takes longer).
lp2p=1;
if lp2p
    plot(pewf,pnsf);
else
    for i=1:length(pewf)
        xpl=[0 pewf(i)];
        ypl=[0 pnsf(i)];
        plot(xpl,ypl,'k-')
    end
end
axis equal
grid on
set(gca,'dataaspectratio',[1,1,1])
title([title5 ' (Horiz)']);
xlabel('EWF*OM v^2')
ylabel('NSF*OM v^2')
%
% get energy in components
sumew=sum(pewf);
sumns=sum(pnsf);
%sumzch=sum(pzchf);
% form weighted sums
vsin=mean(pewf);
vcos=mean(pnsf);
%vzch=mean(pzchf);
sigsin=std(pewf);
sigcos=std(pnsf);
%sigzch=std(pzchf);
bscat=bnorm(57.29577951*atan2(vsin,vcos));
%elscat=atan(vzch/sqrt(vsin^2+vcos^2))/rpdg;
%tscat=toc   %measure basic time for scatterplot
%if bscat < 0
%   bscat=bscat+360.;
%end
sigbscat=57.29577951*(vcos*sigsin-vsin*sigcos)/(vcos*vcos+vsin*vsin);
Sumtxt=['Average Grid Brg bscat=' num2str(bscat,4) ]  %'   Elev=' num2str(elscat,4)];
plabelfs(.05,.15,0,Sumtxt,12);
b=bscat;
%
%
% set up and call Chris Nations bearing scatter routine
tic     %start stopwatch
xy(:,1)=pewf;
xy(:,2)=pnsf;
%save xy.mat xy
[mu,kappa,sd]=get_vmests(xy,500);
%mu=bnorm(mu)
%kappa
%sd
%tchris=toc    %how long for Chris's routines
% rotate data until average direction coincides with +Y axis
tic    %start stopwatch for Bill's scatter
rangle=rpdg*bscat;
cosb=cos(rangle);
sinb=sin(rangle);
%c=[ cosb  -sinb; sinb  cosb];
xp = cosb*pewf -sinb*pnsf;
yp = sinb*pewf +cosb*pnsf;
ybar=mean(yp);
xbar=mean(xp);
xsig=std(xp);
ysig=std(yp);
qs=ysig/xsig;
qm=ybar/xsig;
bspread=180/pi*atan(xsig/ybar);
tbill=toc;      %stop watch for Bill's scatter
plabelfs(.05,.15,1,['Spread=' num2str(bspread,4) 'Deg'],12)
plabelfs(.05,.15,2,['kappa=' num2str(kappa,4)],12)
plabelfs(.05,.15,3,['sd=' num2str(sd,4)],12)
% % make similar scatter plot for y and z axis
%  figure(22)
%  plot(yp,pzchf)
%  axis equal
%  grid on
%  title([title5 ' (Vert)']);
%  xlabel('Horiz*OM V^2');
%  ylabel('PZCH*om V^2');
% % make 3d scatter plot
% figure(23);
% plot3(pewf,pnsf,pzchf);
%

% %if moreplots
% 	figure(23)
% 	clf
% 	plot(xp,yp,'b-',pzchf,yp,'r-')
% 	axis equal
%     grid on
% 	v=axis;
% 	dx=v(2)-v(1);
% 	dy=v(4)-v(3);
% 	text(v(1)+.1*dx,v(3) +.9*dy,['qs=' num2str(qs)])
% 	text(v(1)+.1*dx,v(3)+.85*dy,['qm=' num2str(qm)])
% 	text(v(1)+.1*dx,v(3)+.80*dy,['spread=' num2str(bspread) ' deg'])
% 	title(['Rotated Scatterplot for ' title5])
%     %end
%
% do further bearing processing if wanted
%
%do CSD in interval of unfiltered data
n=128;
if moreplots
    figure(6)
	psd(om(i1:i2),n,Fs,hanning(n),n/2) %plots PSD
	title(title5);
end
ww = ones(n,1);
OX = csd(om(i2:i3),ew,n,Fs,ww);
OY = csd(om(i2:i3),ns,n,Fs,ww);
rb0=bnorm(57.2957*atan2(real(OX),real(OY)));
rb1=dnorm(57.2957*atan2(real(OX),real(OY)));
sig0 = std(rb0);
sig1 = std(rb1);
bscatp = bscat;
if sig0 > sig1
    rb0 = rb1;
    bscatp = dnorm(bscat);
end
lox=length(OX);
df=.5*Fs/lox;
fp=df*(1:lox);
% compute bcsd as average angle in frequency span
jflo = ceil(e(2)/500*lox)+3;
jfhi = floor(e(3)/500*lox)-3;
 ifmax=round(fmax/500*lox);
 bcsd=bnorm(mean(rb0(jflo:jfhi)));
 %
%Put in range 0-360
%bbscat=0.;
%if rb0(ifmax)<0 |  rb0(ifmax)>360
%    rb0=bnorm(rb0);
%end
%for i=1:lox;
%    if rb0(i)<0;
%       rb1(i)=rb0(i)+360;
%    else
%       rb1(i)=rb0(i);
%    end
% end
 % final version has smallest sigma
 %rb = rb1;
 %if std(rb1) > std(rb0)
 %   rb=rb0;
 %   bbscat=-360.;
 %end
 disp(['Bcsd= ' num2str(bcsd)]);
 %
% plot bearing vs frequency
if moreplots
    figure(7)
	plot(fp,rb0,'b',[e(2) e(3)],[bscatp bscatp],'r-o',fp(jflo:jfhi),rb0(jflo:jfhi),'k.');
	title(title5)
	xlabel('Frequency (Hz)')
	ylabel('UTM Grid Bearing (degrees)')
    text(fp(10),rb0(10),['Mean bcsd= ' num2str(bcsd)])
    grid on
end
% return scatter plot to front
figure(21)
