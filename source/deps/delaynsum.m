

%%%%%%%%%%%%%%%%delaynsum.m%%%%%%%%%%%%%%%%
%  Aaron Thode
%  October 17, 1996
%  Given an angle theta, element spacing d,
%  delay and sum an incoming signal
%  function xtot=delaynsum(x,thta,space,Fs,goodel,cc)
% Input:
%  x is [nsample nchan] matrix
%  theta is in degrees, where broadside is zero.  A scalar value
%   space is a scalar in units of meters
%   Fs is sampling rate in Hz
%   goodel is vector of indicies of x that are good
%   cc is sound speed in meters/sec
%   rows of x are samples, columns are channels
% Output:
%   xtot: beamformed time series, averaged by number of elements
%   Mel: reference element.

function [xtot,Mel]=delaynsum(x,thta,space,Fs,goodel,cc)
%disp('pause');pause
nx=size(x,1);
thta=thta*pi/180;
if ~exist('cc','var')
   cc=1500;
end

xtot=[];
if length(thta)~=1
    disp('Error, delaynsum requires scalar input for thta');
    return
elseif length(goodel)~=size(x,2)
    disp('Error, length(goodel) is not the number of channels in x');
    return
end

dpt=Fs*space*sin(thta)/cc; %sample shift per element.  Must be a fraction
fprintf('Fs is %6.2f\n',Fs);
Mel=ceil(.5*(max(goodel)+min(goodel)));  %Reference element: number that experiences zero shift
fprintf('Median element is %i\n',Mel);
maxn=ceil(max(abs(goodel-Mel))*abs(dpt))+1;  %%Maximum element sample shift possible
nx_out=nx-2*maxn+1;  %%%Number of samples in trimmed signal
xtot=zeros(nx_out,1);

for I=1:length(goodel) %From deep to shallow
	n=ceil(dpt*(Mel-goodel(I))); %n decreases toward surface
	index=maxn+n+(0:(nx-2*maxn));  %%%sample shift is zero at reference element (Mel)
	xtot=xtot+x(index,goodel(I));
end


%%%Average in order to keep same units as Mel
xtot=xtot/length(goodel);

%Insert appropriate zeros-lost a day because of this!
%%% Transisent in beamformed series now occr at same sample indicies as Mel
%%% reference

xtot=[zeros(maxn,size(xtot,2)); xtot];

%Append zeros to end so that time series the same length as original
xtot=[xtot; zeros(nx-size(xtot,1),1)];

return


%Makes a simulated plane wave from the surface plus reflection, plot results
if 1==0
  Fs=1500;space=1.875;goodel=1:48;
   % goodel=[1:24 26:31 33:36 38:41 43:48]
     a1=45
    a2=-30
    f1=800
    f2=900
    A1=1;
    A2=1;
    Ioff=48;  %Element around which original signal is centered
    d1=100; %delay in onsets of signals
    d2=300; %delay in onset of signal 2
    x=zeros(4098,48);
    time=2*pi*f1*(1:size(x,1))/Fs;
    puls=1;
    if puls==0
        for I=1:size(x,2)
            time0=time+(I-Ioff)*2*pi*f1*sin(a1*pi/180)*space/1500;
            x(d1:end,I)=A1*sin(time0(d1:end))';
            time1=2*pi*f2*(d2:size(x,1))/Fs+(I-Ioff)*2*pi*f2*sin(a2*pi/180)*space/1500;
            x(d2:end,I)=x(d2:end,I)+A2*sin(time1');
        end
    else
        for I=1:size(x,2)
            time0=d1+Fs*(Ioff-I)*sin(a1*pi/180)*space/1500;
            x(round(time0),I)=A1;
            time1=d2+Fs*(Ioff-I)*sin(a2*pi/180)*space/1500;
            x(round(time1),I)=-A2;
        end
    end
    
    figure
    subplot(3,1,1);
    plot(x(:,25));title(int2str(Ioff));grid;axis([0 500 -2 2]);hold on
    plot(x(:,1),'--');plot(x(:,48),'k');
    subplot(3,1,2);
    y1=delaynsum(x,a1,space,Fs,goodel,1500);
    plot(y1);hold;%plot(x(:,Ioff),'k');
    grid;axis([0 500 -2 2])
    subplot(3,1,3);
    y2=delaynsum(x,a2,space,Fs,goodel,1500);
    plot(y2);hold;%plot(x(:,Ioff),'k');
    grid;axis([0 500 -2 2])

end



