rawfile='/Volumes/macmussel1/Arctic_2008/Site_05/S508A0/S508A0T20080821T000000.gsi';
[y,t,head]=readfile(rawfile,0,4,1,'ctime','no');
 %y=y+0.5*(2^16);  %Remove DC bias in A/D converter
   
disp('First values of y');
disp(y(1:10));


 %%This has a 10 Hz high pass filter
    filt.a=[1.000000000000000e+00    -2.911197067426073e+00     2.826172902227507e+00    -9.149758348014339e-01];
    filt.b=[5.140662826979191e-01    -9.510226229911504e-01     3.598463978885433e-01     7.710994240468787e-02];
    
    %     % oml = oml - mean(oml); % uPa     get rid of DC. not needed, filter has zero @ DC
    % x = (2.5/65535)* (10^(134/20))*filter(filt.b,filt.a,xin); % uPa     equalize
    amp_Scale = (2.5/65535)*(10^(148.8/20));
    %amp_Scale=1;
    x = filter(filt.b,filt.a,y);

    disp('filtered values of x')
    disp(x(1:10));

    x = amp_Scale*x;

    disp('scaled values of x')
    disp(x(1:10));
    
    x(1:1024)=hamming(1024).*x(1:1024);

    disp('windowed output')
    disp(x(1:50));
    
    X=fft(x,1024);
    disp ('FFT:')
    disp(X(1:5))