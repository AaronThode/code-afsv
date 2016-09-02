function [ctimes,tol]=load_debug_ctimes(Icase)
%Icase=1;
ctimes=[];
tol=2;
switch Icase,
    case 1
        %ctimes=1219359776.399784564971924;
        %ctimes=1188352193.598999977112;
        %ctimes=1219276833.999996662139893;
        %ctimes=1219276882.963999986649;
    %ctimes=(1.219277057838725e+09)-0.5;
%     ctimes=1219277397.671640396118164;  %Airgun and whale downsweep!
%     ctimes=1219276905.498782396316528;
%     ctimes=1219348845.977227210998535;
%     ctimes= 1.219896110443037e+09;
    %ctimes=datenum(2008,10,1,3,11,39);
    %ctimes=mat2c_tm(ctimes);
    case 2
        %ctimes=86400*(datenum(2008,8,21,00,1,06)-datenum(1970,1,1,0,0,0));
        %ctimes= 1.219276993175728e+09;
        %ctimes=-3.0+1.219277114778138e+09;
    case 3
        %ctimes=1188349031.954999923706;
end