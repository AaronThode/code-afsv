function bearings=bnorm(bearings)
Ibig=find(bearings>=360);
bearings(Ibig)=bearings(Ibig)-360;

Ilow=find(bearings<0);
bearings(Ilow)=bearings(Ilow)+360;

end