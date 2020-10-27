function bearings=bnorm(bearings)

while any(bearings(:)>=360)
    Ibig=find(bearings>=360);
    bearings(Ibig)=bearings(Ibig)-360;
end

while any(bearings(:)<0)
    Ilow=find(bearings<0);
    bearings(Ilow)=bearings(Ilow)+360;
end

end