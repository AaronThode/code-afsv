%%%Review_changes.m

fnames=dir('*.m')
for I=1:length(fnames),
    eval(sprintf('!diff %s Old/%s >> tmp.txt',fnames(I).name,fnames(I).name));
end