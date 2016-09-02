function tbai=ctime2byteoffset(cwant,did,tsamp,bpsin,bbo)
%byteoff-number of bytes offset from start of raw file acquistion, computed
%from cwant
if did.ctdstart<=cwant && cwant<did.ctdend && ~strcmp(did.use,'F')
    % compute time correction to sample clock
    dtoff=did.tint + (cwant-did.ctref)*did.tdrift/86400;
    ctbsamp=cwant+dtoff;
    %compute offset from file start rounded to nearest ms
    toffb=ctbsamp-did.ctstart;    %seconds and fractions
    tbai=round(tsamp*toffb)*bpsin +bbo;   %byte address in all raw files
    
end
end
