function topography = readImg(fileName, numCols)
% Read 16 bit SRTM30+ file, and keep it int16 to save memory.
fid = fopen(fileName,'r');
[topography, cnt] = fread(fid, inf, 'int16=>int16');
fclose(fid);
% make image a rectangle. SRTM30+ has 4800 columns north of 60S
% but 7200 columns south of 60S. So user has say how many?
numRows = cnt / numCols;
topography = reshape(topography, numCols, numRows)';
% On a ?big endian? CPU, (e.g., Intel Mac), swap bytes.
topography = swapbytes(topography);
end
