function [row col]=findcoords(lat,lon,latarr,lonarr)

errmat=(latarr-lat).^2+(lonarr-lon).^2;
minMatrix = min(errmat(:));
[row,col] = find(errmat==minMatrix);
