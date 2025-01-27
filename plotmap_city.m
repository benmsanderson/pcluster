function plotmap_city(lats,lon,aa,bb,cmap,filename,n,colmap,lon_cl,lat_cl)

figure(n)
clf
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
faceColors = makesymbolspec('Polygon',...
    {'INDEX', [1 numel(states)], 'FaceColor','none'}); %NOTE - colors are random
 geoshow(lats(bb),lon(aa),cmap, 'DisplayType', 'texturemap');

hold on
geoshow(ax, states, 'DisplayType', 'polygon', ...
   'SymbolSpec', faceColors);
 load coast
geoshow(flipud(lat),flipud(long),'DisplayType','polygon','FaceColor','white');
colormap(colmap);

plotm(lat_cl,loncl,'k*')

set(gcf, 'PaperPosition', [0 0 6 5]);
set(gcf, 'PaperSize', [6 5]);

print(gcf,'-dpng','-r600','-painters',filename);


