function plotmap(lats,lon,aa,bb,cmap,filename,n,colmap,lon_cl,lat_cl,dt_clust,citynm)

fig = figure(n)
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
if nargin>8
plotm(lat_cl,lon_cl,'ko','markerfacecolor','w')
for i=1:numel(lat_cl)
tt=textm(lat_cl(i),lon_cl(i),{citynm{i},datestr(dt_clust(i))});
set(tt,'fontsize',5)
end
end
set(gcf, 'PaperPosition', [0 0 6 5]);
set(gcf, 'PaperSize', [6 5]);

print(gcf,'-dpng','-r600','-painters',filename);


