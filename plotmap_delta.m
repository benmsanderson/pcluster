function plotmap_delta(lats,lon,aa,bb,cmap,filename,n,colmap,p_rtn)


for i=1:size(p_rtn,2)
  cmap_ret(i,:,:)=cmap;
for j=1:size(p_rtn,1)
  ccrd=find(cmap==j);
  cmap_ret(i,ccrd)=p_rtn(j,i,3);
end
end

load('colwet.mat')
fig = figure(n)
clf
psn={[0.,0.55,0.5,0.5],...
    [0.5,0.55,0.5,0.5],...
    [0.,0.3,0.5,0.5],...
    [0.5,0.3,0.5,0.5]};

tits={'(a) 1.5C','(b) 2.0C','(c) 3.0C','(d) 4.0C'};
for i=1:size(p_rtn,2)-1
subplot(2,2,i)
hma=gca;
        
ax =usamap('conus');
set(ax,'Position',psn{i})
        setm(ax, 'Frame', 'off', 'Grid', 'off',...
	                'ParallelLabel', 'off', 'MeridianLabel', 'off')
delete(hma) 
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
faceColors = makesymbolspec('Polygon',...
    {'INDEX', [1 numel(states)], 'FaceColor','none'}); %NOTE - colors are random
 geoimg=geoshow(lats(bb),lon(aa),log10(squeeze(cmap_ret(i+1,:,:))), 'DisplayType', 'texturemap');
geoimg.AlphaDataMapping = 'none'; 
geoimg.FaceAlpha = 'texturemap'; 
alpha(geoimg,double(~isnan(squeeze(cmap_ret(i+1,:,:)))));

hold on
geoshow(ax, states, 'DisplayType', 'polygon', ...
   'SymbolSpec', faceColors);
 load coast
%geoshow(flipud(lat),flipud(long),'DisplayType','polygon','FaceColor','white');
colormap(colwet);
caxis([log10(50),3])
xl=get(gca,'XLim');
yl=get(gca,'YLim');
t=text((xl(1)+xl(2))/2,yl(2),tits{i});
set(t,'horizontalalignment','center','fontweight','bold')

end

cb=colorbar('horiz')
set(cb,'position',[0.2,0.4,0.6,0.025])
set(get(cb,'xlabel'),'string','Return Period of Present day 1000 year event (years)')

set(cb,'Ticks',[log10(50),2,log10(200),log10(500),3])
set(cb,'Ticklabels',[50,100,200,500,1000])

set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);

print(gcf,'-dpng','-r600','-painters',filename);


unix(['convert ' filename '.png -trim ' filename '.png'])