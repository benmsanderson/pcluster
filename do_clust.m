%datadir='/glade/p/rda/data/ds728.3/p1d/netcdf'
%load GPCC data https://rda.ucar.edu/datasets/ds497.0/
datadir='/glade/p/rda/data/ds497.0/full_data_daily_v1'

%switch - fit distribution or just use percentiles (default)
fitdist=0;


%read daily data
fils=dir([datadir '/*.nc']);
pmat=[];
tmat=[];

disp('Reading files...');
upd = textprogressbar(numel(fils),'updatestep',1);

for i=1:numel(fils)

upd(i);

tmp=ncread([datadir '/' fils(i).name],'p');
ttmp=ncread([datadir '/' fils(i).name],'time');
if i==1
    lon=ncread([datadir '/' fils(i).name],'lon');
    lats=ncread([datadir '/' fils(i).name],'lat');
end
pmat=cat(3,pmat,tmp);
tmat=cat(1,tmat,ttmp);
end


%%pmat (mm)  360 (lon)        180 (lat)       9497 (days)
year=cumsum([0;diff(tmat)<0]);

%%create date vector in datenum format (begins 1988)
date0=datenum([1988+year,ones(size(year)),ones(size(year))])+tmat;
%date0=datenum([1988,1,1])+tmat;

%full matlab datew vector
datevv=datevec(date0);

%create 3 day rain accumulation matrix
p3daysum=cumsum(pmat,3);
mm3day=cat(3,p3daysum(:,:,4:end)-p3daysum(:,:,1:end-3),NaN([size(p3daysum,1),size(p3daysum,2),3]));


%list of years
fyears=unique(datevv(:,1));



%%%Create matrix of monthly maximum and mean values, by year, for extreme value analysis
maxmat=NaN(size(p3daysum,1),size(p3daysum,2),1,12);
disp('Getting extreme values')
upd = textprogressbar(numel(fyears),'updatestep',1);
for i=1:numel(fyears)
upd(i);
for m=1:12
usedates=find(and(datevv(:,1)==fyears(i),datevv(:,2)==m));
if numel(usedates)>0
maxmat(:,:,i,m)=max(mm3day(:,:,usedates),[],3);
meanmat(:,:,i,m)=mean(pmat(:,:,usedates),3);
else
maxmat(:,:,i,m)=NaN;
meanmat(:,:,i,m)=3;
end
end
end

disp('Making parameter matrix')
upd = textprogressbar(size(pmat,1),'updatestep',1);
for i=1:size(pmat,1)
upd(i);
%loop over longitude
for j=1:size(pmat,2)
    %loop over latitude
if ~isnan(mean(pmat(i,j,:),3));
for m=1:12
    %find dates in month m
usedates=find(datevv(:,2)==m);
%find all days when it rained more than 0.1mm
rainydays=usedates(find(squeeze(pmat(i,j,usedates))>0.1));
%determined rainy day frction
pfrac=numel(rainydays)/numel(pmat(i,j,usedates));
%form distribution from rainy days
raindist=squeeze(pmat(i,j,rainydays));
%%(for fitdist) define 2 parameter exponential function

pdf_expmixture = @(x,p,mu1,mu2) ...
                         p*exppdf(x,mu1) + (1-p)*exppdf(x,mu2);


%guess for mixing ratio of 2 exponential distributions
pStart = .5;
%define the two exponential parameters as 45th and 85th percentile
muStart = quantile(raindist,[.45 .85]);

if isnan(mean(muStart))
    %default to very low values for case when it never rained in the whole dataset (Sahara, some months)
muStart = [0.05,0.1];
end

%%for fitdist options - make the first parameter less then 60th percentile, 2nd parameter >60th percentile
muBnd = quantile(raindist,[.6 .99]);
%lower and upper parameter bounds
lb = [0 0 muBnd(1)];
ub = [1 muBnd(1) muBnd(2)];
%starting vector (mixing ratio, parameter first guesses
start=[pStart muStart ];
if fitdist
    %fitdist option, just cluster based on mixed exponential parameter values (SLOW!!!)
paramEsts = mle(raindist, 'pdf',pdf_expmixture, 'start',start, ...
                          'lower',lb, 'upper',ub);
else
    %with no fitdist option, just cluster based on percentile values
paramEsts=start;
end

%full parameter matrix  - 4 for each gridcell - rain fraction, mixing fraction (always 0.5 if fitdist=0), exponential parameters x2
prmmat(i,j,m,:)=[pfrac,paramEsts];
end
else
    %if the original p data is NaN (ocean, probably), make this Nan too
prmmat(i,j,m,:)=NaN(1,4);
end
end
end



%%make 360x180 lat/lon vectors
latarr=repmat(lats',numel(lon),1);
lonarr=repmat(lon,1,numel(lats));

%% Now define inter-point distances
%first - phyiscal distance
distfun=@(x1,x2) distance(latarr(x1(1),x1(2)),lonarr(x1(1),x1(2)),latarr(x2(1),x2(2)),lonarr(x2(1),x2(2)));
%now - distance in parameter 1 (rain fraction)
p2fun=@(x1,x2) sum((prmmat(x1(1),x1(2),:,1)-prmmat(x2(1),x2(2),:,1)).^2,3);

%to calibrate parameter normalization - find distance between each point and its neighbour
for i=1:size(prmmat,1)-1
for j=1:size(prmmat,2)
tmp2(i,j)=p2fun([i,j],[i+1,j]);
end
end
%normalization factor - 5th percentile of this distribution 
D2=prctile(tmp2(:),5);

%now - distance in parameter 3 (percentile value or exp parameter 1)
p3fun=@(x1,x2) sum((prmmat(x1(1),x1(2),:,3)-prmmat(x2(1),x2(2),:,3)).^2,3);
%to calibrate parameter normalization - find distance between each point and its neighbour
for i=1:size(prmmat,1)-1
for j=1:size(prmmat,2)
tmp3(i,j)=p3fun([i,j],[i+1,j]);
end
end
%normalization factor - 20th percentile of this distribution 
D3=prctile(tmp3(:),20);

%now - distance in parameter 4(percentile value or exp parameter 2)
p4fun=@(x1,x2) sum((prmmat(x1(1),x1(2),:,4)-prmmat(x2(1),x2(2),:,4)).^2,3);
%to calibrate parameter normalization - find distance between each point and its neighbour
for i=1:size(prmmat,1)-1
for j=1:size(prmmat,2)
tmp4(i,j)=p4fun([i,j],[i+1,j]);
end
end
%normalization factor - 20th percentile of this distribution 
D4=prctile(tmp4(:),20);

%combine all parametric distances into single function
combidist=@(x1,x2) p2fun(x1,x2)/D2+p3fun(x1,x2)/D3+p4fun(x1,x2)/D4;

%define boundaries for US domain
[minlon,minlat]=findcoords(20.1,-132.1,latarr,lonarr);
[maxlon,maxlat]=findcoords(50.1,-60.1,latarr,lonarr);

uselats=[minlat:maxlat];
uselons=[minlon:maxlon];
%aa and bb are the indices of points in the US domain
[aa bb]=ndgrid(uselons,uselats);
%make vector (isgd) of all non-nan points in the US domain
 isgd=find(~isnan(pmat(uselons,uselats,1)));

clear tmp tmpd


disp('Calculating inter-point distances...');
upd = textprogressbar(numel(isgd),'updatestep',1);
for ii=1:numel(isgd)
upd(ii);
%define inter-point distances for all point pairs
for jj=1:numel(isgd)
    %parametric distance
tmp(ii,jj)=combidist([aa(isgd(ii)),bb(isgd(ii))],[aa(isgd(jj)),bb(isgd(jj))]);
%physical distance
tmpd(ii,jj)=distfun([aa(isgd(ii)),bb(isgd(ii))],[aa(isgd(jj)),bb(isgd(jj))]);
end
end


%define normalization parameter for physical distance - lower values mean clausters will be more spatially confined
D1=0.0025;

%produce a combined distance vector (both physical and parametric)
tmpmap=NaN(size(aa(:),1),size(aa(:),1));
for i=1:numel(isgd)
tmpmap(isgd(i),isgd)=tmp(i,:)+tmpd(i,:)/D1;
end



%%CLUSTER ANALYSIS
%initial number of clusters (nc - free parameter)
nc=60;

%compute linkage and initial clustering
Z = linkage((tmp+tmpd/D1).^0.5,'complete');
c = cluster(Z,'maxclust',nc);



%define minimum number of gridcells in a cluster
minsize=40;
%duplicate cluster vector before loop
c_org=c;

%%Form coalesced clusters
%initilize minimum cluster varaible
minclust=0;

%continue loop while there exists a cluster smaller than minsize
while minclust<minsize
    %create vector of unique clusters remaining
cred=unique(c_org);

clear inclust numelc
for i=1:numel(cred)
    %make list of cells in each cluster
inclust{i}=(find(c_org==cred(i)));
%count members in each cluster
numelc(i)=numel(inclust{i});
end

%find smallest cluster - size is minclust, end if minclust>minsize
[minclust ac]=min(numelc);

%work with smallest cluster, at index ac
elimc=ac;
%list elements in that cluster
cands=inclust{elimc};
%total clusters vector
alc=[1:numel(cred)];
clear ctarg
for j=1:numel(cands)
    %loop through gridcells in the cluster
clear ctmp
for k=1:numel(cred)
    %find the distance between the candidate gridcell and the mean position of each remaining cluster
ctmp(k)=mean(mean(tmp(cands(j),inclust{k})+tmpd(cands(j),inclust{k})/D1,1),2);
end
%list the clsuters which are not the elimination cluster
allbut=find(alc~=elimc);
%find the cluster which is closest to the candidate point (outside of the original elimination cluster)
[duf ttmp]=min(ctmp(alc(allbut)));
%set the candidate's destination cluster
ctarg(j)=allbut(ttmp);
end
%remap all the elements of the candidate cluster to their new destination
c_org(cands)=cred(ctarg);
end


%list remaining clusters
uni=unique(c_org);

%find mean precipitation in each cluster
 us_mean=mean(mean(meanmat(uselons,uselats,:,:),3),4);
 %map onto landmasked isgd vector
us_mngd=us_mean(isgd);


for i=1:numel(uni)
cmean(i)=mean(us_mngd(find(c_org==uni(i))));
end
[duf drysort]=sort(cmean)
%sort clusters from dryest to wettest region
for i=1:numel(uni)
c_new(find(c_org==uni(drysort(i))))=i;
end


%map onto lat/long grid for plotting
 cmap=NaN(size(pmat(uselons,uselats,1)));
cmap(isgd)=c_new;



%plot clusters
plotmap(lats,lon,aa,bb,cmap,'pmaps.png')

%plot demonstration distance map
[blon,blat]=findcoords(40.01,-105.3,lats(bb),lon(aa));
gdmap=NaN(size(pmat(uselons,uselats,1)));
gdmap(isgd)=1:numel(isgd);
bgood=gdmap(blon,blat);
dmap=reshape(log(1./tmpmap(isgd(bgood),:)),size(aa));
plotmap(lats,lon,aa,bb,dmap,'ptest.png')


%%Create time evolving cluster distributions
