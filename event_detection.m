load('data_out.mat')
a=tdfread('cities_NA.txt');


%Create yearly maximum distribution for each cluster (Cell array);
%cred is index of clusters
for i=1:numel(cred)
  % n_years*n_gridcells array of yearly max events, for each cluster i.
indist{i}=max_in(:,isgd(inclust{i}));
%work with array for cluster i
M=indist{i};
%find maximum event in region for the whole timeseries
[C,I] = max(M(:));
%find reference for the event in the "M" array
[I1,I2] = ind2sub(size(M),I);
M(I1,I2);
%make array n_years*n_gridcells of dates of yearly max events
dttmp=day_in(:,isgd(inclust{i}));
%find date of the cluster overall max event
dt_clust(i)=dttmp(I);
%find reference for the event in the "aa" lat/lon  array
[llo(i),lla(i)]=ind2sub(size(aa),isgd(inclust{i}(I2)));
%find reference for the year
lly(i)=I1;
%save actual lat/lon of the event
lat_cl(i)=clats(lla(i));
lon_cl(i)=clons(llo(i));
%calculate distance to each city in the database
dd=distance(lat_cl(i),lon_cl(i),a.lat,a.lon);
%find city index nearest to event
 [duf ddmin]=min(dd);
%construct city name
 citynm{i}=[strtrim(a.name1(ddmin,:)) ', ' strtrim(a.g(ddmin,:)), ...
            ', ' strtrim(a.e(ddmin,:))];
%save event in inches
c_lvl_in(i)=C*0.0393701;
%save event in mm
c_lvl(i)=C;
%fit GEV to overall distribution for the cluster, over the whole timeseries
 [parmhat(i,:),parmci] = gevfit(M(:));
%save the percentile of the event in the distribution
 c_prctl(i) = gevcdf(c_lvl(i),parmhat(i,1),parmhat(i,2),parmhat(i,3));
end
%plot the events and the regions
plotmap(lats,lon,aa,bb,cmap,'p_events.png',1,flipud(parula),lon_cl,lat_cl,dt_clust,citynm)

%load the CMIP data from read_cmip.m
load('/project/cmip5/ETH/ens_pr.mat')

%set counter to zero
nn=0;
%set temperature thresholds
t_thresh=[1.0,1.5,2,3,4];
clear prmat_*
for i=1:numel(ens)
  %use only models with the required date range
if and(min(ens(i).rcp85.prmx.year)<=1850,max(ens(i).rcp85.prmx.year)>=2100)
%increment index
  nn=nn+1;
  ensmap(nn)=i;
%  Save model name
  mdl{nn}=ens(i).name;
  %find time index of start/end points for past and future sims
  mnyr=find(ens(i).rcp85.prmx.year==1988);
  mxyr=find(ens(i).rcp85.prmx.year==2013);
  
    
  %save present day mat
  prmat_cmip(:,:,:,nn)=ens(i).rcp85.prmx.val(:,:,mnyr:mnyr+25);
  
  %get GMTS (20 yr smoothing)
  mdl_gmts=smooth(ens(i).rcp85.tas.val,20);
%get Pre-industrial average
mdl_pi=mean(mdl_gmts(1:30));

  for j=1:numel(t_thresh)
    %find first date of exceedance of GM temp thresh
  idx=min(find(mdl_gmts-mdl_pi>t_thresh(j)));
  if isempty(idx)
    idx=NaN;
  end
  
  rng=[idx-13:idx+12];
if rng(end)<numel(mdl_gmts)
  prmat_fut(:,:,:,nn,j)=ens(i).rcp85.prmx.val(:,:,rng);
  cntnan=find(~isnan(squeeze(mean(mean(prmat_fut(:,:,:,nn,j),1),2))));
  if numel(cntnan)==numel(rng)-1
  prmat_fut(:,:,:,nn,j)=ens(i).rcp85.prmx.val(:,:,[rng(1)-1,rng(cntnan)]);  
  end
  
else
  prmat_fut(:,:,:,nn,j)=NaN;
end
end
end
end
clear cmip_parmhat* indist_fut cmip_thresh cmip_fut_prctl
for i=1:numel(cred)
  for j=1:numel(mdl)
    pst_in=permute(prmat_cmip(:,:,:,j),[3,1,2,4]);
    fut_in=permute(prmat_fut(:,:,:,j,:),[5,3,1,2,4]);
    % n_years*n_gridcells array of yearly max events, for each cluster i.
    indist_cmip{i,j}=pst_in(:,isgd(inclust{i}));
    for k=1:numel(t_thresh)
    indist_fut{i,j,k}=squeeze(fut_in(k,:,inclust{i}));
    end
    %work with array for cluster i
    M=indist_cmip{i,j};
    [cmip_parmhat(i,j,:),cmip_parmci] = gevfit(M(:));
%map threshold for extreme event in this cluster (i) in this model (j)
%cmip_thresh(i,j)=gevinv(0.99,cmip_parmhat(i,j,1),cmip_parmhat(i,j,2),cmip_parmhat(i,j,3));
%cmip_fut_prctl(i,j,1)=gevcdf(cmip_thresh(i,j),cmip_parmhat(i,j,1),cmip_parmhat(i,j,2),cmip_parmhat(i,j,3));
    
    for k=1:numel(t_thresh)
    M=indist_fut{i,j,k};
    if ~isnan(mean(M))
    [cmip_parmhat_fut(i,j,k,:),cmip_parmci] = gevfit(M(:));
    if k==1
      cmip_thresh(i,j)=gevinv(c_prctl(i),cmip_parmhat_fut(i,j,1,1),cmip_parmhat_fut(i,j,1,2),cmip_parmhat_fut(i,j,1,3));
cmip_thresh_1000(i,j)=gevinv(0.999,cmip_parmhat_fut(i,j,1,1),cmip_parmhat_fut(i,j,1,2),cmip_parmhat_fut(i,j,1,3));
end
    
  else
    cmip_parmhat_fut(i,j,k,:)=NaN;
  end
  
if ~isnan(cmip_thresh(i,j))	  
    cmip_fut_prctl(i,j,k)=gevcdf(cmip_thresh(i,j),cmip_parmhat_fut(i,j,k,1),cmip_parmhat_fut(i,j,k,2),cmip_parmhat_fut(i,j,k,3));
    cmip_fut_1000(i,j,k)=gevcdf(cmip_thresh_1000(i,j),cmip_parmhat_fut(i,j,k,1),cmip_parmhat_fut(i,j,k,2),cmip_parmhat_fut(i,j,k,3));
  else
    cmip_fut_prctl(i,j,k)=NaN;
    cmip_fut_1000(i,j,k)=NaN;
  end
  
    end
    
end
end

%wlevel=[1,t_thresh];
alllvls=find(~isnan(squeeze(mean(mean(cmip_fut_prctl,1),3))));


for i=1:numel(cred)
  for k=1:numel(t_thresh)
    p_dist=cmip_fut_prctl(i,:,k);
    p_dist_fut(i,k,:)=prctile(p_dist,[5,25,50,75,95]);
    p_dist_1000=cmip_fut_1000(i,:,k);
    p_dist_fut_1000(i,k,:)=prctile(p_dist_1000,[5,25,50,75,95]);
  end
end


p_rtn=1./(1-p_dist_fut);
p_rtn_1000=1./(1-p_dist_fut_1000);
plot_cmip
plotmap_delta(lats,lon,aa,bb,cmap,'delta',4,parula,p_rtn_1000)


   