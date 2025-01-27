function [c_new cmap]=make_clusters(nc,D1,minsize,tmp,tmpd,meanmat,uselats,uselons,isgd,aa)


%define normalization parameter for physical distance - lower values mean clausters will be more spatially confined
%D1=0.0025;

%produce a combined distance vector (both physical and parametric)
tmpmap=NaN(size(aa(:),1),size(aa(:),1));
for i=1:numel(isgd)
tmpmap(isgd(i),isgd)=tmp(i,:)+tmpd(i,:)/D1;
end

%%CLUSTER ANALYSIS                                                                                                                                                                                                                          
%initial number of clusters (nc - free parameter)
%nc=60;

%compute linkage and initial clustering
Z = linkage((tmp+tmpd/D1).^0.5,'complete');
c = cluster(Z,'maxclust',nc);



%define minimum number of gridcells in a cluster
%minsize=40;
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
[duf drysort]=sort(cmean);
%sort clusters from dryest to wettest region
for i=1:numel(uni)
c_new(find(c_org==uni(drysort(i))))=i;
end


%map onto lat/long grid for plotting
 cmap=NaN(size(meanmat(uselons,uselats,1)));
cmap(isgd)=c_new;


