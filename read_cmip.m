datadir='/project/cmip5/ETH/cmip5/historical/day/pr'
t_datadir='/project/cmip5/ETH/cmip5/historical/Amon/tas'

scen={'rcp85','rcp45'};




load('data_out.mat')
mdls=dir([datadir '/*']);
longrid=lon(aa);
latgrid=lats(bb);


for i=16:numel(mdls)
disp(mdls(i).name);
  runs=dir([datadir '/' mdls(i).name]);
  for j=1:numel(scen)
datadir_fut=['/project/cmip5/ETH/cmip5/' scen{j} '/day/pr'];
t_datadir_fut=['/project/cmip5/ETH/cmip5/' scen{j} '/Amon/tas'];
try
    files=dir([datadir '/' mdls(i).name '/' runs(3).name]);
prcon=[];
tcon=[];
    for k=3:numel(files)
      fn=[datadir '/' mdls(i).name '/' runs(3).name '/' files(k).name];
      prtmp=ncread(fn,'pr');
      lattmp=ncread(fn,'lat');
      lontmp=ncread(fn,'lon');
      ttmp=ncread(fn,'time');
      torg=ncread(fn,'time');
      tunits=ncreadatt(fn,'time','units');
      tcal=ncreadatt(fn,'time','calendar');
      dates=convert_time(ttmp,tunits,tcal);
       [loni,lati] = ndgrid(mod(lontmp,360),lattmp);
       [lono,lato] = ndgrid(mod(clons,360),clats);
     clear vq
       for ii=1:size(prtmp,3)
       vq(:,:,ii) = interpn(loni,lati,prtmp(:,:,ii),lono,lato,'linear');
     end
     
      prcon=cat(3,prcon,vq);
      tcon=cat(1,tcon,dates);
    end
    pstname=runs(3).name;
catch
disp('error reading historical')
prcon=[];
tcon=[];
pstname='';
end

    %%now look for matching future runs
  runs_fut=dir([datadir_fut '/' mdls(i).name]);
   try
    files2=dir([datadir_fut '/' mdls(i).name '/' runs_fut(3).name]);
  
    for k=3:numel(files2)
      fn=[datadir_fut '/' mdls(i).name '/' runs_fut(3).name '/' files2(k).name];
      prtmp=ncread(fn,'pr');
      lattmp=ncread(fn,'lat');
      lontmp=ncread(fn,'lon');
      ttmp=ncread(fn,'time');
      torg=ncread(fn,'time');
      tunits=ncreadatt(fn,'time','units');
      tcal=ncreadatt(fn,'time','calendar');
      dates=convert_time(ttmp,tunits,tcal);
       [loni,lati] = ndgrid(mod(lontmp,360),lattmp);
       [lono,lato] = ndgrid(mod(clons,360),clats);
clear vq
       for ii=1:size(prtmp,3)
	 
       vq(:,:,ii) = interpn(loni,lati,prtmp(:,:,ii),lono,lato,'linear');
       end
      prcon=cat(3,prcon,vq);
      tcon=cat(1,tcon,dates);
    end	
  rfname=runs_fut(3).name;
    catch
    disp('error getting future precip')
rfname='';
end
  p3daysum=cumsum(prcon,3);
  mm3day=cat(3,p3daysum(:,:,4:end)-p3daysum(:,:,1:end-3),NaN([size(p3daysum,1),size(p3daysum,2),3]));
  
  %save extreme values
      datevc=datevec(tcon);
      cyrs=unique(datevc(:,1));
prmax=[];
      for k=1:numel(cyrs)
	usedates=find(datevc(:,1)==cyrs(k));
	prmx(:,:,k)=max(mm3day(:,:,usedates),[],3);
      end
      ens(i-2).(scen{j}).prmx.val=prmx;
      ens(i-2).(scen{j}).prmx.year=cyrs;
      ens(i-2).(scen{j}).runs=pstname;
      ens(i-2).(scen{j}).runs_fut=rfname;
	    
      ens(i-2).name=mdls(i).name;
    
try

          %%now look for historical temperature
   files=dir([t_datadir '/' mdls(i).name '/' runs(3).name]);
   prcon=[];
   tcon=[];
        for k=3:numel(files)
      fn=[t_datadir '/' mdls(i).name '/' runs(3).name '/' files(k).name];
      prtmp=ncread(fn,'tas');
      lattmp=ncread(fn,'lat');
      lontmp=ncread(fn,'lon');
      ttmp=ncread(fn,'time');
      torg=ncread(fn,'time');
      tunits=ncreadatt(fn,'time','units');
      tcal=ncreadatt(fn,'time','calendar');
      dates=convert_time(ttmp,tunits,tcal);
      larr_tmp=cos(repmat(lattmp(:)',numel(lontmp),1)/180*pi);
      latwgt=larr_tmp/mean(larr_tmp(:));
      prwgt=squeeze(mean(mean(prtmp.*repmat(larr_tmp,1,1,size(prtmp,3)),1),2));
      
      prcon=cat(1,prcon,prwgt);
      tcon=cat(1,tcon,dates);
    end
catch
disp('error getting past temperture')
prcon=[];
tcon=[];
tsmn=[];
end

try    
              %%and future temperature
   files=dir([t_datadir_fut '/' mdls(i).name '/' runs_fut(3).name]);
        for k=3:numel(files)
      fn=[t_datadir_fut '/' mdls(i).name '/' runs_fut(3).name '/' files(k).name];
      prtmp=ncread(fn,'tas');
      lattmp=ncread(fn,'lat');
      lontmp=ncread(fn,'lon');
      ttmp=ncread(fn,'time');
      torg=ncread(fn,'time');
      tunits=ncreadatt(fn,'time','units');
      tcal=ncreadatt(fn,'time','calendar');
      dates=convert_time(ttmp,tunits,tcal);
      larr_tmp=cos(repmat(lattmp(:)',numel(lontmp),1)/180*pi);
      latwgt=larr_tmp/mean(larr_tmp(:));
      prwgt=squeeze(mean(mean(prtmp.*repmat(larr_tmp,1,1,size(prtmp,3)),1),2));
      
      prcon=cat(1,prcon,prwgt);
      tcon=cat(1,tcon,dates);
    end

  
  
 %save temperature values
      datevc=datevec(tcon);
      cyrs=unique(datevc(:,1));
clear tasmn
      for k=1:numel(cyrs)
	usedates=find(datevc(:,1)==cyrs(k));
	tasmn(k)=mean(prcon(usedates));
      end
      ens(i-2).(scen{j}).tas.val=tasmn';
      ens(i-2).(scen{j}).tas.year=cyrs;

    catch
    disp('error getting future temperature');
  end    
    
    
    
    end
  end
  
  save('/project/cmip5/ETH/ens_pr.mat','ens')