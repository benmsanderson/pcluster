datadir='/project/cmip5/ETH/cmip5/historical/day/pr'
t_datadir='/project/cmip5/ETH/cmip5/historical/Amon/tas'

scen={'rcp85','rcp45'};




load('data_out.mat')
mdls=dir([datadir '/*']);
longrid=lon(aa);
latgrid=lats(bb);


for i=3:3%numel(mdls)
  runs=dir([datadir '/' mdls(i).name]);
  for j=1:numel(scen)
datadir_fut=['/project/cmip5/ETH/cmip5/' scen{j} '/day/pr'];
t_datadir_fut=['/project/cmip5/ETH/cmip5/' scen{j} '/Amon/tas'];

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
       [loni,lati,timi] = ndgrid(mod(lontmp,360),lattmp,ttmp);
       [lono,lato,timo] = ndgrid(mod(clons,360),clats,ttmp);
       vq = interpn(loni,lati,timi,prtmp,lono,lato,timo,'linear');
      prcon=cat(3,prcon,vq);
      tcon=cat(1,tcon,dates);
    end

    %%now look for matching future runs
  runs_fut=dir([datadir_fut '/' mdls(i).name]);
  
    files2=dir([datadir_fut '/' mdls(i).name '/' runs_fut(3).name]);
   try
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
       [loni,lati,timi] = ndgrid(mod(lontmp,360),lattmp,ttmp);
       [lono,lato,timo] = ndgrid(mod(clons,360),clats,ttmp);
       vq = interpn(loni,lati,timi,prtmp,lono,lato,timo,'linear');
      prcon=cat(3,prcon,vq);
      tcon=cat(1,tcon,dates);
    end	
    catch
    disp('error getting future precip')
  end
  
  %save extreme values
      datevc=datevec(tcon);
      cyrs=unique(datevc(:,1));
      for k=1:numel(cyrs)
	usedates=find(datevc(:,1)==cyrs(k));
	prmx(:,:,k)=max(prcon(:,:,usedates),[],3);
      end
      ens(i-2).(scen{j}).prmx.val=prmx;
      ens(i-2).(scen{j}).prmx.year=cyrs;
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

    catch
    disp('error getting future temperature')
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
    
    
    
    
    end
end



    
    
	    
     