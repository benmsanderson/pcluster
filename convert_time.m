function outnum=convert_time(time,units,calendar)
%units
conyear=regexp(units,'\d\d\d\d','match');
if isempty(conyear)==1
    conyear={'1'};
  end
   conyear=str2num(conyear{:});
   condat=regexp(units,'\d+.\d+.\d+','match');
   startdate = sscanf(condat{1},'%f-%f-%f');
   
   if length(startdate)>1
     conmon=startdate(2);
   else
     conmon=1
   end
   
   if length(startdate)>2
     conday=startdate(3);
   else
     conday=1
   end
   
   
   
   
   
   if any([strmatch(calendar,'noleap')==1,strmatch(calendar,'365_day')==1])
     n_days=365;
     smon=[31,28,31,30,31,30,31,31,30,31,30,31];
     yearf=((time+conyear*n_days+sum(smon(1:conmon))-smon(conmon)+conday)/n_days);
     year=floor(yearf);
     monthf=((yearf-year)*12)+1;
     month=floor(monthf);
     day=floor((monthf-month)*30+conday);
     outtime=[year,month,day,time*0,time*0,time*0];
   end
   
   if (strmatch(calendar,'360_day')==1)
     n_days=360;
     yearf=((time+conyear*n_days+(conmon-1)*30+conday)/n_days);
     year=floor(yearf);
     monthf=((yearf-year)*12)+1;
     month=floor(monthf);
     day=floor((monthf-month)*30+conday);
     outtime=[year,month,day,time*0,time*0,time*0];
   end
   
   if  any([strmatch(calendar,'standard')==1,strmatch(calendar,'gregorian')==1,strmatch(calendar,'proleptic_gregorian')==1])
     
     %  disp('Gregorian');
     starttime=datenum(conyear,1,1);
     abstime=time+starttime;
     outtime=datevec(abstime);
   end
   
   outnum=datenum(outtime);
   return
   