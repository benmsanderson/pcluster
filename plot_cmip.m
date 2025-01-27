
figure(3)
clf
for i=1:numel(cred)
  subplot(3,5,i)
   a1=area(t_thresh,[p_rtn_1000(i,:,2);p_rtn_1000(i,:,4)-p_rtn_1000(i,:,2)]');ylim([0,1/(1-0.998)*2]);
   set(a1(1),'visible','off')
   set(a1(2),'linestyle','none')
   set(a1(2),'facecolor',[0.8,0.8,0.8])
   hold on
   plot(t_thresh,p_rtn_1000(i,:,3),'k-','linewidth',2);

   plot(t_thresh,p_rtn_1000(i,:,1),'k:','linewidth',2);
   plot(t_thresh,p_rtn_1000(i,:,5),'k:','linewidth',2);
   plot(t_thresh,p_rtn(i,:,3),'r--','linewidth',2);
%   plot([1,4],[1,1]./(1-c_prctl(i)),'r:','linewidth',2);
title({['(' char(96+i) ') ' citynm{i}],datestr(dt_clust(i))});
ylabel('Return Period')
xlabel('\DeltaT_{global}')
set(gca,'YScale','log')
   ylim([10,30000]);
set(gca,'ytick',[10,100,1000,10000])
set(gca,'yticklabel',([10,100,1000,10000]))
grid on
xlim([1,4])
end
 
set(gcf, 'PaperPosition', [0 0 12 6]);
set(gcf, 'PaperSize', [12 6]);

print(gcf,'-dpdf','cmip_proj.pdf');

				
				