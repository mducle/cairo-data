function hf = plot_bfo_in20(spec)

ym=30; % y-max
cm=5;
cfm=1/4;
logplot=0;

if nargin > 0
  sw_plotspec(sw_egrid(sw_neutron(spec),'Evect',0:0.2:100),'mode','color','dE',1); 
end
hold all;
legend('off'); title(''); ylabel(''); xlabel('');
load 'in20_fitted_peaks.mat';
plotted=zeros(1,length(setindex{3}));
for ist=1:length(setindex{3})
  hkl=sscanf(setindex{3}{ist},'(%f %f %f)');
  ihkl=sum(abs(cumhkl(:,1:3)-(repmat(hkl',length(cumhkl),1)))'); 
  ihkl=find(ihkl==min(ihkl));
  if(sum(abs(cumhkl(ihkl(1),1:3)-hkl'))<1e-1)
   for iu=unique(cumhkl(ihkl,4))'
    for ie=1:(length(fits{ist}{1})-1)/3
     if(abs(fits{ist}{1}((ie-1)*3+1))<1e7 && fits{ist}{1}((ie-1)*3+2)<ym)
      plotted(ist)=1;
      plot(iu,fits{ist}{1}((ie-1)*3+2),'ok','MarkerSize',abs(fits{ist}{1}((ie-1)*3+1))/2e5);
      errorbar(iu,fits{ist}{1}((ie-1)*3+2),fits{ist}{1}((ie-1)*3+3),'.k');
     end
    end
   end
  end
end
%pp = parula; pp(1,:) = [1 1 1]; colormap(pp)
set(gca,'XTickLabel',[]);

xlim([0 3.514]);                              % (h-0.5,k,-0.5)
plot([1 1]*0.5232,get(gca,'YLim'),'-k');      % (h-0.5,k,0)
plot([1 1]*0.9206,get(gca,'YLim'),'-k');      % (h,k,0)
plot([1 1]*1.293,get(gca,'YLim'),'-k');       % (h,k-0.5,0)
plot([1 1]*1.691,get(gca,'YLim'),'-k');       % (h-0.5,k-0.5,0)
plot([1 1]*2.236,get(gca,'YLim'),'-k');       % (h,k,0)
plot([1 1]*2.991,get(gca,'YLim'),'-k');       % (h-0.5,k-0.5,-0.5)
set(gca,'XTick',[0 0.5232 0.9206 1.293 1.691 2.236 2.991]); % (h-0.5,k-0.5,0)
set(gca,'XTickLabel',[]);
set(gca,'FontSize',14);
box on

ffo=ym/50; ofs = -3.8*ffo; oft = ofs-7*ffo; 
text(0     ,ofs,'$(h-\frac{1}{2},k,-\frac{1}{2})$', 'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(0.5232,oft,'$(h-\frac{1}{2},k,0)$',            'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(0.9206,ofs,'$(h,k,0)$',                        'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(1.293 ,oft,'$(h,k-\frac{1}{2},0)$',            'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(1.691 ,ofs,'$(h-\frac{1}{2},k-\frac{1}{2},0)$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(2.236 ,oft,'$(h,k,0)$',                        'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(2.991 ,ofs,'$(h-\frac{1}{2},k-\frac{1}{2},-\frac{1}{2})$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(3.514 ,oft,'$(h-\frac{1}{2},k-\frac{1}{2},0)$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);

ylb=ylabel('Energy Transfer (meV)','FontSize',18); yps=get(ylb,'Position'); set(ylb,'Position',[yps(1) ym yps(3)]); 
set(gca,'FontSize',14);
