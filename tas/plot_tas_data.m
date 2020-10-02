close all;

[thisfolder, ~, ~] = fileparts(mfilename('fullpath'));

if ~exist('illdata'); addpath([thisfolder '/load']); end
if ~exist('cividis'); addpath([thisfolder '/../calculations/matlab']); end

% Some parameters for the plots
ym = 30;   % Max energy (y-axis)
cfm = 1/4; % Max intensity ratio (color axis is to this fraction of max intensity)
logplot=0; % Plot color log scale?

if ~exist('eigdat', 'var') || ~exist('in20dat', 'var')
    % Load Eiger data
    [x, y, e, setindex] = eigermerge(4213:4348,'eiger_data');
    % Rebins data to 0.5meV steps from 1 to 30meV
    en = 1:0.5:29.5; 
    for ie=1:length(en); 
        for ist=3:22; 
            ief=find(abs(x{ist}-(en(ie)-0.25))<1e-3); 
            ieg=find(abs(x{ist}-en(ie))<1e-3); 
            ys{ist}(ie)=sum(y{ist}([ief ieg])); 
            ms{ist}(ie)=sum(e{ist}([ief ieg])); 
        end;
    end
    % Applies the Bose factor correction
    for ist=3:22;
        ysb{ist} = (ys{ist}./ms{ist}*1e6) .* (1-exp(-en*(11.604519339/setindex{6}(ist))));
        esb{ist} = (sqrt(ys{ist})./ms{ist}*1e6) .* (1-exp(-en*(11.604519339/setindex{6}(ist))));
    end
    % Store in a struct
    eigdat.x = x; eigdat.y = ysb; eigdat.e = esb; eigdat.s = setindex;
    
    % Load IN20 data
    [x, y, mn, setindex] = flexmerge(91215:91415, [thisfolder '/in20_data'], 3e5);
    in20dat.x = x; in20dat.y = y; in20dat.e = e; in20dat.s = setindex;
    
    % Loads fitted positions
    fits_in20 = load([thisfolder '/../calculations/matlab/in20_fitted_peaks.mat']);
    
    % Load Bi2Fe4O9 fitted spin wave calculations
    load([thisfolder '/../calculations/longpowspec.mat']);

    % Loads fitted positions from Eiger data
    fits_eiger = load([thisfolder '/fits_eiger.mat']);
end

% Plots the calculated spectrum at the top
hf = figure; 
s1 = subplot(221); hold all;
sw_plotspec(sw_egrid(sw_neutron(spec),'Evect',0:0.2:100),'mode','color','dE',1); 
hold all;
legend('off'); colorbar off; title(''); ylabel(''); xlabel('');
fits = fits_in20.fits; setindex = fits_in20.setindex; cumhkl = fits_in20.cumhkl;
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

xlim([0 3.514]); ylim([0 ym]);%box on;        % (h-0.5,k,-0.5)
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

% Polarised neutron data points
lpa = {[2.5 3 0 4 21] [2.5 2.5 0 15] [2.5 2.5 -0.5 4]};
lpax = []; lpaid = [];

% Plots the data below
s2 = subplot(223); hold all;
ma = 0;
xsb = in20dat.x; ysb = in20dat.y;
for ist=1:length(in20dat.s{3})
  hkl=sscanf(in20dat.s{3}{ist},'(%f %f %f)');
  if isempty(hkl); continue; end
  ihkl=sum(abs(cumhkl(:,1:3)-(repmat(hkl',length(cumhkl),1)))'); 
  ihkl=find(ihkl==min(ihkl));
  if(sum(abs(cumhkl(ihkl(1),1:3)-hkl'))<1e-1)
   for iu=unique(cumhkl(ihkl,4))'
     plotted(ist)=1;
     try
      if(logplot==1)
        pcolor([-0.025 0.025]+iu,xsb{ist},repmat(log10(ysb{ist}+1),2,1)'); shading flat;
        if(ma<max(log10(ysb{ist}+1))); ma=max(log10(ysb{ist}+1)); end
      else
        pcolor([-0.025 0.025]+iu,xsb{ist},repmat(ysb{ist},2,1)'); shading flat;
        if(ma<max(ysb{ist})); ma=max(ysb{ist}); end
      end
     end
     % Plots polarised data points
     hkld = cellfun(@(hkle)sum(abs(hkl-hkle(1:3)')), lpa);
     if min(hkld) < 1e-2
        lpax = [lpax iu]; lpaid = [lpaid find(hkld==min(hkld))];
     end
   end
  end
end
if(logplot==1)
    caxis([0.8 3])
else
    caxis([1 ma*cfm]);
end

for idd = 1:numel(lpax)
  plot(lpax(idd), lpa{lpaid(idd)}(4), 'ro', 'MarkerSize', 15, 'LineWidth', 2)
  if numel(lpa{lpaid(idd)}) > 4
    plot(lpax(idd), lpa{lpaid(idd)}(5), 'ro', 'MarkerSize', 15, 'LineWidth', 2)
  end
end


xlim([0 3.514]); ylim([0 ym]);%box on;        % (h-0.5,k,-0.5)
plot([1 1]*0.5232,get(gca,'YLim'),'-k');      % (h-0.5,k,0)
plot([1 1]*0.9206,get(gca,'YLim'),'-k');      % (h,k,0)
plot([1 1]*1.293,get(gca,'YLim'),'-k');       % (h,k-0.5,0)
plot([1 1]*1.691,get(gca,'YLim'),'-k');       % (h-0.5,k-0.5,0)
plot([1 1]*2.236,get(gca,'YLim'),'-k');       % (h,k,0)
plot([1 1]*2.991,get(gca,'YLim'),'-k');       % (h-0.5,k-0.5,-0.5)
set(gca,'XTick',[0 0.5232 0.9206 1.293 1.691 2.236 2.991]); % (h-0.5,k-0.5,0)
set(gca,'XTickLabel',[]);

ffo=ym/50; ofs = -3.8*ffo; oft = ofs-7*ffo; 
%text(0     ,ofs,'$(h-\frac{1}{2},k,-\frac{1}{2})$', 'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(0.5232,oft,'$(h-\frac{1}{2},k,0)$',            'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(0.9206,ofs,'$(h,k,0)$',                        'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(1.293 ,oft,'$(h,k-\frac{1}{2},0)$',            'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(1.691 ,ofs,'$(h-\frac{1}{2},k-\frac{1}{2},0)$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(2.236 ,oft,'$(h,k,0)$',                        'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(2.991 ,ofs,'$(h-\frac{1}{2},k-\frac{1}{2},-\frac{1}{2})$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(3.514 ,oft,'$(h-\frac{1}{2},k-\frac{1}{2},0)$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(0     ,ofs,'$(\frac{5}{2},3,-\frac{1}{2})$',    'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(0.5232,oft,'$(\frac{5}{2},3,0)$',               'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(0.9206,ofs,'$(3,3,0)$',                         'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(1.293 ,oft,'$(3,\frac{5}{2},0)$',               'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(1.691 ,ofs,'$(\frac{5}{2},\frac{5}{2},0)$',     'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(2.236 ,oft,'$(3,3,0)$',                         'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(2.991 ,ofs,'$(\frac{5}{2},\frac{5}{2},-\frac{1}{2})$','HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(3.514 ,oft,'$(\frac{5}{2},\frac{5}{2},0)$',     'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);

ylb=ylabel('Energy Transfer (meV)','FontSize',18); yps=get(ylb,'Position'); set(ylb,'Position',[yps(1) ym yps(3)]); 
set(gca,'FontSize',14);
box on

% Plots the Eiger spin wave calculation
s3 = subplot(222);
sw_plotspec(sw_egrid(sw_neutron(eigspec),'Evect',0:0.2:100),'mode','color','dE',1); 
fits = fits_eiger.fits; setindex = fits_eiger.setindex;
hold all;
kboth=[1.5 1.55 1.6:0.1:2 1:0.1:1.4];
for ist=[18:-1:14 7:13]
  for ie=1:(length(fits{ist}{1})-1)/3
    if(abs(fits{ist}{1}((ie-1)*3+1))<1e4 && fits{ist}{1}((ie-1)*3+2)<ym)
      xx = kboth(ist-6) - 1; if xx < 0.5; xx = xx * 0.5421/0.5; else xx = (xx - 0.5)*(1.066-0.5421)/0.5 + 0.5421; end
      plot(xx,fits{ist}{1}((ie-1)*3+2),'ok','MarkerSize',abs(fits{ist}{1}((ie-1)*3+1)));
      errorbar(xx,fits{ist}{1}((ie-1)*3+2),fits{ist}{1}((ie-1)*3+3)/2.35,'.k');
    end
  end
end
legend('off'); colorbar off; title(''); ylabel(''); xlabel('');
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
plot([1 1]*0.5421, get(gca,'YLim'),'-k');
xlim([0 1.067]); ylim([0 ym]);

% Plots Eiger data on right panel
s4 = subplot(224);
xsb = eigdat.x; ysb = eigdat.y;
hold all; shading flat;
for ist=[18:-1:14 7:13]
  if(ist < 10); ff = 1; else; ff = 2; end
  pcolor([-0.02 0.02]*ff+kboth(ist-6),en,repmat(ysb{ist},2,1)'); shading flat; caxis([5 35]);
end
xlim([0.98 2.02]); ylim([0 30]); box on;
plot([1 1]*1.5, get(gca,'YLim'),'-k');
set(gca, 'XTick', [1 1.5 2]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
%text(1     ,ofs,'$(h-1,k-1,\frac{3}{2})$',                     'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(1.5   ,oft,'$(h-\frac{3}{2},k-\frac{3}{2},\frac{3}{2})$', 'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
%text(2.0   ,ofs,'$(h-\frac{3}{2},k-\frac{3}{2},2)$',           'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(1     ,ofs,'$(2,2,\frac{3}{2})$',                     'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(1.5   ,oft,'$(\frac{3}{2},\frac{3}{2},\frac{3}{2})$', 'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);
text(2.0   ,ofs,'$(\frac{3}{2},\frac{3}{2},2)$',           'HorizontalAlignment','Center','interpreter', 'latex','FontSize',12);

% Rearranges the panels so that the IN20 data covers 3/4 of the space
s1p = get(s1, 'Position');
s2p = get(s2, 'Position');
s3p = get(s3, 'Position');
s4p = get(s4, 'Position');
xpp = (s4p(3) + s4p(1) - s1p(1)) * 0.85; 
ypp = (s1p(4) + s1p(2) - s2p(2)) / 2 + s2p(2);
gap = 0.02;
set(s1, 'Position', [s1p(1) ypp+gap xpp-s1p(1) s1p(4)+s1p(2)-(ypp+gap)]);
set(s2, 'Position', [s2p(1) s2p(2) xpp-s1p(1) ypp-gap-s2p(2)]);
set(s3, 'Position', [xpp+gap ypp+gap s3p(3)+s3p(1)-(xpp+gap) s3p(4)+s3p(2)-(ypp+gap)]);
set(s4, 'Position', [xpp+gap s4p(2) s4p(3)+s4p(1)-(xpp+gap) ypp-gap-s4p(2)]);

cmap = flipud(cividis); cmap(1,:) = [1 1 1];
colormap(cmap); 
