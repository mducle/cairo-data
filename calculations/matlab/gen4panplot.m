function hf = gen4panplot(spec, bfopowspec)

hf = figure; 
s{1}=subplot(141); 
try
    sw_plotspec(sw_egrid(sw_neutron(spec),'Evect',0:0.2:100),'mode','color','dE',1); 
catch
    sw_plotspec(spec); 
end
caxis([0 1])
%tt=title(sprintf('J_{c1}=%4.2f J_{c2}=%4.2f J_{ab1}=%4.2f J_{ab2}''=%4.2f J_{d}=%4.2f',squeeze(spec.obj.matrix.mat(3,3,1:5))));
pars = squeeze(spec.obj.matrix.mat(3,3,:))
if numel(pars)==6
    tt=title(sprintf('J_{c1}=%4.2f J_{c2}=%4.2f J_{ab1}=%4.2f J_{ab2}''=%4.2f J_{d}=%4.2f D=%4.2f',pars));
elseif numel(pars)==7
    tt=title(sprintf('J_{c1}=%4.2f J_{c2}=%4.2f J_{ab1}=%4.2f J_{ab2}''=%4.2f J_{d}=%4.2f D_{12}=%4.2f D_{3}=%4.2f',pars));
else
    error('bad pars');    
end
set(tt, 'Position', [5.7517,100.9617,100.0000]);
colorbar off
if nargin > 1
    s{2}=subplot(142); %sw_plotspec(bfopowspec,'dE',1); caxis([0 5]); title('');
    %ispe = sw_instrument(bfopowspec, 'dE', 2.9, 'dQ', 0.01, 'Ei', 66.4, 'ThetaMin', 3.5);
    eis = [16.6 66.4];
    des = [0.65 2.9];
    eis = [16.33 66 160];
    %des = [0.37 1.8 7.4];
    %des = [0.8 1.8 7.4];
    des = [1.1 2.8 7.4];
    iqm = [1 1 1];
    iqx = [1.5 4 4];
    ispe = sw_instrument(bfopowspec, 'dE', des(2), 'dQ', 0.01, 'Ei', 166.4, 'ThetaMin', 3.5);
    sw_plotspec(ispe); caxis([0 0.5]); title('');
    colorbar off
    s{3}=subplot(143); hold all;
    load('bfof_cut_dat');
    mar160 = load('mar_ei160_cut'); errorbar(mar160.x, mar160.y*100, mar160.e*100, 'or');
    mar66 = load('mar_ei66_cut'); errorbar(mar66.x, mar66.y*100+0.2, mar66.e*100, 'sb');
    %dt = cutinx(out2p22, 'x', 1, 1.5); errorbar(dt.x, dt.y*143, dt.e*143, '-k'); 
    sf = 0.2; ofs = [0.025 0]; %0.05;
    lg = {};
    cc = 'krb';
    eis = [160 66 16.6]; des = [7.4 1.8 1.1]; iqm = [1 1 1]; iqx = [4 4 1.5]; cc = 'rbk';
    bgp = {[1.4832 10.5033 117.5923 0.9293] [0.5974 2.4936 142.83 0.669]};
    bgp{1} = [1.4936 10.4833 117.9228 0.9180  20.6089 30.3615 31.972/(29.16/4) 0.66  47.3719 53.4852 27.6894/(29.16/4) 0.14];
    bgp{2} = [0.5974 2.4936*0.8 142.83 0.669  9.267 7.457 3.4521/(11.56/3.5) 0.24  20.056 14.648 10.9651/(11.56/3.5) 0.99  38.0319 22.9696 9.2782/(11.56/3.5) 0.612];
    %bgp{1} = bgp{2};
    %iqm = [1 2 1]; iqx = [1.5 2.5 1.5]; eis = [66 66 16.6]; des = [1.8 1.8 1.1];
    iqm = [1 2 1]; iqx = [4 2.5 1.5]; eis = [166 66 16.6]; des = [7.4 1.8 1.1];
    for ii = 1:2 % 2:3;
        ispe = sw_instrument(bfopowspec, 'dE', des(ii), 'dQ', 0.01, 'Ei', eis(ii), 'ThetaMin', 3.5);
        iq0=max(find(ispe.hklA<=iqm(ii)));
        iq1=max(find(ispe.hklA<=iqx(ii)));
        sqw=ispe.swConv; sqw=sqw(:,iq0:iq1);
        isqw=sqw; isqw(find(~isnan(sqw)))=1; isqw(find(isnan(sqw)))=0; sqw(find(isnan(sqw)))=0; 
        xx = (ispe.Evect(2:end)+ispe.Evect(1:end-1))/2;
        yy = sum(sqw')./sum(isqw')+(ii-1)*sf+ofs(ii);
        bg = fvoigt(xx, bgp{ii})';
        plot(xx, yy+bg, cc(ii))
        plot(xx, bg, ['--' cc(ii)]);
        plot(xx, yy, [':' cc(ii)]);
    end
    ylim([0 2]); xlim([0 100]);
    box on; xlabel('Energy Transfer (meV)'); ylabel('Intensity (arb. units)');
    %legend({'\lambda=2.21', '\lambda=1.11'}, 'location', 'NorthEast');
    s{4} = subplot(144); hold all;
    dt = cutinx(out2p22, 'x', 1, 1.5); errorbar(dt.x, dt.y*143, dt.e*143, '-k'); 
    ii = 3; sf = 0; ofs = 0.05;
        ispe = sw_instrument(bfopowspec, 'dE', des(ii), 'dQ', 0.01, 'Ei', eis(ii), 'ThetaMin', 3.5);
        iq0=max(find(ispe.hklA<=iqm(ii)));
        iq1=max(find(ispe.hklA<=iqx(ii)));
        sqw=ispe.swConv; sqw=sqw(:,iq0:iq1);
        isqw=sqw; isqw(find(~isnan(sqw)))=1; isqw(find(isnan(sqw)))=0; sqw(find(isnan(sqw)))=0; 
        plot((ispe.Evect(2:end)+ispe.Evect(1:end-1))/2, (sum(sqw')./sum(isqw'))*1.7*1.27+(ii-1)*sf+ofs, cc(ii))
    box on;
    xlim([0 12]); ylim([0 5])
end
cm = viridis; cm(1,:) = [1 1 1]; 
colormap(s{1}, viridis);
colormap(s{2}, viridis)
set(s{1}, 'position', [0.13  0.11 0.213 0.815]);
if nargin > 1
set(s{2}, 'position', [0.411 0.11 0.213 0.815]);
set(s{3}, 'position', [0.692 0.11 0.213 0.815]);
set(s{4}, 'position', [0.755 0.44 0.125 0.382]);
end
set(gcf, 'PaperOrientation', 'landscape')
set(gcf, 'PaperPosition', [0.6 0.6 28.5 19.7])
