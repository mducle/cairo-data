%% Some initial parameters

force_recalculate = false; % Recalculate spectra? (takes 2-3 days)
do_plot = true;            % Make Matlab plots?

%% Sets paths

if ~exist('spinw', 'file')
    error('These scripts need SpinW installed. Please go to https://www.spinw.org and install it or use the "Add-Ons" menu');
end

if ~exist('bi2fe4o9_spinw', 'file') || ~exist('bi4fe5o13f_spinw', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/matlab']);
end

%% Does the long powder average SpinW calculations
% For Bi2Fe4O9
if force_recalculate || ~exist('longpowspec.mat', 'file')
    [bfo, spec] = bi2fe4o9_spinw([0.8 1.0 3.5 7.5 19]);
    bfopowspec = bfo.powspec(linspace(0,5,400),'Evect',0:0.2:100,'nRand',5000,'fibo',true,'hermit',false,'formfact',true,'optmem',50);
    save('longpowspec.mat', 'bfo', 'spec', 'bfopowspec')
end

% For Bi4Fe5O13F
if force_recalculate || ~exist('longpowspec.mat', 'file')
    [bfof, spec_bfof] = bi4fe5o13f_spinw([2.9984 0.4457 3.3665 9.9118 13.9041 0.1008 0.5664 -0.1647],2);
    bfofpowspec = bfof.powspec(linspace(0,5,100),'Evect',0:0.2:100,'nRand',5000,'fibo',true,'hermit',false,'formfact',true,'optmem',200);
    save('bfof_fit_3sia_phase2.mat', 'bfof', 'spec_bfof', 'bfofpowspec')
end

%% Generate data for Bi2Fe4O9 Python plots.
clearvars -except do_plot
load('longpowspec.mat')
powspec = sw_instrument(bfopowspec, 'dE', 11.5, 'dQ', 0.05, 'Ei', 300, 'ThetaMin', 3.5);
swConv = powspec.swConv; hklA = powspec.hklA; Evect = powspec.Evect; 
save('bfo_powspec.mat', 'swConv', 'hklA', 'Evect')
if do_plot
    figure; 
    subplot(121); 
    plot_bfo_in20(spec); caxis([0 50]);
    tt=title(sprintf('J_{44}=%4.2f J_c=%4.2f J_{43}=%4.2f J_{43}''=%4.2f J_{33}=%4.2f',squeeze(bfo.matrix.mat(3,3,1:5))));
    set(tt, 'Position', [5.7517,100.9617,100.0000]);
    subplot(122); %sw_plotspec(bfopowspec,'dE',1); caxis([0 1]); title('');
    sw_plotspec(powspec); caxis([0 0.2]); title('')
    set(gcf, 'PaperOrientation', 'landscape'); set(gcf, 'PaperPosition', [1 1 28 20]);
    %print('-dpdf', 'longpow.pdf');
end

eis = [25 38 62 120 180 300];
des = [0.7 1.1 2 5.6 9.2 11.4];
lg = {};
bfo_powcut = {};
for ii = 1:numel(eis)
    ispe = sw_instrument(bfopowspec, 'dE', des(ii), 'dQ', 0.01, 'Ei', eis(ii), 'ThetaMin', 3.5);
    iq=max(find(ispe.hklA<=3));
    sqw=ispe.swConv; sqw=sqw(:,1:iq);
    isqw=sqw; isqw(find(~isnan(sqw)))=1; isqw(find(isnan(sqw)))=0; sqw(find(isnan(sqw)))=0; 
    x = (ispe.Evect(2:end)+ispe.Evect(1:end-1))/2;
    y = sum(sqw')./sum(isqw');
    if ii < numel(eis)
        lg{ii} = ['Merlin Ei=' num2str(eis(ii))];
    else
        lg{ii} = ['MAPS Ei=' num2str(eis(ii))];
    end
    bfo_powcut{ii} = {x, y, lg{ii}};
end
save('bfo_powcut.mat', 'bfo_powcut');
if do_plot
    sf = 0.05;
    figure; hold all;
    for ii = 1:numel(eis)
        plot(bfo_powcut{ii}{1}, bfo_powcut{ii}{2}+(ii-1)*sf);
    end
    box on; xlabel('Energy Transfer (meV)'); ylabel('Intensity (arb. units)');
    legend(lg, 'location', 'SouthEast');
    %print('-dpdf', 'longpow_cuts.pdf')
end

%% Generate data for Bi4Fe5O13F Python plots
clearvars -except do_plot
load('bfof_fit_3sia_phase2.mat')
eis = [166, 66, 16.6];
des = [7.4 1.8 1.1];
powspec = sw_instrument(bfofpowspec, 'dE', des(2), 'dQ', 0.01, 'Ei', 166.4, 'ThetaMin', 3.5);
swConv = powspec.swConv; hklA = powspec.hklA; Evect = powspec.Evect; 
save('bfof_powspec.mat', 'swConv', 'hklA', 'Evect');

if do_plot
    figure; sqw = sw_plotspec(powspec); caxis([0 0.5]); title('');
end

bfof_powcuts = {};

sf = 0; 
ofs = [0.025 0 0.05];
lg = {};
cc = 'rbk';
bgp{1} = [1.4936 10.4833 117.9228 0.9180  20.6089 30.3615 31.972/(29.16/4) 0.66  47.3719 53.4852 27.6894/(29.16/4) 0.14];
bgp{2} = [0.5974 2.4936*0.8 142.83 0.669  9.267 7.457 3.4521/(11.56/3.5) 0.24  20.056 14.648 10.9651/(11.56/3.5) 0.99  38.0319 22.9696 9.2782/(11.56/3.5) 0.612];
bgp{3} = [0 0 0 0];
iqm = [1 2 1]; iqx = [4 2.5 1.5];
lg = {'MARI Ei=160 meV', 'MARI Ei=66 meV', 'IN4 $\lambda = 2.22\mathrm{\AA}$ (Ei=16 meV)'};
for ii = 1:3
    ispe = sw_instrument(bfofpowspec, 'dE', des(ii), 'dQ', 0.01, 'Ei', eis(ii), 'ThetaMin', 3.5);
    iq0=max(find(ispe.hklA<=iqm(ii)));
    iq1=max(find(ispe.hklA<=iqx(ii)));
    sqw=ispe.swConv; sqw=sqw(:,iq0:iq1);
    isqw=sqw; isqw(find(~isnan(sqw)))=1; isqw(find(isnan(sqw)))=0; sqw(find(isnan(sqw)))=0; 
    xx = (ispe.Evect(2:end)+ispe.Evect(1:end-1))/2;
    yy = sum(sqw')./sum(isqw')+(ii-1)*sf+ofs(ii);
    bg = fvoigt(xx, bgp{ii})';
    bfof_powcuts{ii} = {xx, yy, bg, lg{ii}};
end
save('bfof_powcut.mat', 'bfof_powcuts');
if do_plot
    figure; hold all;
    load('bfof_cut_dat.mat');
    mar160 = load('mar_ei160_cut'); errorbar(mar160.x, mar160.y*100, mar160.e*100, 'or');
    mar66 = load('mar_ei66_cut'); errorbar(mar66.x, mar66.y*100, mar66.e*100, 'sb');
    dt = cutinx(out2p22, 'x', 1, 1.5); errorbar(dt.x, dt.y*66, dt.e*66, '^k'); 
    for ii = 1:3
        xx = bfof_powcuts{ii}{1};
        yy = bfof_powcuts{ii}{2};
        bg = bfof_powcuts{ii}{3};
        plot(xx, yy+bg, cc(ii))
        plot(xx, bg, ['--' cc(ii)]);
        plot(xx, yy, [':' cc(ii)]);
    end
    ylim([0 2]); xlim([0 100]);
    box on; xlabel('Energy Transfer (meV)'); ylabel('Intensity (arb. units)');
end
