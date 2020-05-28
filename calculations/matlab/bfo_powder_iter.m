function chi2 = bfofit(Jpars, do_plot)

if nargin < 2
    do_plot = false;
end

if nargin < 1
    load('longpowspec.mat')
else
    bfo = bi2fe4o9_spinw(Jpars);
    bfopowspec = bfo.powspec(linspace(0,3,60), 'Evect', 0:0.2:100, 'nRand', 100, 'fibo', true, ...
                            'hermit', false, 'formfact', true, 'optmem', 50);
end

powspec = sw_instrument(bfopowspec, 'dE', 11.5, 'dQ', 0.05, 'Ei', 300, 'ThetaMin', 3.5);
swConv = powspec.swConv; hklA = powspec.hklA; Evect = powspec.Evect; 

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

bkg_cuts = load('bfo_cuts_bkg.mat');
fn = fieldnames(bkg_cuts);
for ii = 1:numel(fn)
    dat_cuts.(fn{ii}) = load([fn{ii} '.mat']);
end

scale_fac = [1000, 1000, 1000, 1000, 1000, 2500];
xmin = 10;
sf = 100;
if do_plot; figure; hold all; end
for ii = 1:numel(eis)
    lev = (ii-1)*sf;
    id_dat = find(dat_cuts.(fn{ii}).x > xmin);
    id_cal = find(bfo_powcut{ii}{1} > xmin);
    xbkg = bkg_cuts.(fn{ii}){1};
    xbkg = (xbkg(2:end)+xbkg(1:end-1))/2.;
    id_bkg = find(xbkg > xmin);
    ybkg = interp1(xbkg(id_bkg), bkg_cuts.(fn{ii}){2}(id_bkg), bfo_powcut{ii}{1}(id_cal));
    ycal = interp1(bfo_powcut{ii}{1}(id_cal), bfo_powcut{ii}{2}(id_cal)*scale_fac(ii)+ybkg, dat_cuts.(fn{ii}).x(id_dat));
    ydat = dat_cuts.(fn{ii}).y(id_dat);
    edat = dat_cuts.(fn{ii}).e(id_dat);
    id_nan = find(isnan(ycal));
    ydat(id_nan) = [];
    edat(id_nan) = [];
    ycal(id_nan) = [];
    %chi2(ii) = sum((ydat - ycal).^2 ./ ycal) / numel(ycal);  % runs 1-3
    chi2(ii) = sum((ydat - ycal).^2 ./ (edat.^2)) / numel(ycal);
    if do_plot
        errorbar(dat_cuts.(fn{ii}).x(id_dat), dat_cuts.(fn{ii}).y(id_dat)+lev, dat_cuts.(fn{ii}).e(id_dat), 'o');
        plot(bfo_powcut{ii}{1}(id_cal), bfo_powcut{ii}{2}(id_cal)*scale_fac(ii)+lev, '--');
        plot(bfo_powcut{ii}{1}(id_cal), bfo_powcut{ii}{2}(id_cal)*scale_fac(ii)+lev+ybkg, '-');
    end
end
if do_plot
    box on; xlabel('Energy Transfer (meV)'); ylabel('Intensity (arb. units)');
    legend(lg, 'location', 'SouthEast');
    ylim([0 700])
    %print('-dpdf', 'longpow_cuts.pdf')
end
%chi2
