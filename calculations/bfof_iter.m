function [chi2,bfopowspec] = bfof_iter(pars, doplt, save_file_name)
% Runs one iteration of the fit of INS data on Bi4Fe5O13F.
%
% To run a fit, use:
%
% x0 = [1.3 3.7 2.9 6.3 24.0 0.1 0.5 -0.1];  % Parameters in Beauvois, Simonet et al.
% x_opt = fminsearch(@bfof_iter, x0)

if nargin < 3
    save_file_name = 'fit_par_tmp.mat';
end

% If temp save file not present write it with current parameters
if exist(save_file_name, 'file')
    x0 = load(save_file_name);
else
    x0 = struct;
end

if ~exist('bi4fe5o13f_spinw', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/matlab']);
end

mar66 = load('mar_ei66_cut');
mar160 = load('mar_ei160_cut');
load('bfof_cut_dat');
dt = cutinx(out2p22, 'x', 1, 1.5);
xd = {{mar160.x(1:150), mar160.y(1:150)*100, mar160.e(1:150)*100} ...               % MARI Ei=160meV data
      {mar66.x(1:120), mar66.y(1:120)*100, mar66.e(1:120)*100} ...                  % MARI Ei=66meV
      {dt.x(420:522), dt.y(420:522)*143/1.7/1.27/3, dt.e(420:522)*143/1.7/1.27}};

if nargin < 2
    doplt = 0;
end

if doplt
    figure; hold all;
end

if iscell(pars)
    bfopowspec = pars;  % Spectrum already calculated, just calculate chi and/or plot.
    for ii=1:5; Jvar(ii) = pars{1}.obj.matrix.mat(1,1,ii); end
    bfo = pars{1}.obj;
else 
    bfo = bi4fe5o13f_spinw(pars, 2); % Runs SpinW model
    Jvar = pars(1:5);
end

ofs = [0.025 0 0.05/3];  % Flat background level
ifs = [1 1 1/3];         % Intensity scaling factor
lg = {};
% Parameters for background function (series of pseudo-voigts, fitted from high-Q data)
bgp{1} = [1.4936 10.4833 117.9228 0.9180 ...
          20.6089 30.3615 31.972/(29.16/4) 0.66 ...
          47.3719 53.4852 27.6894/(29.16/4) 0.14];
bgp{2} = [0.5974 2.4936*0.8 142.83 0.669 ...
          9.267 7.457 3.4521/(11.56/3.5) 0.24 ...
          20.056 14.648 10.9651/(11.56/3.5) 0.99 ...
          38.0319 22.9696 9.2782/(11.56/3.5) 0.612];
bgp{3} = [0 1 0 0];
eis = [166 66 16.6];
des = [7.4 1.8 1.1];
iqm = [1 2 1];
iqx = [4 2.5 1.5];
cc = 'rbkr';
pp = 'os^+';
nsp=[10 10 10];
nr = [50 50 100];

chi2 = 0;
for ii = 1:3
    if doplt
        subplot(4,1,ii); 
        hold all;
        errorbar(xd{ii}{1}, xd{ii}{2}, xd{ii}{3}, [pp(ii) cc(ii)]);
        if ii == 2
            title(sprintf('J_{c1}=%4.2f J_{c2}=%4.2f J_{ab1}=%4.2f J_{ab2}=%4.2f J_{d}=%4.2f', Jvar(1:5)));
        end
    end
    if ~iscell(pars)
        bfopowspec{ii} = bfo.powspec(linspace(iqm(ii),iqx(ii),nsp(ii)), 'Evect',0:0.2:100, 'nRand', nr(ii), ...
                                     'fibo',true, 'hermit',false, 'formfact',true, 'optmem',0, 'imagChk', -1);
    end
    ispe = sw_instrument(bfopowspec{ii}, 'dE', des(ii), 'dQ', 0.1, 'Ei', eis(ii), 'ThetaMin', 3.5);
    iq0=max(find(ispe.hklA<=iqm(ii)));
    iq1=max(find(ispe.hklA<=iqx(ii)));
    sqw=ispe.swConv*ifs(ii); sqw=sqw(:,iq0:iq1);
    isqw=sqw; isqw(find(~isnan(sqw)))=1; isqw(find(isnan(sqw)))=0; sqw(find(isnan(sqw)))=0; 
    xx = (ispe.Evect(2:end)+ispe.Evect(1:end-1))/2;
    yy = ( sum(sqw')./sum(isqw') );
    % Simulations only down to zero meV, for negative energies need to extend x-axis for background
    xn = -20:0.01:0; xx = [xn xx]; yy = [xn*0+ofs(ii) yy];
    bg = fvoigt(xx, bgp{ii})' + ofs(ii);
    if doplt
        plot(xx, yy + bg, cc(ii))
        plot(xx, bg, ['--' cc(ii)]);
        plot(xx, yy, [':' cc(ii)]);
        ylim([0 1]);
        box on; xlabel('Energy Transfer (meV)'); ylabel('Intensity (arb. units)');
        if ii == 3
            ylim([0 0.5]); xlim([0 20]);
        end
    end
    y1 = interp1(xx, yy + bg, xd{ii}{1}, 'linear', 0);
    vchi2(ii) = (sum((xd{ii}{2} - y1).^2 ./ xd{ii}{3}.^2) / numel(y1)) * (1 + bfopowspec{ii}.ioMax);
end

if doplt
    subplot(414);
    q0 = [3 3 0];
    if ~iscell(pars)
        bfopowspec{4} = bfo.spinwave({[-0.5 0 -0.5]+q0 [-0.5 0 0]+q0 [0 0 0]+q0 [0 -0.5 0]+q0 [-0.5 -0.5 0]+q0 ...
                              [0 0 0]+q0 [-0.5 -0.5 0.5]+q0 [-0.5 -0.5 0]+q0 20},'hermit',false,'optmem',10);
    end
    try
        sw_plotspec(sw_egrid(sw_neutron(bfopowspec{4}),'Evect',0:0.2:100.0),'mode','color','dE',1)
    catch
        sw_plotspec(bfopowspec{4});
    end
    ylim([0 100]); caxis([0 5]);
end

chi2 = mean(vchi2);

disp_struct = struct;
disp_struct.x = pars;
disp_struct.chi2_vec = vchi2;
disp_struct.chi2 = chi2;
disp_struct.ang13 = min(abs((acos(dot(bfo.magstr.S(:,[4]), bfo.magstr.S(:,[17]))/5) * 180 / pi) + [0 -180]));
disp_struct.angs = bfo.cache.angs;
disp_struct

chi2 = chi2 / 10;

if isfield(x0, 'xvec')
    xvec = [x0.xvec disp_struct];
else
    xvec = [disp_struct];
end
save(save_file_name, 'xvec');
