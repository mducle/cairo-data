function fitStr = fit_in20_spinw(Jvals, Jvary, do_plot)

if nargin < 3
    do_plot = false;
end

% Generate the text file for fitspec if it doesn't exist
if ~exist('in20_fitted_modes.txt', 'file')
    gen_in20_spinw_dat
end

Jvals = [1.3 3.7 2.9 6.3 24.0]; % Initial parameters from Beauvois, Simonet et al. PRL 124 127202 '19

% Sets up model
bfo = bi2fe4o9_spinw(Jvals);
Jnames = bfo.matrix.label(cellfun(@(ss) strncmp(ss, 'J', 1), bfo.matrix.label));

% If certain parameters are to be fixed
xmax = [5 5 10 10 20 1];
if nargin > 1
    J0 = Jvals(find(Jvary));
    parnames = Jnames(find(Jvary));
    xmax = [xmax(find(Jvary)) 1];
else
    J0 = Jvals;
    parnames = Jnames;
end

% Sets up fitting
par_fit           = struct;
par_fit.datapath  = 'in20_fitted_modes.txt';
par_fit.Evect     = linspace(5,30,51);
par_fit.func      = @(obj, p) matparser(obj, 'param', p, 'mat', {parnames{:} 'K(3,3)'}, 'init', true);
par_fit.xmin      = zeros(1, numel(J0) + 1);
par_fit.xmax      = xmax;
par_fit.x0        = [J0 0.02];
par_fit.plot      = do_plot;
par_fit.hermit    = false;
par_fit.optimizer = 'simplex';
par_fit.maxiter   = 200;
%par_fit.optimizer = 'pso';
%par_fit.maxiter   = 20;
par_fit.optmem    = 1;
par_fit.nrun      = 1;
%par_fit.imagChk   = 'penalize';

fitStr = bfo.fitspec(par_fit);
