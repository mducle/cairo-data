function chi2 = bfo_fit_iter(Jpars, do_plot, save_file_name)
% Runs one iteration of the fit of INS data on Bi2Fe4O9.
% The iteration first fits the single-crystal IN20 data with J_c, J_43 and J'_43 parameters
% Then it evaluates the powder spectra with the other parameters.
%
% To run a fit, use:
%
% x0 = [1.3 3.7 2.9 6.3 24.0];  % Parameters in Beauvois, Simonet et al.
% save('fit_par_tmp.mat', x0);
% x_opt = fminsearch(@bfo_fit_iter, x0([1 5]))
%
% Where we only use fminsearch to fit the J44 and J33 parameters which determine the high
% energy modes seen in the powder data (the other parameters are fitted internally on the
% IN20 data).

if nargin < 3
    save_file_name = 'fit_par_tmp.mat';
end

% If temp save file not present use Beauvois parameters
if exist(save_file_name, 'file')
    x0 = load(save_file_name);
else
    x0 = struct;
    x0.x0 = [Jpars(1) 3.7 2.9 6.3 Jpars(2) 0.021];
end

% Add SIA initial value if it is not there
if numel(x0.x0) == 5;
    x0.x0(6) = 0.021;
end

if nargin < 2
    do_plot = false;
end

if nargin < 1
    Jpars = [1.3 3.7 2.9 6.3 24.0];
else
    Jpars = [Jpars(1) x0.x0(1:3) Jpars(2) x0.x0(4)];
end

if ~exist('bi2fe4o9_spinw', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/matlab']);
end

fit_xtl = fit_in20_spinw(Jpars, [0 1 1 1 0], do_plot);

disp_struct = struct;
disp_struct.x = [Jpars(1) fit_xtl.x(1:3) Jpars(5) fit_xtl.x(4)];

chi2_pow = bfo_powder_iter(disp_struct.x, do_plot);
%chi2_pow = chi2_pow ./ [1 1 10 1 5 5];
chi2_vec = [chi2_pow fit_xtl.redX2];
if ~isfield(x0, 'norm')
    norm = chi2_vec;
else
    norm = x0.norm;
end

disp_struct.chi2_vec = chi2_vec;
disp_struct.chi2_norm = (chi2_vec ./ norm) .* 10;
disp_struct.chi2 = mean(disp_struct.chi2_norm);
disp_struct

chi2 = disp_struct.chi2;

if isfield(x0, 'xvec')
    xvec = [x0.xvec disp_struct];
else
    xvec = [disp_struct];
end
x0 = disp_struct.x;
save(save_file_name, 'x0', 'xvec', 'norm');
