function out = fvoigt(xdat,pars)
% Calculates a pseudo-voigt function with given parameters at given ordinate.
%
% Syntax:  out  = fvoigt(xdat,pars)
%
% Inputs:  xdat - vector - is the ordinate at which to calculate
%          pars - vector - [centre fwhm area lfrac] - lfrac is lorentzian fraction
%
% Outputs: out  - vector - is the abscisa which is calculated.

% Reference: pvoigt.m file from the fitting routines for ID20 by SPC, SBW.

% Makes equation looks better
%c = pars(1);
%w = pars(2);
%a = pars(3);
%f = pars(4);
x = xdat(:);

%out = (a/w/ (f*pi/2 + (1-f)*sqrt(pi/4/log(2)) )) ...
%      .* ( f./(1 + 4*((x-c)/w).^2) + (1-f)*exp(-4*log(2)*((x-c)/w).^2) ); 

np = numel(pars);
if mod(np, 4) ~= 0
    error('Input must be N sets of parameters [centre fwhm area frac] where N is number of peaks');
end
np = np / 4;

out = x * 0;

for n = 1:np
    c = pars((n-1)*4 + 1);
    w = pars((n-1)*4 + 2);
    a = pars((n-1)*4 + 3);
    f = pars((n-1)*4 + 4);
    
    sigma = w / (2*sqrt(2*log(2)));
    ag = (2/w) * sqrt(log(2)/pi);
    bg = 4 * log(2) / w^2;
    
    g = ag .* exp(-bg * (x-c).^2);
    l = ((1/pi) * (w/2)) ./ ((x-c).^2 + (w/2)^2);
    
    out = out + a .* (f.*g + (1-f).*l);
end
