clear all; clc; 

if ~exist('illdata', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/load']);
end

if ~exist('mfit', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/mfit4']);
    addpath([wd '/nllsq']);
end

% Load a previous fit if this exist.
if exist('fits_eiger.mat', 'file')
    load('fits_eiger.mat');
    s0 = setindex;
    f0 = fits;
else
    s0 = cell(1, 7);
    f0 = [];
end

% Loads Eiger data 
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
    y{ist} = (ys{ist}./ms{ist}*1e6) .* (1-exp(-en*(11.604519339/setindex{6}(ist))));
    e{ist} = (sqrt(ys{ist})./ms{ist}*1e6) .* (1-exp(-en*(11.604519339/setindex{6}(ist))));
end

% Parse to see if we've seen this data (and/or fit it) before.
for ii = 1:length(x)
    fits{ii} = [];
    if (ii <= length(s0{3}))
        isk = strcmp(s0{3}{ii}, setindex{3}{ii});
        if (~isk)
            i2 = find(strcmp(setindex{3}{ii}, s0{3}));
            if ~isempty(i2) 
              fits{ii} = f0{i2}; 
            end
        else
            fits{ii} = f0{ii};
        end
    end
end

% Plots the energy scans. If right click, will push data to mfit
hf = figure;
for ii = 3:22;
    set(0, 'CurrentFigure', hf);
    clf; hold all;
    title(sprintf('%d - %s', ii, setindex{3}{ii}));
    errorbar(en, y{ii}, e{ii}, 'o');
    fitit = 0;
    if(~isempty(fits{ii}))
        mx=linspace(0,70,200); plot(mx,ngauss(mx,fits{ii}{1}),'-');
    else
        fitit = 0;
    end
    waitforbuttonpress;
    [mx, my, mb] = ginput(1); 
    if(mb==3 || fitit)
        ynn = y{ii}; enn = e{ii}; xnn = en; idn = isnan(ynn); ynn(idn) = []; enn(idn) = []; xnn(idn) = [];
        tomfit(xnn, ynn, enn); 
        mf_rbbox('zoomreset'); 
        title(setindex{3}{ii});
        uiwait(findobj('Tag', 'mf_DataWindow')); 
        [xx, yy, ee, selected, fit, pc, dc, fixed] = fromfit;
        fits{ii} = {pc dc};
    end
end

%save('fits_eiger.mat', 'setindex', 'fits');
