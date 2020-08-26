clear all; clc; 

if ~exist('illdata', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/load']);
end

if ~exist('mfit', 'file')
    wd = fileparts(mfilename('fullpath'));
    addpath([wd '/mfit4']);
end

% Load a previous fit if this exist.
if exist('fits.mat', 'file')
    load('fits.mat');
    s0 = setindex;
    f0 = fits;
else
    s0 = cell(1, 7);
    f0 = [];
end

% Loads IN20 data (first run is 91198, last run undetermined)
[x, y, mn, setindex] = flexmerge(91198:99999, 'in20_data', 3e5);

% Only keep energy scans
id_en = cellfun(@(x)strcmp(x, 'EN'), setindex{1});
x = x(id_en);
y = y(id_en);
mn = mn(id_en);
for ii = 1:7
    setindex{ii} = setindex{ii}(id_en);
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
figure;
for ii = 1:length(x)
    clf; hold all;
    title(sprintf('%d - %s', ii, setindex{3}{ii}));
    errorbar(x{ii}, y{ii}.*3e5./mn{ii}, sqrt(y{ii}).*3e5./mn{ii}, 'o');
    fitit = 0;
    if(~isempty(fits{ii}))
        %mx=linspace(min(x{ii}),max(x{ii}),200); plot(mx,ngauss(mx,fits{ii}{1}),'-');
        mx=linspace(0,70,200); plot(mx,ngauss(mx,fits{ii}{1}),'-');
    else
        fitit = 0;
    end
    [mx, my, mb] = ginput(1); 
    if(mb==3 || fitit)
        tomfit(x{ii}, y{ii}.*3e5./mn{ii}, sqrt(y{ii}).*3e5./mn{ii}); 
        mf_rbbox('zoomreset'); 
        title(setindex{3}{ii});
        uiwait(findobj('Tag', 'mf_DataWindow')); 
        [xx, yy, ee, selected, fit, pc, dc, fixed] = fromfit;
        fits{ii} = {pc dc};
    end
end

%save('fits.mat', 'setindex', 'fits');
