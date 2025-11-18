close all

%% Configuration
% List the images you want to analyse. The script works with a single image or
% a list of images (png/jpg/mat files). Each image is processed
% independently.
imageFiles = { ...
    'smoothed_image2.png'  % Example input
    % 'screen_1024bins_0002.mat'
};

% Set threshold percentiles between 0 (black) and 255 (white)
warmPercentile = 97;  % warm: top 3%
coldPercentile = 3;   % cold: bottom 3%
areasToKeep    = 20;  % number of largest regions to keep

%% Load all requested images as grayscale double matrices
numImages = numel(imageFiles);
loadedImages = cell(numImages, 1);
labels = cell(numImages, 1);
for i = 1:numImages
    [loadedImages{i}, labels{i}] = loadSquareImage(imageFiles{i});
end

%% Analyse each image independently
for imgIdx = 1:numImages
    X = loadedImages{imgIdx};
    label = labels{imgIdx};

    vals = X(:);
    vals = vals(~isnan(vals));

    thresh    = prctile(vals, warmPercentile);
    invthresh = prctile(vals, coldPercentile);

    %% Warm areas
    Xw = X > thresh;
    CCw = bwconncomp(Xw);
    areasWarm = cell2mat(struct2cell(regionprops(CCw, "Area")));

    [sortedWarm, orderWarm] = sort(areasWarm, 'descend');
    bigWarmIdx = orderWarm(1:min(areasToKeep, numel(orderWarm)));
    smallWarmIdx = setdiff(orderWarm, bigWarmIdx, 'stable');

    XWarmLargest = cc2bw(CCw, ObjectsToKeep=bigWarmIdx);
    XWarmSmall   = cc2bw(CCw, ObjectsToKeep=smallWarmIdx);

    plotAreaSizes(sortedWarm, label, sprintf('Warm areas (>%d percentile)', warmPercentile));

    warmFig = figure('Name', sprintf('Warm areas - %s', label));
    warmFig.Position = [133 549 1200 700];
    t = tiledlayout(warmFig, 2, 2, 'Padding','compact','TileSpacing','compact');
    nexttile(t), plotSquareImage(X, sprintf('Original (%s)', label));
    nexttile(t), plotSquareImage(Xw, sprintf('Mask > %0.2f', thresh));
    nexttile(t), plotSquareImage(XWarmLargest, sprintf('Only %d largest warm areas', numel(bigWarmIdx)));
    nexttile(t), plotSquareImage(XWarmSmall, sprintf('Without %d largest warm areas', numel(bigWarmIdx)));

    %% Cold areas
    Xb = X < invthresh;
    CCb = bwconncomp(Xb);
    areasCold = cell2mat(struct2cell(regionprops(CCb, "Area")));

    [sortedCold, orderCold] = sort(areasCold, 'descend');
    bigColdIdx = orderCold(1:min(areasToKeep, numel(orderCold)));

    XColdLargest = cc2bw(CCb, ObjectsToKeep=bigColdIdx);

    plotAreaSizes(sortedCold, label, sprintf('Cold areas (<%d percentile)', coldPercentile));

    coldFig = figure('Name', sprintf('Cold areas - %s', label));
    coldFig.Position = [183 349 1200 600];
    t = tiledlayout(coldFig, 1, 3, 'Padding','compact','TileSpacing','compact');
    nexttile(t), plotSquareImage(X, sprintf('Original (%s)', label));
    nexttile(t), plotSquareImage(1 - Xb, sprintf('Mask < %0.2f', invthresh));
    nexttile(t), plotSquareImage(1 - XColdLargest, sprintf('Only %d largest cold areas', numel(bigColdIdx)));
end

%% Helper functions
function [img, label] = loadSquareImage(source)
%LOADSQUAREIMAGE Load a file as a grayscale double image.
%   Supports image formats recognised by IMREAD as well as MAT files that
%   contain a struct named `screen` with a field `image`.

    label = source;
    [~, ~, ext] = fileparts(source);
    ext = lower(ext);
    switch ext
        case '.mat'
            data = load(source);
            if isfield(data, 'screen') && isfield(data.screen, 'image')
                img = double(data.screen.image);
                label = sprintf('%s (screen.image)', source);
            else
                error('MAT file %s must contain a struct named screen.image.', source);
            end
        otherwise
            img = imread(source);
            if ndims(img) == 3
                img = rgb2gray(img);
            end
            img = double(img);
    end
end

function plotSquareImage(img, titleText)
%PLOTSQUAREIMAGE Display an image with square aspect ratio.
    imagesc(img);
    axis image off;
    title(titleText, 'Interpreter','none');
    colormap gray;
end

function plotAreaSizes(areaSizes, label, plotTitle)
%PLOTAREASIZES Plot the distribution of area sizes on a semilog graph.
    if isempty(areaSizes)
        warning('No connected components found for %s (%s).', label, plotTitle);
        return;
    end
    figure('Name', sprintf('Area sizes - %s', label));
    semilogy(1:length(areaSizes), areaSizes, 'LineWidth', 1.5);
    xlabel('Area index');
    ylabel('Area (pixels)');
    title(sprintf('%s - %s', plotTitle, label), 'Interpreter','none');
    grid on;
end
