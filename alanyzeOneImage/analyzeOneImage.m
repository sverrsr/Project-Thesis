%% Single Screen Smoothing Analysis (black & white)
% Works on ONE screen_*.mat file at a time.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------------------------------
% Folder and file
%inDir   = 'C:\Users\sverr\Documents\NTNU\Prosjekt\Project-Thesis\tenSampledSurfaces_bin';
fileName = 'screen_1024bins_0002.mat';   % <--- change to any file you want

% Pre-processing options
useLogScale     = true;      % log(img+1) to reveal weak structures
useContrastStretch = true;   % percentile-based contrast stretch
lowPct          = 1;         % lower percentile for contrast (0–10)
highPct         = 99;        % upper percentile for contrast (90–100)

useGaussianBlur = true;      % apply Gaussian smoothing
gaussSigma      = 3;       % blur strength (0 = off)

useHistEq       = false;     % adaptive histogram equalization (CLAHE)

% Colormap / display
useManualCLim   = false;     % set true to force fixed caxis
manualCLim      = [0 4];     % only used if useManualCLim = true
% ------------------------------------------------------------------------


%% Load screen and extract image
data   = load(fileName) %load(fullfile(inDir, fileName));
screen = data.screen;              % Screen object
imgRaw = double(screen.image);     % make sure it's double

fprintf('Loaded %s, size = [%d, %d]\n', fileName, size(imgRaw,1), size(imgRaw,2));

%% Step 1: basic contrast info
rawMin = min(imgRaw(:));
rawMax = max(imgRaw(:));
fprintf('Raw intensity range: [%.3e, %.3e]\n', rawMin, rawMax);

%% Step 2: contrast stretching (percentile based)
imgCS = imgRaw;

if useContrastStretch
    pLow  = prctile(imgRaw(:), lowPct);
    pHigh = prctile(imgRaw(:), highPct);
    if pHigh <= pLow
        warning('Percentile range collapsed; disabling contrast stretch.');
    else
        imgCS = imgRaw;
        imgCS(imgCS < pLow)  = pLow;
        imgCS(imgCS > pHigh) = pHigh;
        imgCS = (imgCS - pLow) ./ (pHigh - pLow);  % scale to [0,1]
    end
else
    % normalize to [0,1] using min/max
    imgCS = (imgRaw - rawMin) / max(rawMax - rawMin, eps);
end

%% Step 3: log scaling
if useLogScale
    imgLog = log1p(imgCS);                % log(1+x), x in [0,1]
    imgLog = imgLog / max(imgLog(:));     % normalize to [0,1]
else
    imgLog = imgCS;
end

%% Step 4: Gaussian smoothing
if useGaussianBlur && gaussSigma > 0
    imgSmooth = imgaussfilt(imgLog, gaussSigma);
else
    imgSmooth = imgLog;
end

%% Step 5: adaptive histogram equalization (optional)
if useHistEq
    imgFinal = adapthisteq(imgSmooth);
else
    imgFinal = imgSmooth;
end

%% Display: 2×3 B/W overview
figure;
tiledlayout(2,3, "TileSpacing","compact", "Padding","compact");

% helper for consistent display
function showimg(im, ttl, useManualCLim)
        imagesc(im);
        axis image off;
        colormap(gray);
        if useManualCLim
            caxis(manualCLim);
        end
        title(ttl, 'Interpreter','none');
    end

nexttile; showimg(imgRaw,   'Raw', useManualCLim);
nexttile; showimg(imgCS,    'Contrast stretched', useManualCLim);
nexttile; showimg(imgLog,   'Log-scaled', useManualCLim);

nexttile; showimg(imgSmooth,'Smoothed', useManualCLim);
nexttile; showimg(imgFinal, 'Final (after all steps)', useManualCLim);

% Histogram of final image
nexttile;
histogram(imgFinal(:), 256);
xlabel('Intensity'); ylabel('Count');
title('Histogram (final)');
grid on;

%% Done
disp('Single-image analysis complete.');
