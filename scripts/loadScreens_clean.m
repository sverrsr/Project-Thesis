close all; clear;

inDir  = 'C:\Users\sverr\Documents\NTNU\Prosjekt\Project-Thesis\tenSampledSurfaces_bin';
outDir = 'C:\Users\sverr\Documents\NTNU\Prosjekt\Project-Thesis\tenSampledSurfaces_bin_screens';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

files = dir(fullfile(inDir, 'screen_*.mat'));
[~, idx] = sort({files.name});
files = files(idx);

figure;
colormap gray;

for k = 1:numel(files)
    % Load image
    data = load(fullfile(inDir, files(k).name));
    img  = data.screen.image;

    % Optional blur
    %img  = imgaussfilt(img, 1);

    % Log scale
    %img = log10(img + eps);

    % Show with automatic colorbar per frame
    imagesc(img);
    axis image;
    set(gca, 'YDir', 'normal');
    colorbar;
    title(sprintf('Frame %d / %d', k, numel(files)));
    drawnow;

    % Save PNG (needs [0,1] range)
    img_norm = mat2gray(img);   % auto scale this frame to [0,1]

    [~, base, ~] = fileparts(files(k).name);
    pngName = sprintf('%s_log_%03d.png', base, k);
    imwrite(img_norm, fullfile(outDir, pngName));
end
