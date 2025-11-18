% -----------------------------------------------------------------------------
% Creates an MP4 or AVI animation from .mat image files.
% The Images must be screen objects. See OPTOMETRIKA Toolbox
% Loads each frame, applies Gaussian blur, updates the plot, 
% and writes frames to a video file with fixed color scale and colormap.
%
% Options:
%   outDir       – output folder for the video
%   vidName      – video file name
%   v.FrameRate  – playback speed (fps)
%   imgaussfilt(img, 2) – blur strength (higher = smoother)
%   caxis([0 0.05]) – color range for all frames
%
% Output:
%    figure
%    mp4 or avi file 
% -----------------------------------------------------------------------------

close all; clear all;

%GaussFiltVal = 0.5; %0.5 er bra med lim
%caxis([0, 4]);

% Create output folder
outDir = 'C:\Users\sverr\Documents\NTNU\Prosjekt\Project-Thesis\tenSampledSurfaces_bin_screens';
inDir = 'C:\Users\sverr\Documents\NTNU\Prosjekt\Project-Thesis\tenSampledSurfaces_bin';
vidName = 'DNS_SCREENS_150k_dist3pi_binAnalysis.mp4';

% outDir = '\\sambaad.stud.ntnu.no\sverrsr\Documents\DNS_SCREENS_150k_dist3pi';
% inDir = '\\sambaad.stud.ntnu.no\sverrsr\Documents\DNS_SCREENS_150k_dist3pi';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% AVI file:
%vidName = 'blurred_animation.avi'
%vidPath = fullfile(outDir, vidName);
%v = VideoWriter(vidPath);

% MP4 file:
vidPath = fullfile(outDir, vidName);
v = VideoWriter(vidPath, 'MPEG-4');


files = dir(fullfile(inDir, 'screen_*')); %.mat 
files = {files.name};
files = sort(files);

v.FrameRate = 20; % fps
open(v);

% --- Load first frame to initialize ---
data = load(fullfile(inDir, files{1}));
img = data.screen.image;

%img = imgaussfilt(img, GaussFiltVal);

% --- Set up figure once ---
figure;
hImg = imagesc(img);
axis image;
set(gca, 'YDir', 'normal');
colormap("gray"); %flipud(sky)
colorbar;

% Compute global limits across all frames
minVal = inf; maxVal = -inf;

for k = 1:numel(files)
    data = load(fullfile(inDir, files{k}));
    img = data.screen.image;
    minVal = min(minVal, min(img(:)));
    maxVal = max(maxVal, max(img(:)));
end

fprintf('Global intensity range: [%.3e, %.3e]\n', minVal, maxVal);

% adjust range to match your data
% Must be adjusted to Gauss filtering. More gauss, lower axis

%caxis([minVal, maxVal]);
colorbar; hold on;
caxis([0, 4]);

%%
% --- Main animation loop ---
for k = 1:numel(files)
    data = load(fullfile(inDir, files{k}));

    img = data.screen.image;

    %img = imgaussfilt(img, GaussFiltVal);
    set(hImg, 'CData', img);  % update only the image data
    view(180,90)
    title(sprintf('Frame %d / %d', k, numel(files)));
    drawnow;
    hold on;

    originalFilename = files{k};
    frame = getframe(gcf);
    %writeVideo(v, frame);
    % Save each screen as a .png file

    % Optional blur
    
    img = imgaussfilt(img, 1.5);
    
    imwrite(img, fullfile(outDir, sprintf('screen_filtered_1.5_%03d_%s.png', k, originalFilename)));

end

close(v);
fprintf('Saved animation in: %s\n', vidPath);
    
    