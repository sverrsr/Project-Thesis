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

GaussFiltVal = 0.5; %0.5 er bra med lim
caxis([1, 2]);

% Create output folder
outDir = 'D:\';
inDir = 'D:\DNS_SCREENS_150k_dist3pi_fullSIM';
vidName = 'DNS_SCREENS_150k_dist3pi.mp4';

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

if isfield(data, 'screen_image')
    img = data.screen_image;
elseif isfield(data, 'screen')
    img = data.screen.image;
else
    fns = fieldnames(data);
    img = data.(fns{1});
end

img = imgaussfilt(img, GaussFiltVal);

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
    if isfield(data, 'screen_image')
        img = data.screen_image;
    elseif isfield(data, 'screen')
        img = data.screen.image;
    else
        fns = fieldnames(data);
        img = data.(fns{1});
    end
    minVal = min(minVal, min(img(:)));
    maxVal = max(maxVal, max(img(:)));
end
fprintf('Global intensity range: [%.3e, %.3e]\n', minVal, maxVal);

% adjust range to match your data
% Must be adjusted to Gauss filtering. More gauss, lower axis
caxis([0.5, 4]);

%%
% --- Main animation loop ---
for k = 1:numel(files)
    data = load(fullfile(inDir, files{k}));

    if isfield(data, 'screen_image')
        img = data.screen_image;
    elseif isfield(data, 'screen')
        img = data.screen.image;
    else
        fns = fieldnames(data);
        img = data.(fns{1});
    end

    img = imgaussfilt(img, GaussFiltVal);
    set(hImg, 'CData', img);  % update only the image data
    view(180,90)
    title(sprintf('Frame %d / %d', k, numel(files)));
    drawnow;

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Saved animation in: %s\n', vidPath);
