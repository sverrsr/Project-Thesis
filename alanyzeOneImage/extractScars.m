function analyzeSurfaceScars(pathToMatFile)

fileName = 'screen_1024bins_0005.mat';   % <--- change to any file you want

%% LOAD IMAGE
data = load(fileName);
img = double(data.screen.image);

%% USER PARAMETERS
params.sigmaSmooth = 4;        % Gaussian smoothing
params.sigmaHigh   = 1.5;      % DoG high
params.sigmaLow    = 8;        % DoG low
params.sigmaHess   = 2;        % Hessian scale
params.thr         = 0.15;     % threshold for ridge image
params.minObjSize  = 50;       % remove small blobs

%% STEP 1 — SMOOTH
imgSmooth = imgaussfilt(img, params.sigmaSmooth);
%% SAVE SMOOTHED IMAGE
imwrite(uint8(imgSmooth), 'smoothed_image2.png');

%% STEP 2 — BANDPASS
imgHigh = imgaussfilt(imgSmooth, params.sigmaHigh);
imgLow  = imgaussfilt(imgSmooth, params.sigmaLow);
imgBand = imgHigh - imgLow;

%% STEP 3 — HESSIAN-BASED RIDGE DETECTION
ridge = ridgeHessian(imgBand, params.sigmaHess);

%% STEP 4 — THRESHOLD + MORPHOLOGY
bw = ridge > params.thr;
bw = bwareaopen(bw, params.minObjSize);
bw = imclose(bw, strel('disk',3));
bw = imopen(bw, strel('disk',1));

%% STEP 5 — OPTIONAL SKELETON
skel = bwskel(bw,'MinBranchLength',20);

%% PLOT ALL STEPS (same layout you liked)
figure('Position',[100 100 1800 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; imagesc(img); colormap gray; axis image off;
title('Raw');

nexttile; imagesc(imgSmooth); colormap gray; axis image off;
title('Smoothed');

nexttile; imagesc(imgBand); colormap gray; axis image off;
title('Bandpass (DoG)');

nexttile; imagesc(ridge); colormap gray; axis image off;
title('Ridge Strength (Hessian)');

nexttile; imagesc(bw); colormap gray; axis image off;
title('Binary scars');

nexttile;
histogram(ridge(:),200,'FaceColor',[0.3 0.3 0.3]);
xlabel('Intensity'); ylabel('Count'); grid on;
title('Histogram (ridge response)');

%% Skeleton figure
figure;
imagesc(skel); axis image off; colormap gray;
title('Skeleton (curve tracing)');

end

function R = ridgeHessian(I, sigma)
% Computes ridge response using Hessian eigenvalues (2D)

I = double(I);
g = fspecial('gaussian', ceil(6*sigma), sigma);
Ig = imfilter(I, g, 'replicate');

% Derivatives
[Ix, Iy] = gradient(Ig);
[Ixx, Ixy] = gradient(Ix);
[~, Iyy] = gradient(Iy);

% Hessian eigenvalues
tmp = sqrt((Ixx - Iyy).^2 + 4*Ixy.^2);
lambda1 = 0.5 * (Ixx + Iyy + tmp);
lambda2 = 0.5 * (Ixx + Iyy - tmp);

% Ridge strength = negative eigenvalue magnitude
ridge = max(-lambda1, -lambda2);

% Normalize
R = ridge / max(ridge(:) + eps);
end
