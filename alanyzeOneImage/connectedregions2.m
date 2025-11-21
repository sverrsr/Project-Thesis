clear all; clc; close all;

fileName = 'screen_1024bins_0001.mat';   % <--- change to any file you want

% Create a new folder to save processed files
outputFolder = fullfile('C:\Users\sverr\Documents\NTNU\Prosjekt\Project-Thesis', 'tenSampledSurfaces_smoothed_filtered');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created folder: %s\n', outputFolder);
end

%% LOAD IMAGE
data = load(fileName);
img = double(data.screen.image);

%% USER PARAMETERS
params.sigmaSmooth = 4;        % Gaussian smoothing, 2 seemed ok
params.sigmaHigh   = 1.5;      % DoG high
params.sigmaLow    = 8;        % DoG low
params.sigmaHess   = 2;        % Hessian scale
params.thr         = 0.1;     % threshold for ridge image
params.minObjSize  = 50;       % remove small blobs


%% STEP 1 — SMOOTH
imgSmooth = imgaussfilt(img, params.sigmaSmooth);

%% SHOW SMOOTHED IMAGE
figure;
imagesc(imgSmooth); axis image off; colormap gray;
title('Smoothing Result');

%% STEP 1.5 — INCREASE CONTRAST
imgContrast = imadjust(imgSmooth);


%% SHOW CONTRAST ENHANCED IMAGE
figure;
imagesc(imgContrast); axis image off; colormap gray;
title('Contrast Enhanced Result');

%% SAVE SMOOTHED IMAGE
%imwrite(uint8(imgSmooth), 'screen_1024bins_0002_smoothed.jpg');
[~, baseName, ~] = fileparts(fileName);   % baseName = "screen_1024bins_0002"
outFile = fullfile(outputFolder, baseName + "_smoothed.mat");

% save(outFile, "imgSmooth");
% fprintf('Saved: %s\n', outFile);

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
%figure('Position',[100 100 1800 900]);
figure();
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; imagesc(imgSmooth); colormap gray; axis image off;
title('Smoothed');

nexttile; imagesc(imgContrast); colormap gray; axis image off;
title('Adjusted contrast');

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


%% Connected Areas
X1 = imgContrast; %Setting it to bw can seem to make it better

%%%%%%%%%%%% CONNECTED AREAS ANALYSIS %%%%%%%%%%%%%%
% Provided by Simen Ådnøy Ellingsen

%Set threshhold between 0 (black) and 1 (white)
thresh = 0.7;
invthresh = 0.7;
areastokeep = 50;

%Make binary images
X1w = X1*0 + (X1>thresh);       %For light-area
%X2w = X2*0 + (X2>thresh);

CC1w = bwconncomp(X1w);
%CC2w = bwconncomp(X2w);

areas1w = cell2mat(struct2cell(regionprops(CC1w,"Area")));
%areas2w = cell2mat(struct2cell(regionprops(CC2w,"Area")));

%Lists of sizes of areas
[Asz1,order1] = sort(areas1w,'descend');
%[Asz2,order2] = sort(areas2w,'descend');
n1 = length(order1); 
%n2 = length(order2);


%Separate out the largest areas
bigidx1 = order1(1:areastokeep); 
%bigidx2 = order2(1:areastokeep);
smallidx1 = order1(areastokeep+1:n1); 
%smallidx2 = order2(areastokeep+1:n2); 

X1tBig = cc2bw(CC1w,ObjectsToKeep=bigidx1);
%X2tBig = cc2bw(CC2w,ObjectsToKeep=bigidx2);
X1tSmall = cc2bw(CC1w,ObjectsToKeep=smallidx1);
%X2tSmall = cc2bw(CC2w,ObjectsToKeep=smallidx2);

f=figure;
%set(f,"Position",[ 133         549        1594         730])
colormap gray
 
t=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile, imagesc(X1), axis image off,title('Original');
%nexttile, imagesc(X2), axis image off, title('Original, ');
nexttile, imagesc(X1w), axis image off,title(sprintf('Threshold: %d',thresh));
%nexttile, imagesc(X2w), axis image off, title(sprintf('. Threshold: %d',thresh));
nexttile, imagesc(X1tBig), axis image off, title(sprintf('Threshold: %d, Only %d largest warm areas',thresh,areastokeep));
%nexttile, imagesc(X2tBig), axis image off, title(sprintf('Threshold: %d, Only %d largest warm areas',thresh,areastokeep));
nexttile, imagesc(X1tSmall), axis image off, title(sprintf('Threshold: %d, Without %d largest warm areas',thresh,areastokeep));
%nexttile, imagesc(X2tSmall), axis image off, title(sprintf('Threshold: %d, Without %d largest warm areas',thresh,areastokeep));


% Save the binary image of the largest warm areas
[~, baseName, ~] = fileparts(fileName);
outImg = fullfile(outputFolder, baseName + "_largest_warm_areas.jpg");

imwrite(X1tBig, outImg);
fprintf('Saved: %s\n', outImg);


%%%%%%%%%%%%%% SAME FOR DARK AREAS %%%%%%%%%%%%%% 
X1b = X1*0 + (X1<invthresh);    %For dark-area analysis
%X2b = X2*0 + (X2<invthresh);

CC1b = bwconncomp(X1b);
%CC2b = bwconncomp(X2b);

areas1b = cell2mat(struct2cell(regionprops(CC1b,"Area")));
%areas2b = cell2mat(struct2cell(regionprops(CC2b,"Area")));

%Lists of sizes of areas
[Asz1,order1] = sort(areas1b,'descend');
%[Asz2,order2] = sort(areas2b,'descend');


figure
%semilogy(1:length(Asz1),Asz1,1:length(Asz2),Asz2),
semilogy(1:length(Asz1), Asz1)
%legend('No waves', '');
%ylim([5 max([Asz1(1) Asz2(1)])])
ylim([5 max([Asz1(1)])])
title('Area sizes, hot areas areas');
ylabel('Area (pixels)');
xlabel('Area index');

figure
%semilogy(1:length(Asz1),Asz1,1:length(Asz2),Asz2)
semilogy(1:length(Asz1),Asz1)
%legend('No waves', '');
%ylim([5 max([Asz1(1) Asz2(1)])])
ylim([5 max([Asz1(1)])])
title('Area sizes, cold areas');
ylabel('Area (pixels)');
xlabel('Area index');


%Just the largest areas
n1 = length(order1);
bigidx1 = order1(1:n1); 
%bigidx2 = order2(1:areastokeep);
X1tBig = cc2bw(CC1b,ObjectsToKeep=bigidx1);
%X2tBig = cc2bw(CC2b,ObjectsToKeep=bigidx2);


f=figure;
%set(f,"Position",[ 183         349        1594         620])
colormap gray
 
t=tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
nexttile, imagesc(X1), axis image off,title('Original');
%nexttile, imagesc(X2), axis image off, title('Original, ');
nexttile, imagesc(1-X1b), axis image off,title(sprintf('Threshold: %d',invthresh));
%nexttile, imagesc(1-X2b), axis image off, title(sprintf('. Threshold: %d',invthresh));
nexttile, imagesc(1-X1tBig), axis image off, title(sprintf('Threshold: %d, Only %d largest cold areas',invthresh,areastokeep));
%nexttile, imagesc(1-X2tBig), axis image off, title(sprintf('Threshold: %d, Only %d largest cold areas',invthresh,areastokeep));


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
