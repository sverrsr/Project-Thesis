clear all; clc; %close all;


% load processed_surfElev_333.34.mat %02
% load processed_surfElev_416.68.mat %03
% load processed_surfElev_499.96.mat %04
% load processed_surfElev_583.30.mat %05
% load processed_surfElev_666.64.mat %06
% load processed_surfElev_749.98.mat %07
% load processed_surfElev_833.26.mat %08
% load processed_surfElev_916.60.mat %09
% load processed_surfElev_999.94.mat %10
 load processed_surfElev_1000.00.mat %01


 % Image to compare with
img = double(imread('screen_1024bins_0001_largest_warm_areas.jpg'));

Z = surfElev;
clear surfElev;

nx = 256;
ny = 256;
nt = 12500;

dx = 2 * pi / nx;
dy = 2 * pi / ny;

dt = 0.06;

lx = 2 * pi;
ly = 2 * pi;

nu = 1 / 2500;
overflatespenning = 0;
g = 10;

% Create a new mesh grid based on the specified dimensions
% Create a new mesh grid based on the specified dimensions
% Mesh from -pi to pi
[X, Y] = meshgrid( ...
    single(linspace(-pi, pi, nx)), ...
    single(linspace(-pi, pi, ny)) ...
);


% --- 2D height map ---
figure;
imagesc(X(1,:), Y(:,1), Z);
axis image;
%set(gca, 'YDir', 'normal');
colormap(turbo);
colorbar;
title('Surface elevation (2D height map)');
xlabel('X');
ylabel('Y');

% Nice pi ticks
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})


% Curvature
[K,H,Pmax,Pmin] = curvature(X,Y,Z); %[K, H, Pmax, Pmin] = surfature(X, Y, Z);

% Find correlation of H and surface elevation
R_H = corr(H(:), Z(:), 'rows', 'complete');

unique(img(:))

BW_resized = imresize(img, [size(H,1), size(H,2)], 'nearest');

figure;
imagesc(BW_resized); axis image off; colormap gray;
title('Resized image');

%%
figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(Z); axis image; colormap turbo; colorbar;
title('Z (høyde)');

nexttile;
imagesc(H); axis image; colormap turbo; colorbar;
title('H (kurvatur)');

%% Correlation between image and surface elevation
M = BW_resized;      % mask after resize
Z0 = Z(:);

% no flip
c_noflip = corr(M(:), Z0, 'rows','complete');

% flip up-down
M_ud = flipud(M);
c_flipud = corr(M_ud(:), Z0, 'rows','complete');

% flip left-right
M_lr = fliplr(M);
c_fliplr = corr(M_lr(:), Z0, 'rows','complete');

% rotate 90
M_r90 = rot90(M);
c_rot90 = corr(M_r90(:), Z0, 'rows','complete');

[c_noflip c_flipud c_fliplr c_rot90]

%% Correlation between image and curvature
M = BW_resized;      % mask after resize
H0 = H(:);

% no flip
cH_noflip = corr(M(:), H0, 'rows','complete');

% flip up-down
M_ud = flipud(M);
cH_flipud = corr(M_ud(:), H0, 'rows','complete');

% flip left-right
M_lr = fliplr(M);
cH_fliplr = corr(M_lr(:), H0, 'rows','complete');

% rotate 90
M_r90 = rot90(M);
cH_rot90 = corr(M_r90(:), H0, 'rows','complete');

[cH_noflip cH_flipud cH_fliplr cH_rot90]

%%
mask = logical(BW_resized);          % 1 = scar, 0 = background

H_scar = H(mask);           % curvature values on scars
H_bg   = H(~mask);          % curvature values off scars

figure;
histogram(H_scar); hold on;
histogram(H_bg);
legend('scar','background');
title('Mean curvature distribution');

figure;
boxplot([H_scar; H_bg], [ones(size(H_scar)); 2*ones(size(H_bg))])
%set(gca,'XTickLabel',{'scar','background'})
ylabel('Mean curvature');

mean_H_scar = mean(H_scar);
mean_H_bg   = mean(H_bg);
std_H_scar  = std(H_scar);
std_H_bg    = std(H_bg);

[~,p] = ttest2(H_scar, H_bg);    % two-sample t-test

scar_indicator = double(mask(:));   % 1 = scar, 0 = no scar
H_flat = H(:);

R = corr(scar_indicator, H_flat);   % Pearson correlation

figure;
imagesc(H); axis image; colorbar; hold on;
contour(mask, [0.5 0.5], 'w', 'LineWidth', 1);   % white contour of scars
title('Mean curvature with scar locations');
