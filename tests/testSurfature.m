clear all; clc; close all;

% load processed_surfElev_499.96.mat
% load processed_surfElev_333.34.mat %02
% load processed_surfElev_416.68.mat
% load processed_surfElev_499.96.mat
% load processed_surfElev_583.30.mat
% load processed_surfElev_666.64.mat
% load processed_surfElev_749.98.mat
% load processed_surfElev_833.26.mat
% load processed_surfElev_916.60.mat
% load processed_surfElev_999.94.mat 
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
set(gca, 'YDir', 'normal');
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



%


unique(img(:))

BW_resized = imresize(img, [size(H,1), size(H,2)], 'nearest');

figure;
imagesc(BW_resized); axis image off; colormap gray;
title('Resized image');


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
set(gca,'XTickLabel',{'scar','background'})
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
