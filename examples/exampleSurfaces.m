clear all; close all; clc;

load surfaceData1200.mat 
load surfMesh.mat

X = xMesh; 
Y = yMesh; 
Z = surfaceData1200;

%% --- Surface plot at t = 1200---
figure();                     % create figure first
h = surf(X, Y, Z);             % plot surface
set(h, 'EdgeColor', 'none');   % hide black grid lines
shading interp;                % smooth shading
colormap parula;               % color map
colorbar;                      % color bar
zlim([-2 2]);                  % adjust Z axis
xlabel('Y'); ylabel('X'); zlabel('Z');
title('Surface Plot of surfaceData1200');

%% --- Gaussian blob ---
A = 50;   % 15 is good in ray
a = 100;   
G = A * exp(- (X.^2 + Y.^2) / a^2);

figure();                      % new figure for blob
h = surf(X, Y, G);
set(h, 'EdgeColor', 'none');
shading interp;
colormap parula;
colorbar;
xlabel('Y'); ylabel('X'); zlabel('Z');
title('Gaussian blob');

%% --- 45 deg surface ---
H = X;

figure();                      % new figure for blob
h = surf(X, Y, H);
set(h, 'EdgeColor', 'none');
shading interp;
colormap parula;
colorbar;
xlabel('Y'); ylabel('X'); zlabel('Z');
title('Gaussian blob');

%% --- Pyramid ---
% Define pyramid shape (45° slopes)
% Each side rises linearly toward the center
Z = -max(abs(X), abs(Y));          % inverted pyramid
Z = Z - min(Z(:));                 % shift so bottom is at z = 0
%Z = max(abs(X), abs(Y));  % upright pyramid (peak)

% Plot
figure;
h = surf(X, Y, Z);
set(h, 'EdgeColor', 'none');
shading interp;
colormap parula;
colorbar;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Pyramid Surface (45° Slopes)');
