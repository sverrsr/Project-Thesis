%%
clear all; close all; clc;

%%
load surfaceData1200.mat 
load surfMesh.mat
%% Plots only the surface at s = 1200
X = xMesh; Y = yMesh; Z = surfaceData1200;

fprintf('X: %dx%d\n', size(X,1), size(X,2));
fprintf('Y: %dx%d\n', size(Y,1), size(Y,2));
fprintf('Z: %dx%d\n', size(Z,1), size(Z,2));


% Read data reads data as (Y, X, Z)

h = surf(X, Y, Z);

set(h, 'EdgeColor','none')   % hide black grid lines
shading interp;      % optional – smooths the surface
colormap parula;        % optional – sets color map
colorbar;            % optional – adds color bar
zlim([-2 2]); %Adjust Z reasonably
xlabel('Y');
ylabel('X');
zlabel('Z');
title('Surface Plot of surfaceData1200');

%% Plotting normals
[Nx, Ny, Nz] = surfnorm(X, Y, Z);
h_norm = surfnorm(Y, X, Z);

figure(2);surfnorm(Y, X, Z);

%% Convert between meshgrid and ngrid formats
% https://se.mathworks.com/help/matlab/ref/ndgrid.html
% MESHGRID: X, Y
% NGRID: Xn, Yn

% Meshgrid -> ngrid
% Needs to be sorted to get the same surface
[xa, ix] = sort(X(1,:));
[ya, iy] = sort(Y(:,1));

[Xn, Yn] = ndgrid(xa, ya);

%% Interpolating and Compare surface normals

% Evaluate Z on Ngrid
Za = Z(iy, ix);                % reorder Z to match sorted axes

% Interpolant built on ngrid
F = griddedInterpolant({xa, ya}, Za.', 'linear', 'none');  % F(X,Y)

% Evaluate Interpolated surface
Zi = F(X', Y');  % should equal Zo (check)
Zi = Zi';

% Orginal surform
[Nx0, Ny0, Nz0] = surfnorm(X, Y, Za);   % "original" normals

%Interpolated surfnorm
[NxI, NyI, NzI] = surfnorm(X, Y, Zi);   % "interpolated" normals

% Align possible global sign flip
dotp = Nx0.*NxI + Ny0.*NyI + Nz0.*NzI;
if median(dotp,'all') < 0
    NxI = -NxI; NyI = -NyI; NzI = -NzI;
end

% Differences (should be ~machine eps)
maxerr = [max(abs(NxI(:)-Nx0(:))), ...
          max(abs(NyI(:)-Ny0(:))), ...
          max(abs(NzI(:)-Nz0(:)))]

% 5) (Optional) quick check plot
figure; surf(X, Y, Zi);
shading interp; colormap parula; colorbar;
title('Interpolated surface');
zlim([-2 2]); %Adjust Z reasonably
xlabel('Y');
ylabel('X');
zlabel('Z');


%% Error
% Assume you've already sorted to xa, ya and reordered Za
[xa, ix] = sort(X(1,:));                % Sort the first row of X and get sorted indices
[ya, iy] = sort(Y(:,1));                % Sort the first column of Y and get sorted indices
Za = Z(iy, ix);                          % Reorder Z based on sorted indices, resulting in Ny×Nx
F  = griddedInterpolant({xa, ya}, Za.', 'linear', 'none');  % Create a gridded interpolant for linear interpolation
[Xn, Yn] = ndgrid(xa, ya);               % Create a grid of points for interpolation (Nx×Ny)
Zhat = F(Xn, Yn);                        % Evaluate the interpolant at the grid points (Nx×Ny)
E = Zhat - Za.';                         % Calculate the error between interpolated and actual values
rmse = sqrt(mean(E(:).^2));              % Compute the root mean square error
mae  = mean(abs(E(:)));                  % Compute the mean absolute error
maxe = max(abs(E(:)));                   % Compute the maximum absolute error
fprintf('RMSE=%.3e, MAE=%.3e, Max=%.3e\n', rmse, mae, maxe);  % Display the error metrics

%% At the nearest grid node. Don't use this. Just to check
x0 = 0; y0 = 0;                          % your point

[~, ix0] = min(abs(xa - x0));
[~, iy0] = min(abs(ya - y0));
n = [Nx(iy0,ix0), Ny(iy0,ix0), Nz(iy0,ix0)];  % already unit length


%% Single points. A. At an arbitrary point (interpolate, then renormalize)


Fx = griddedInterpolant({xa, ya}, Nx.', 'linear', 'none');
Fy = griddedInterpolant({xa, ya}, Ny.', 'linear', 'none');
Fz = griddedInterpolant({xa, ya}, Nz.', 'linear', 'none');

x0 = 0.12; y0 = -0.07;
n = [Fx(x0,y0); Fy(x0,y0); Fz(x0,y0)];
n_a = n ./ norm(n);                        % re-unitize after interpolation

%% New surface



