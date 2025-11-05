%%
clear all; close all; clc;
load surfaceData1200.mat 
load surfMesh.mat
%% Plots only the surface at s = 1200
X = xMesh; Y = yMesh; Z = surfaceData1200;
clear xMesh yMesh surfaceData1200

%% Convert between meshgrid and ngrid formats
% Gridded interpolant uses NGRID format, so it's necessary to convert 
% https://se.mathworks.com/help/matlab/ref/ndgrid.html
% MESHGRID: X, Y
% NGRID: Xn, Yn

% Meshgrid -> ngrid
% Needs to be sorted to get the same surface
[xa, ix] = sort(X(1,:));
[ya, iy] = sort(Y(:,1));

% Unsorted NGRID is made like this, but is not used in this code
% [Xn, Yn] = ndgrid(xa, ya);

%% Interpolating and Compare surface normals
% To evaluate the surface at all points it is necessary to interpolate the
% surfae

% Evaluate Z on Ngrid
Za = Z(iy, ix);

% Interpolant built on ngrid
F = griddedInterpolant({xa, ya}, Za.', 'linear', 'none');  % F(X,Y)
clear Za

% Zi is teh interpolated surface
% Evaluate Interpolated surface
% Converting to Ngrid yield better performance. Keep '
Zi = F(X', Y')';  isequal(Zi, Z); % should equal Z (check).

%% SINGLE POINT EVALUATION
% It is now possible to evaluate both the height and the surface normal in
% an arbitrary point

% griddedInterpolant uses NDGRID order (rows→ya, cols→xa), 
% so vectors are {ya, xa} to match dZdx/dZdy layout
% Evaluate in that order: Fdx(yq, xq) and Fdy(yq, xq), 
% which matches {ya, xa}.
[dZdx, dZdy] = gradient(Zi, xa, ya);          % X spacing first, then Y
Fdx = griddedInterpolant({ya, xa}, dZdx, 'linear','none');   % or 'makima'/'spline'
Fdy = griddedInterpolant({ya, xa}, dZdy, 'linear','none');

%% To be evaluated inside the lens function
xq = 0.0; yq = 0.0; % Example points

% Evaluate surface height at a single point (xq, yq)
z0 = F(xq, yq);

% --- Evaluate surface normal at a single point (xq, yq)
zx = Fdx(yq, xq);
zy = Fdy(yq, xq);

n0 = [-zx, -zy, 1];  n0 = n0 / norm(n0);

fprintf('NORM(%g, %g) = %.6g\n', xq, yq, n0);
fprintf('Z(%g, %g) = %.6g\n', xq, yq, z0);

