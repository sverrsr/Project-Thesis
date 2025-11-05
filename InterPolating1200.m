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

%% Interpolating and Compare surface normals

% Orginal surform
[Nx0, Ny0, Nz0] = surfnorm(X, Y, Z);   % "original" normals

%Interpolated surfnorm
[NxI, NyI, NzI] = surfnorm(X, Y, Zi);   % "interpolated" normals


%% Finding surface normals 
% This code calculates the gradient of the interpolated surface Zi, 
% normalizes the resulting vectors, 
% and checks for consistency in their direction.

tic
[dZdx, dZdy] = gradient(Zi, xa, ya);          % X spacing first, then Y

NxB = -dZdx;  NyB = -dZdy;  NzB = ones(size(Zi));
L = sqrt(NxB.^2 + NyB.^2 + NzB.^2);
NxB = NxB./L; NyB = NyB./L; NzB = NzB./L;

dpB = Nx0.*NxB + Ny0.*NyB + Nz0.*NzB;

if median(dpB(:)) < 0
    NxB = -NxB; NyB = -NyB; NzB = -NzB;
    dpB = -dpB;                                  % keep angles consistent
end

maxerrB = [max(abs(NxB(:)-Nx0(:))), max(abs(NyB(:)-Ny0(:))), max(abs(NzB(:)-Nz0(:)))];
maxAngB = max(acosd(max(-1,min(1,dpB(:)))));
tB = toc;

%% SINGLE POINT EVALUATION
% It is now possible to evaluate both the height and the surface normal in
% an arbitrary point

[dZdx, dZdy] = gradient(Zi, xa, ya);          % X spacing first, then Y
Fdx = griddedInterpolant({ya, xa}, dZdx, 'linear');   % or 'makima'/'spline'
Fdy = griddedInterpolant({ya, xa}, dZdy, 'linear');

xq = 0.0; yq = 0.0; % Example points

% Evaluate surface height at a single point (xq, yq)
z0 = F(xq, yq);

% --- Evaluate surface normal at a single point (xq, yq)
zx = Fdx(yq, xq);
zy = Fdy(yq, xq);

n0 = [-zx, -zy, 1];  n0 = n0 / norm(n0);

fprintf('NORM(%g, %g) = %.6g\n', xq, yq, n0);
fprintf('Z(%g, %g) = %.6g\n', xq, yq, z0);



%% Report
fprintf('Method B (gradient):   time = %.4f s, max |err| = [%.3g %.3g %.3g], max angle = %.3g deg\n', ...
         tB, maxerrB(1), maxerrB(2), maxerrB(3), maxAngB);

%% Removing end nodes
I = 2:size(Zi,1)-1; J = 2:size(Zi,2)-1;

% angles already from dpB after sign-fix
angB = acosd(max(-1,min(1,dpB)));
maxAng_int = max(angB(I,J), [], 'all');
p95Ang = prctile(angB(I,J), 95);

fprintf('Interior max angle: %.3g deg, 95th pct: %.3g deg\n', maxAng_int, p95Ang);

% component errors
ex = abs(NxB-Nx0); ey = abs(NyB-Ny0); ez = abs(NzB-Nz0);
fprintf('Interior max |err|: [%.3g %.3g %.3g]\n', ...
    max(ex(I,J),[],'all'), max(ey(I,J),[],'all'), max(ez(I,J),[],'all'));

%%
% Your current query normal n0 at (xq,yq) is already computed above.

% 1) Get node normals like MATLAB's surf does
[Xg, Yg] = meshgrid(xa, ya);                 % size matches Zi: rows=ya, cols=xa
[Nsx, Nsy, Nsz] = surfnorm(Xg, Yg, Zi);      % node normals

% 2) Interpolate those node normals to the same query point
FNsx = griddedInterpolant({ya, xa}, Nsx, 'linear');   % 'makima'/'spline' also fine
FNsy = griddedInterpolant({ya, xa}, Nsy, 'linear');
FNsz = griddedInterpolant({ya, xa}, Nsz, 'linear');

ns = [FNsx(yq, xq), FNsy(yq, xq), FNsz(yq, xq)];      % surf-style normal at (xq,yq)

% 3) Normalize and align orientation (so we don't compare a flipped vector)
ns = ns / norm(ns);
if dot(ns, n0) < 0, ns = -ns; end

% 4) Compare: angle and component-wise diff
ang_deg = acosd(max(-1, min(1, dot(n0, ns))));
comp_err = abs(n0 - ns);

fprintf('Angle difference: %.6g deg\n', ang_deg);
fprintf('Component |err|: [%.3g  %.3g  %.3g]\n', comp_err(1), comp_err(2), comp_err(3));







