% Load the data
u = h5read('wave.h5', '/u');   % shape: [Ny, Nx, Nt]

% Define spatial grid (same as before)
Ny = 501; Nx = 408; Nt = 1201;

load surfaceData1200.mat 
load surfMesh.mat


X = xMesh; Y = yMesh; %Z = surfaceData1200;

A = 50;

% Set up figure
figure
hSurf = surf(X, Y, u(:,:,1));
shading interp
colormap parula
axis tight
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) -A A])
xlabel('x'); ylabel('y'); zlabel('Amplitude');

% Animate
for k = 1:Nt
    set(hSurf, 'ZData', u(:,:,k));
    title(sprintf('Frame %d / %d', k, Nt));
    drawnow
end
