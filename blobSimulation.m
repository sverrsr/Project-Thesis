% Load the data
u = h5read('wave.h5', '/u');   % shape: [Ny, Nx, Nt]

% Define spatial grid (same as before)
Ny = 501; Nx = 408; Nt = 1201;

Lx = 10; Ly = 12;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[X, Y] = meshgrid(x, y);

A = 50;

% Set up figure
figure
hSurf = surf(X, Y, u(:,:,1));
shading interp
colormap parula
axis tight
axis([min(x) max(x) min(y) max(y) -A A])
xlabel('x'); ylabel('y'); zlabel('Amplitude');

% Animate
for k = 1:Nt
    set(hSurf, 'ZData', u(:,:,k));
    title(sprintf('Frame %d / %d', k, Nt));
    drawnow
end