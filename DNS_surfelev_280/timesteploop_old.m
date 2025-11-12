clear all; close all; clc;

%% Set Mesh and Surface Properties
load surfElev_280.00.mat;

Z = double(surfElev);
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
[X, Y] = meshgrid((linspace(0, lx, nx)), (linspace(0, ly, ny)));

outDir = 'DNS_SCREENS_1k';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

Nt = 1;
for k = 1:Nt  % or however many surfaces you have
    fprintf('Tracing frame %d of %d ...\n', k, Nt);

    % --- your setup and tracing code here ---
    [screen, rays_out, bench, surf] = DNS(X, Y, Z);

    filename = fullfile(outDir, sprintf('screen_%04d.mat', k));
    screen_image = screen.image; %#ok<NASGU>
    save(filename, 'screen');

    fprintf('Saved %s\n', filename);
    % Visualize

    bench.draw(rays_out, 'lines', 1, 1.5);
    axis equal;
    grid on;
    view(50, 50);
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
    camlight('headlight'); camlight('left'); camlight('right');
    lighting gouraud;
    title('GeneralLens using surface_lens (interpolated surface)', 'Color','w');
    
    figure('Name','surface_lens screen capture','NumberTitle','Off');
    imagesc(screen.image); axis image; colormap hot; colorbar;
    set(gca,'YDir','normal');
    title('Illumination after surface_lens'); xlabel('Screen Y bins'); ylabel('Screen Z bins');
end

