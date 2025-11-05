% Load the time series
u = h5read(fullfile(pwd, 'blob', 'wave.h5'), '/u');   % size: [Ny, Nx, Nt]
[Ny, Nx, Nt] = size(u);

% Build the mesh (must match u(:,:,k) indexing)
Lx = 10; Ly = 12;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[X, Y] = meshgrid(x, y);

% Optional: pre-allocate storage for screen images
store_images = true;
if store_images
    % Get one frame to size the buffer
    [screen0, ~, ~, ~, state] = validation_Blob_simulation(X, Y, u(:,:,1), [], false);
    scr = screen0.image;
    screen_stack = zeros([size(scr), Nt], 'like', scr);
    screen_stack(:,:,1) = scr;
else
    state = [];
    validation_Blob_simulation(X, Y, u(:,:,1), state, true); % show plots if you want
end

% Process remaining frames
for k = 2:Nt
    Zk = u(:,:,k);
    [screen, ~, ~, ~, state] = validation_Blob_simulation(X, Y, Zk, state, false);
    if store_images
        screen_stack(:,:,k) = screen.image;
    end
    % If you want a live view, flip do_plot=true in the call above.
end

% Example: quick playback of the screen intensity (if stored)
if store_images
    figure('Name','Screen intensity over time','NumberTitle','Off');
    him = imagesc(screen_stack(:,:,1)); axis image; colormap hot; colorbar;
    set(gca,'YDir','normal');
    title(sprintf('Frame 1 / %d', Nt)), xlabel('Screen Y bins'), ylabel('Screen Z bins');
    for k = 1:Nt
        set(him, 'CData', screen_stack(:,:,k));
        title(sprintf('Frame %d / %d', k, Nt));
        drawnow;
    end
end
