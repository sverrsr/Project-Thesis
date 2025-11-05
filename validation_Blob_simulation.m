function [screen, rays_out, bench, surf, state] = examplesurface_lensRun(X, Y, Z, state, do_plot)
% examplesurface_lensRun
% X, Y: meshgrid-sized coordinates (Ny x Nx)
% Z   : surface heights at this timestep (Ny x Nx)
% state: pass [] on first call; reuse returned state afterwards
% do_plot: logical, optional (default: false). If true, draws bench & screen.
%
% Returns: screen, rays_out, bench, surf, state (updated)

if nargin < 5, do_plot = false; end
is_init = (nargin < 4) || isempty(state);

% --- Prepare sorted axes & bounds (meshgrid -> ndgrid bookkeeping) ---
[xa, ix] = sort(X(1,:));         % 1 x Nx  (monotonic)
[ya, iy] = sort(Y(:,1));         % Ny x 1  (monotonic)

x_limits = [xa(1), xa(end)];
y_limits = [ya(1), ya(end)];
grid_center = [mean(x_limits), mean(y_limits)];
half_span   = 0.5 * [diff(x_limits), diff(y_limits)];

dx = diff(xa); dy = diff(ya);
dx_min = min(dx(:)); dy_min = min(dy(:));
if isempty(dx_min) || isempty(dy_min) || ~isfinite(dx_min) || ~isfinite(dy_min)
    error('surface_run:InvalidGrid', 'Surface grid must contain at least two unique samples per axis.');
end
min_spacing = min([dx_min, dy_min]);
edge_buffer = 0.5 * min_spacing;

usable_radius = min(half_span) - edge_buffer;
if usable_radius <= 0
    error('surface_run:InvalidSurfaceBounds', 'Surface span insufficient once edge buffer is removed.');
end
usable_half_span = half_span - edge_buffer;

% Reorder Z to match sorted axes (meshgrid order -> Za)
Za = Z(iy, ix);                   % still (Ny x Nx) in meshgrid orientation

% --- Build or update interpolants (F: height, Fdx/Fdy: slopes) ---
if is_init
    % Height interpolant uses ndgrid order {xa,ya} with values transposed
    F   = griddedInterpolant({xa, ya}, Za.', 'linear', 'none');  % F(x,y)
    % Gradients on the meshgrid-aligned Z, spacing: X then Y (cols->xa, rows->ya)
    [dZdx, dZdy] = gradient(Z, xa, ya);
    % Slope interpolants use {ya, xa} and meshgrid-sized arrays
    Fdx = griddedInterpolant({ya, xa}, dZdx, 'linear', 'none');
    Fdy = griddedInterpolant({ya, xa}, dZdy, 'linear', 'none');

    lens_args = {F, Fdx, Fdy, grid_center, x_limits, y_limits};

    % --- Build bench (once) ---
    bench = Bench;

    ap_half_y = usable_half_span(1);
    ap_half_z = usable_half_span(2);
    rect_aperture = [0; 0; 2*ap_half_y; 2*ap_half_z];  % full widths

    surf = GeneralLens([0 0 0], rect_aperture, 'surface_lens', {'mirror','air'}, lens_args{:});
    bench.append(surf);

    % Screen along +X
    aperture = 2 * usable_radius;
    screen_distance = 200;                     % mm
    screen_size = max(aperture * 1.25, 64);    % mm
    screen = Screen([screen_distance 0 0], screen_size, screen_size, 512, 512);
    screen.rotate([1 0 0], pi);                % face back to optic
    bench.append(screen);

    % Source: collimated beam along +X
    nrays = 10;
    source_distance = 300;
    source_pos   = [source_distance 0 0];
    incident_dir = [-1 0 0];

    % Square beam slightly smaller than the larger side of the rectangular aperture
    beam_side = 0.98 * max(2*ap_half_y, 2*ap_half_z);
    rays_in = Rays(nrays, 'collimated', source_pos, incident_dir, beam_side, 'random');

    % Trace once
    rays_out = bench.trace(rays_in);

    % Optional plots (first frame)
    if do_plot
        figure('Name','surface_lens screen capture','NumberTitle','Off');
        imagesc(screen.image); axis image; colormap hot; colorbar;
        set(gca,'YDir','normal');
        title('Illumination after surface_lens'); xlabel('Screen Y bins'); ylabel('Screen Z bins');
    end

    % Pack state to reuse next frames
    state = struct( ...
        'F', F, 'Fdx', Fdx, 'Fdy', Fdy, ...
        'xa', xa, 'ya', ya, 'ix', ix, 'iy', iy, ...
        'grid_center', grid_center, 'x_limits', x_limits, 'y_limits', y_limits, ...
        'bench', bench, 'screen', screen, 'surf', surf, 'rays_in', rays_in, ...
        'usable_half_span', usable_half_span, 'edge_buffer', edge_buffer);

else
    % Update only the values (fast path)
    state.F.Values = Za.';                     % height
    [dZdx, dZdy] = gradient(Z, xa, ya);        % meshgrid orientation
    state.Fdx.Values = dZdx;                   % slopes
    state.Fdy.Values = dZdy;

    bench  = state.bench;
    screen = state.screen;
    surf   = state.surf;

    % Re-trace with the same rays
    rays_out = bench.trace(state.rays_in);

    if do_plot
        % if you created a figure yourself, just refresh it here if needed
        figure(findobj('Name','surface_lens screen capture'));
        imagesc(screen.image); axis image; colormap hot; colorbar;
        set(gca,'YDir','normal');
        title('Illumination after surface_lens'); xlabel('Screen Y bins'); ylabel('Screen Z bins');
        drawnow;
    end
end
end
