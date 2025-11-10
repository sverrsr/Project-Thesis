function generalBench(lens_args)

%% --- Build the bench (same layout as your example) ---
bench = Bench;

aperture = 2 * usable_radius;  % mm, matches the usable surface data
% Use 'air' to 'mirror' for a reflective test (no dispersion setup needed).
% For a transmissive lens, swap 'mirror' -> a glass name present in your material set (e.g., 'bk7').

usable_half_span = half_span - edge_buffer;
ap_half_y = usable_half_span(1);
ap_half_z = usable_half_span(2);

rect_aperture = [0; 0; 2*ap_half_y; 2*ap_half_z];   % FULL widths (required)

surf = GeneralLens([0 0 0], rect_aperture, 'surface_lens', {'mirror','air'}, lens_args{:});
%surf.rotate([0 0 1], pi/4);   % face back toward the optic

bench.append(surf);

% Screen downstream (along +X)
screen_distance = 200;  % mm along +X
screen_size = max(aperture * 1.25, 64);                    % mm
screen = Screen([screen_distance 0 1], screen_size, screen_size, 512, 512);
screen.rotate([1 0 0], pi);   % face back toward the optic
bench.append(screen);

% Collimated beam aimed along +X
nrays = 1000;
source_distance = 300;
source_pos   = [source_distance 0 0];
incident_dir = [-1 0 0];


% beam: square side >= the larger rectangle side
ap_half_y = usable_half_span(1);
ap_half_z = usable_half_span(2);

rect_wy = 2*ap_half_y;
rect_wz = 2*ap_half_z;
beam_side = 0.98 * max(rect_wy, rect_wz);  % slightly smaller than the larger side

rays_in = Rays(nrays, 'collimated', source_pos, incident_dir, beam_side, 'random');

figure('Name','Launch footprint (source plane)');
scatter(rays_in.r(:,2), rays_in.r(:,3), 6, 'filled'); axis equal; grid on;
xlabel('Y at source'); ylabel('Z at source'); title('Rays launch footprint');

fprintf('Tracing rays through surface_lens ...\n');
rays_out = bench.trace(rays_in);

% Print screen geometry
fprintf('Surface radius R: %.3f\n', surf.R);
fprintf('Surface normal n: [%.4f %.4f %.4f]\n', surf.n);
fprintf('Surface position r: [%.4f %.4f %.4f]\n', surf.r);
% Print screen geometry
fprintf('Screen radius R: %.3f\n', screen.R);
fprintf('Screen normal n: [%.4f %.4f %.4f]\n', screen.n);
fprintf('Screen position r: [%.4f %.4f %.4f]\n', screen.r);
% Print screen geometry
fprintf('Beam normal n: [%.4f %.4f %.4f]\n', incident_dir);
fprintf('Beam position r: [%.4f %.4f %.4f]\n', source_pos);

% Visualize
bench.draw(rays_out, 'lines', 1, 1.5);
axis equal;
grid on;
view(35, 20);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
camlight('headlight'); camlight('left'); camlight('right');
lighting gouraud;
title('GeneralLens using surface_lens (interpolated surface)', 'Color','w');

figure('Name','surface_lens screen capture','NumberTitle','Off');
imagesc(screen.image); axis image; colormap hot; colorbar;
set(gca,'YDir','normal');
title('Illumination after surface_lens'); xlabel('Screen Y bins'); ylabel('Screen Z bins');

if nargout >= 1, varargout{1} = screen; end
if nargout >= 2, varargout{2} = rays_out; end
if nargout >= 3, varargout{3} = bench; end
if nargout >= 4, varargout{4} = surf; end
end


end
