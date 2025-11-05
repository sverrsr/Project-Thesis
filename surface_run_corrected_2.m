function varargout = surface_run_corrected_2()
%SURFACE_RUN_CORRECTED Trace rays against the interpolated sampled mirror.
%   This helper mirrors the workflow from SURFACE_RUN but fixes two issues:
%     1) The imported sag is flipped so the surface faces the incoming beam.
%     2) Rays that miss the tabulated footprint no longer clip to the edge;
%        they are marked as lost and stop contributing to the screen image.
%
%   The script loads the mesh and sag data, builds gridded interpolants for
%   the sag and its slopes, feeds them to the SURFACE_LENS_BOUNDED helper,
%   and launches a reflective ray trace.
%
%   Outputs (optional):
%       screen   - Screen object capturing the irradiance map
%       rays_out - Cell array of traced rays
%       bench    - Bench containing the optical elements
%       mirror   - GeneralLens instance representing the sampled mirror
%
%   Copyright: 2024 Optometrika contributors

clearvars; close all; clc;

load surfaceData1200.mat
load surfMesh.mat

X = xMesh; %#ok<NODEF>
Y = yMesh; %#ok<NODEF>
Z = surfaceData1200; %#ok<NODEF>
clear xMesh yMesh surfaceData1200

[xa, ix] = sort( X( 1, : ) );
[ya, iy] = sort( Y( :, 1 ) );

x_limits = [ xa( 1 ), xa( end ) ];
y_limits = [ ya( 1 ), ya( end ) ];
grid_center = [ mean( x_limits ), mean( y_limits ) ];
half_span = 0.5 * [ diff( x_limits ), diff( y_limits ) ];

dx = diff( xa );
dy = diff( ya );
dx_min = min( dx( : ) );
dy_min = min( dy( : ) );
if isempty( dx_min ) || isempty( dy_min ) || ~isfinite( dx_min ) || ~isfinite( dy_min )
    error( 'surface_run_corrected:InvalidGrid', ...
        'Surface grid must contain at least two unique samples per axis.' );
end
min_spacing = min( [ dx_min, dy_min ] );
edge_buffer = 0.5 * min_spacing;

usable_radius = min( half_span ) - edge_buffer;
if usable_radius <= 0
    error( 'surface_run_corrected:InvalidSurfaceBounds', ...
        'Surface data has insufficient span once edge buffer is removed.' );
end

z_limits = [ min( Z( : ) ), max( Z( : ) ) ];
lens_depth = diff( z_limits );

Za = Z( iy, ix );
F = griddedInterpolant( { xa, ya }, Za.', 'linear', 'none' );
Zi = F( X', Y' )';
[ dZdx, dZdy ] = gradient( Zi, xa, ya );
Fdx = griddedInterpolant( { ya, xa }, dZdx, 'linear', 'none' );
Fdy = griddedInterpolant( { ya, xa }, dZdy, 'linear', 'none' );

orientation_settings = struct( ...
    'swap_yz', false, ...
    'flip_y', false, ...
    'flip_z', false, ...
    'invert_sag', true );

lens_args = { F, Fdx, Fdy, grid_center, x_limits, y_limits, orientation_settings };

bench = Bench;
aperture = 2 * usable_radius;
mirror = GeneralLens( [ 0 0 0 ], aperture, 'surface_lens_bounded', { 'air' 'mirror' }, lens_args{:} );
bench.append( mirror );

screen_distance = -max( 20, lens_depth + aperture * 0.75 );
screen_size = max( aperture * 1.25, 64 );
screen = Screen( [ screen_distance 0 1 ], screen_size, screen_size, 512, 512 );
screen.rotate( [ 0 1 0 ], pi );
bench.append( screen );

nrays = 5000;
source_pos = [ -( lens_depth + aperture * 1.5 ) 0 0 ];
incident_dir = [ 1 0 0 ];
beam_diam = aperture * 0.9;
rays_in = Rays( nrays, 'collimated', source_pos, incident_dir, beam_diam, 'hexagonal' );

fprintf( 'Tracing rays through surface\_lens\_bounded ...\n' );
rays_out = bench.trace( rays_in );

bench.draw( rays_out, 'lines', 1, 1.5 );
axis equal;
grid on;
view( 35, 20 );
xlabel( 'X (mm)' ); ylabel( 'Y (mm)' ); zlabel( 'Z (mm)' );
camlight( 'headlight' ); camlight( 'left' ); camlight( 'right' );
lighting gouraud;
title( 'GeneralLens using surface\_lens\_bounded', 'Color', 'w' );

figure( 'Name', 'surface\_lens\_bounded screen capture', 'NumberTitle', 'Off' );
imagesc( screen.image ); axis image; colormap hot; colorbar;
set( gca, 'YDir', 'normal' );
title( 'Illumination after surface\_lens\_bounded' );
xlabel( 'Screen Y bins' ); ylabel( 'Screen Z bins' );

if nargout >= 1, varargout{ 1 } = screen; end
if nargout >= 2, varargout{ 2 } = rays_out; end
if nargout >= 3, varargout{ 3 } = bench; end
if nargout >= 4, varargout{ 4 } = mirror; end

end