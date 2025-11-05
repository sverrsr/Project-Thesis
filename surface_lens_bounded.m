function x = surface_lens_bounded_2( y, z, args, flag )
%SURFACE_LENS_BOUNDED GeneralLens wrapper with explicit footprint handling.
%   This variant extends SURFACE_LENS by (1) enforcing the finite footprint
%   of the tabulated surface so rays that miss the sampled area are rejected
%   instead of clamping to the nearest edge sample, and (2) providing simple
%   orientation controls that correct upside-down surfaces when importing
%   mirror sags measured in an opposing coordinate frame.
%
%   The ARGS cell array must contain the same inputs as SURFACE_LENS:
%       ARGS{1} - griddedInterpolant returning the surface sag (X coordinate)
%       ARGS{2} - griddedInterpolant returning d(sag)/dX (evaluated as Fdx(y,x))
%       ARGS{3} - griddedInterpolant returning d(sag)/dY (evaluated as Fdy(y,x))
%       ARGS{4} - 1x2 vector with the [X Y] coordinates of the grid centre
%       ARGS{5} - 1x2 vector with the minimum/maximum X coordinates sampled
%       ARGS{6} - 1x2 vector with the minimum/maximum Y coordinates sampled
%
%   An optional ARGS{7} structure configures orientation fixes:
%       .swap_yz      - swap incoming Y/Z coordinates before sampling
%       .flip_y       - multiply the (post-swapped) Y coordinates by -1
%       .flip_z       - multiply the (post-swapped) Z coordinates by -1
%       .invert_sag  - mirror the sag about the lens plane (x -> -x)
%
%   When FLAG == 0 the function returns the sag. Otherwise it returns unit
%   surface normals pointing in the +X direction. Queries outside the
%   tabulated domain return Inf for the sag and NaN normals to force the ray
%   tracer to treat them as misses.
%
%   Copyright: 2024 Optometrika contributors

if nargin < 4
    flag = 0;
end

if numel( args ) < 6
    error( 'surface_lens_bounded:MissingArgs', ...
        'Expected interpolants, grid centre and axis limits in args.' );
end

F   = args{ 1 };
Fdx = args{ 2 };
Fdy = args{ 3 };

grid_center = double( args{ 4 } );
if numel( grid_center ) < 2
    grid_center( 2 ) = 0;
elseif numel( grid_center ) > 2
    grid_center = grid_center( 1:2 );
end

x_limits = double( args{ 5 } );
y_limits = double( args{ 6 } );

orientation.swap_yz   = false;
orientation.flip_y    = false;
orientation.flip_z    = false;
orientation.invert_sag = false;

if numel( args ) >= 7 && ~isempty( args{ 7 } )
    user = args{ 7 };
    if isstruct( user )
        fields = fieldnames( orientation );
        for k = 1:numel( fields )
            name = fields{ k };
            if isfield( user, name ) && ~isempty( user.( name ) )
                orientation.( name ) = logical( user.( name ) );
            end
        end
    else
        error( 'surface_lens_bounded:InvalidOrientation', ...
            'Orientation settings must be supplied as a struct.' );
    end
end

% Apply orientation adjustments to the incoming coordinates.
y_in = y;
z_in = z;
if orientation.swap_yz
    tmp = y_in;
    y_in = z_in;
    z_in = tmp;
end
if orientation.flip_y
    y_in = -y_in;
end
if orientation.flip_z
    z_in = -z_in;
end

% Map to the original data coordinates.
x_orig = y_in + grid_center( 1 );
y_orig = z_in + grid_center( 2 );

% Determine which queries fall inside the sampled domain.
inside = x_orig >= x_limits( 1 ) & x_orig <= x_limits( 2 ) & ...
         y_orig >= y_limits( 1 ) & y_orig <= y_limits( 2 );

% Prepare output containers with the requested shape.
out_shape = size( y );

if flag == 0
    sag = Inf( out_shape );
    if any( inside(:) )
        values = F( x_orig( inside ), y_orig( inside ) );
        if orientation.invert_sag
            values = -values;
        end
        sag( inside ) = values;
    end
    x = sag;
    return;
end

% Flag==1: return normals. Outside the grid report NaNs so the caller can
% discard the ray after the intersection routine marks it as a miss.
nx_full = NaN( out_shape );
ny_full = NaN( out_shape );
nz_full = NaN( out_shape );
if any( inside(:) )
    dxd = Fdx( y_orig( inside ), x_orig( inside ) );
    dyd = Fdy( y_orig( inside ), x_orig( inside ) );
    if orientation.invert_sag
        dxd = -dxd;
        dyd = -dyd;
    end
    c = 1 ./ sqrt( 1 + dxd.^2 + dyd.^2 );
    nx = c;
    ny = -dxd .* c;
    nz = -dyd .* c;
    flipmask = nx < 0;
    if any( flipmask )
        nx( flipmask ) = -nx( flipmask );
        ny( flipmask ) = -ny( flipmask );
        nz( flipmask ) = -nz( flipmask );
    end
    nx_full( inside ) = nx;
    ny_full( inside ) = ny;
    nz_full( inside ) = nz;
end

x = [ nx_full( : ), ny_full( : ), nz_full( : ) ];

end
surface_run_corrected.m
Ny
+113-0
function varargout = surface_run_corrected()
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