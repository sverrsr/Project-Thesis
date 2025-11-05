function varargout = example19Peak()
%EXAMPLE19PEAK Demonstration bench using a PEAKS-based GeneralLens surface.
%
%   This example constructs a GeneralLens surface driven by MATLAB's
%   PEAKS profile. A collimated bundle of rays illuminates the free-form
%   surface and the resulting distribution is captured on a screen. The
%   example illustrates how arbitrary analytic surfaces can be wrapped and
%   explored inside the Optometrika bench using GeneralLens.
%
% Copyright: 2024 Optometrika contributors

% Create the optical bench container
bench = Bench;

% Configure the PEAKS surface parameters
aperture = 80;          % clear aperture diameter (mm)
height_scale = 6;       % scale of the sag values (mm)
lateral_scale = 20;     % lateral scale applied before evaluating PEAKS

peaks_surface_element = GeneralLens( [ 0 0 0 ], aperture, 'peaks_surface', ...
    { 'air' 'mirror' }, height_scale, lateral_scale );
bench.append( peaks_surface_element );

% Place a capture screen downstream of the surface
screen_distance = -20;   % mm along +X
screen_size = 180;       % mm side length
screen = Screen( [ screen_distance 0 1 ], screen_size, screen_size, 512, 512 );
screen.rotate( [ 0 1 0 ], pi );
bench.append( screen );

% Generate a collimated bundle of rays aimed at the surface
nrays = 100;
source_pos = [ -120 0 0 ];
incident_dir = [ 1 0 0 ];
beam_diameter = aperture * 0.95;

rays_in = Rays( nrays, 'collimated', source_pos, incident_dir, beam_diameter, 'hexagonal' );

fprintf( 'Tracing rays through the PEAKS-based surface...\n' );
rays_out = bench.trace( rays_in );

% Visualise the bench and ray paths
bench.draw( rays_out, 'lines', [], 1.5 );
title( 'GeneralLens surface shaped by MATLAB PEAKS', 'Color', 'w' );

% Display the intensity distribution collected on the screen
figure( 'Name', 'PEAKS screen capture', 'NumberTitle', 'Off' );
imagesc( screen.image );
axis image;
colormap hot;
colorbar;
set( gca, 'YDir', 'normal' );
title( 'Illumination captured after the PEAKS surface' );
xlabel( 'Screen Y bins' );
ylabel( 'Screen Z bins' );

if nargout >= 1
    varargout{ 1 } = screen;
end
if nargout >= 2
    varargout{ 2 } = rays_out;
end
if nargout >= 3
    varargout{ 3 } = bench;
end
if nargout >= 4
    varargout{ 4 } = peaks_surface_element;
end

end
