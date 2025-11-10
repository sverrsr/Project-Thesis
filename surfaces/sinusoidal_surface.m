function x = sinusoidal_surface( y, z, args, flag )
%SINUSOIDAL_SURFACE General surface function for a bi-axial sinusoidal profile.
%
%   x = SINUSOIDAL_SURFACE( y, z, args, flag ) implements the interface
%   expected by the GeneralLens class. The surface height is defined as a
%   sum of sine waves along the Y and Z directions. The args cell array must
%   contain the amplitudes and spatial periods for each axis:
%       args{1} - amplitude along Y (in the same length units as the bench)
%       args{2} - period along Y
%       args{3} - (optional) amplitude along Z. If omitted, args{1} is used.
%       args{4} - (optional) period along Z. If omitted, args{2} is used.
%
%   When flag == 0 the function returns the surface sag (X coordinate)
%   evaluated at the supplied Y and Z coordinates. Otherwise it returns a
%   matrix of surface normals pointing in the +X direction.
%
%   This surface function is used by example17 to create a reflective
%   sinusoidal free-form surface for FSGD demonstrations.
%
% Copyright: 2024 Optometrika contributors

amp_y = args{ 1 };
per_y = args{ 2 };
if numel( args ) >= 3 && ~isempty( args{ 3 } )
    amp_z = args{ 3 };
else
    amp_z = amp_y;
end
if numel( args ) >= 4 && ~isempty( args{ 4 } )
    per_z = args{ 4 };
else
    per_z = per_y;
end

ky = 0;
if isfinite( per_y ) && per_y ~= 0
    ky = 2 * pi / per_y;
end
kz = 0;
if isfinite( per_z ) && per_z ~= 0
    kz = 2 * pi / per_z;
end

if flag == 0
    x = amp_y * sin( ky * y ) + amp_z * sin( kz * z );
else
    dfd_y = amp_y * ky * cos( ky * y );
    dfd_z = amp_z * kz * cos( kz * z );

    dfd_y = dfd_y( : );
    dfd_z = dfd_z( : );
    ones_vec = ones( numel( dfd_y ), 1 );
    normals = [ ones_vec, -dfd_y, -dfd_z ];
    normals = normals ./ sqrt( sum( normals.^2, 2 ) );
    
    disp(size(normals));
    disp(head(normals));
    x = normals;
end
end
