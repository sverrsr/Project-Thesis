function x = peaks_surface( y, z, args, flag )
%PEAKS_SURFACE GeneralLens surface based on MATLAB's PEAKS function.
%
%   x = PEAKS_SURFACE( y, z, args, flag ) implements the interface required
%   by the GeneralLens class using the analytic expression behind MATLAB's
%   built-in PEAKS function. The surface height is scaled and sampled in the
%   YZ plane and returned as the X coordinate.
%
%   The args cell array supports two optional parameters:
%       args{1} - height scale factor applied to the PEAKS surface (default 6)
%       args{2} - lateral scale factor controlling the footprint size
%                 (default 20). The input Y and Z coordinates are divided by
%                 this factor before evaluating PEAKS.
%
%   When flag == 0 the function returns the surface sag (X coordinate).
%   Otherwise it returns a matrix of surface normals pointing in the +X
%   direction.
%
% Copyright: 2024 Optometrika contributors

if nargin < 3 || isempty( args )
    args = { [] [] };
end

if numel( args ) < 2
    args{ 2 } = [];
end

height_scale = args{ 1 };
if isempty( height_scale )
    height_scale = 6;
end

lateral_scale = args{ 2 };
if isempty( lateral_scale )
    lateral_scale = 20;
end

X = y ./ lateral_scale;
Y = z ./ lateral_scale;

[ F, Fx, Fy ] = local_peaks_eval( X, Y );

if flag == 0
    x = height_scale * F;
else
    dfd_y = height_scale / lateral_scale * Fx;
    dfd_z = height_scale / lateral_scale * Fy;

    dfd_y = dfd_y( : );
    dfd_z = dfd_z( : );
    ones_vec = ones( numel( dfd_y ), 1 );
    normals = [ ones_vec, -dfd_y, -dfd_z ];
    normals = normals ./ sqrt( sum( normals.^2, 2 ) );

    x = normals;
end

end

function [ F, Fx, Fy ] = local_peaks_eval( X, Y )
%LOCAL_PEAKS_EVAL Analytic evaluation of PEAKS and its gradients.
%   This helper reproduces MATLAB's PEAKS surface and provides its partial
%   derivatives with respect to X and Y.

T1 = 3 * ( 1 - X ).^2 .* exp( -( X.^2 ) - ( Y + 1 ).^2 );
T2 = -10 * ( X / 5 - X.^3 - Y.^5 ) .* exp( -X.^2 - Y.^2 );
T3 = -( 1 / 3 ) * exp( -( X + 1 ).^2 - Y.^2 );

F = T1 + T2 + T3;

% Derivatives with respect to X
A = 3 * ( 1 - X ).^2;
B = exp( -( X.^2 ) - ( Y + 1 ).^2 );
A_prime = -6 * ( 1 - X );
B_prime_x = B .* ( -2 .* X );

T1_dx = B .* ( A_prime + A .* ( -2 .* X ) );

C = -10 * ( X / 5 - X.^3 - Y.^5 );
D = exp( -X.^2 - Y.^2 );
C_prime_x = -2 + 30 .* X.^2;
D_prime_x = D .* ( -2 .* X );

T2_dx = D .* ( C_prime_x + C .* ( -2 .* X ) );

T3_dx = -( 1 / 3 ) * exp( -( X + 1 ).^2 - Y.^2 ) .* ( -2 .* ( X + 1 ) );

Fx = T1_dx + T2_dx + T3_dx;

% Derivatives with respect to Y
B_prime_y = B .* ( -2 .* ( Y + 1 ) );
T1_dy = A .* B_prime_y;

C_prime_y = 50 .* Y.^4;
D_prime_y = D .* ( -2 .* Y );
T2_dy = D .* ( C_prime_y + C .* ( -2 .* Y ) );

T3_dy = -( 1 / 3 ) * exp( -( X + 1 ).^2 - Y.^2 ) .* ( -2 .* Y );

Fy = T1_dy + T2_dy + T3_dy;
end
