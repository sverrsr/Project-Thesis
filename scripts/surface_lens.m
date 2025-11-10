function x = surface_lens( y, z, args, flag )
% surface_lens
%
% Computes either the surface height or the normal vector of a general lens
% profile. If flag == 0 the function returns the lens height x for the given position ( y, z ).
% Otherwise, the function returns the lens normal at this position. 
%
% Usage:
%   x = surface_lens(y, z, args, flag)
%
% Inputs:
%   y, z   - Coordinates in the lens plane.
%   args   - Cell array with surface data and optional parameters:
%              args{1} : Function handle for lens height Z(x, y)    
%              args{2} : Function handle for dZ/dx
%              args{3} : Function handle for dZ/dy
%              args{4} : (optional) [x0, y0] grid center
%              args{5} : (optional) [xmin, xmax] limits
%              args{6} : (optional) [ymin, ymax] limits
%   flag   - Mode selector:
%              0 : return surface height (x = Z)
%              1 : return surface normal [nx, ny, nz]
%
% Outputs:
%   x - If flag == 0: scalar lens height.
%       If flag ~= 0: unit normal vector at (y, z), oriented along +x.
%
% Notes:
% - Keeps evaluations within the defined grid limits.
% - Surface normals are normalized and flipped if not pointing toward +x.


F   = args{1};   % height Z(x,y)
Fdx = args{2};   % dZ/dx, called as Fdx(y,x)
Fdy = args{3};   % dZ/dy, called as Fdy(y,x)

arg_idx = 4;
grid_center = [0, 0];
if numel(args) >= arg_idx && ~isempty(args{arg_idx})
    grid_center = double(args{arg_idx}(:).');
    if numel(grid_center) < 2
        grid_center(2) = 0;
    elseif numel(grid_center) > 2
        grid_center = grid_center(1:2);
    end
    arg_idx = arg_idx + 1;
end

x_limits = [-inf, inf];
y_limits = [-inf, inf];
if numel(args) >= arg_idx && ~isempty(args{arg_idx})
    x_limits = double(args{arg_idx}(:).');
    if numel(x_limits) < 2
        x_limits(2) = x_limits(1);
    elseif numel(x_limits) > 2
        x_limits = x_limits(1:2);
    end
    arg_idx = arg_idx + 1;
end
if numel(args) >= arg_idx && ~isempty(args{arg_idx})
    y_limits = double(args{arg_idx}(:).');
    if numel(y_limits) < 2
        y_limits(2) = y_limits(1);
    elseif numel(y_limits) > 2
        y_limits = y_limits(1:2);
    end
end

% Map lens (y_in,z_in) -> original (x_orig,y_orig)
x_orig = y + grid_center(1);      % original x
y_orig = z + grid_center(2);      % original y

% Keep evaluations inside the tabulated domain to avoid NaNs during the
% numerical intersection search.
if all(isfinite(x_limits))
    x_orig = min(max(x_orig, x_limits(1)), x_limits(2));
end
if all(isfinite(y_limits))
    y_orig = min(max(y_orig, y_limits(1)), y_limits(2));
end


if flag == 0
    % --- Return sag along x (lens axis): x = Z(x_orig, y_orig)
    x = F(x_orig, y_orig);   % <-- your F expects (x,y)

else
    % --- Return unit normal [nx, ny, nz] in lens coordinates
    % Slopes of Z at (x_orig, y_orig).
    % NOTE: Fdx/Fdy are built with {ya, xa}, so call them as (y,x):
    Zx = Fdx(y_orig, x_orig);
    Zy = Fdy(y_orig, x_orig);

    c  = 1 ./ sqrt(1 + Zx.^2 + Zy.^2);  % normalization factor

    nx = c;                 % = n_z,orig
    ny = -Zx .* c;          % = n_x,orig
    nz = -Zy .* c;          % = n_y,orig

    x = [nx, ny, nz];


    % Keep orientation (pointing ~ +x). Flip if needed:
    flipmask = (nx < 0);
    if any(flipmask, 'all')
        x(flipmask,:) = -x(flipmask,:);
    end
end

