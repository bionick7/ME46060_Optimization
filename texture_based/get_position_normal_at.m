function [positions, directions] = get_position_normal_at(pt, maps)
    % Retrieve position and orientation given the uv and the maps
    % Can easily be enhanced to have arbitrary rotations stored in pt
    pt_adj = pt(:, [2 1]);  % Flips x and y coordinates
    rotation = squeeze(pt(:, 3)) * pi;
    positions = squeeze(maps.position_interpolant(pt_adj));
    normal = squeeze(maps.normal_interpolant(pt_adj));
    
    % Gram-schmidt process to project onto y plane
    tang1_direction = [0 1 0] - normal .* repmat(normal(:,2), [1 3]); % Find tangent
    tang1_direction = tang1_direction ./ repmat(vecnorm(tang1_direction, 2, 2), [1 3]); % Normalization
    tang2_direction = cross(normal, tang1_direction);
    % normal, tang1_direction and tang2_direction form an orthonormal basis
    directions = tang1_direction .* sin(rotation) + tang2_direction .* cos(rotation);
end
