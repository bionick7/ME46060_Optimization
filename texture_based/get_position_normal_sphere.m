function [position, direction] = get_position_normal_sphere(pt)
    % Retrieve position and orientation given the uv and the maps
    % Can easily be enhanced to have arbitrary rotations stored in pt

    azimuth = (pt(:, 1) - 0.5) * 2.0 * pi;
    elevation = (pt(:, 2) - 0.5) * pi;
    
    normal = [ cos(azimuth).*cos(elevation) sin(azimuth).*cos(elevation) sin(elevation) ];
    position = normal;
    rotation = pt(:,3);

    % Gram-schmidt process to project onto y plane
    tang1_direction = [0 1 0] - normal .* repmat(normal(:,2), [1 3]); % Find tangent
    tang1_direction = tang1_direction ./ repmat(vecnorm(tang1_direction, 2, 2), [1 3]); % Normalization
    tang2_direction = cross(normal, tang1_direction);
    % normal, tang1_direction and tang2_direction form an orthonormal basis
    direction = tang1_direction .* sin(rotation) + tang2_direction .* cos(rotation);
end
