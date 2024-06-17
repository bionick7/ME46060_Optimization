% Store the results in a json file, so that blender can present it again in 3d (to inspect results)
% Also handy to refer to previous results as it saves the seed
function save_results(x, seed, maps)
    state = state_transform(x);
    [rr, nn] = get_position_normal_at(state, maps);
    A = get_rcs_matrix(x, maps);
    y = -double(norm([det(A), 0.1]));
    F = linsolve(A, vertcat(eye(3), zeros(3)));
    result = struct('seed', seed, 'state', state, 'y', y, 'positions', rr, 'normals', nn, 'F', F);
    encoded = jsonencode(result);
    
    fname = sprintf('../outputs/outp_%d.json', seed);
    fid = fopen(fname,'w');
    fprintf(fid,'%s',encoded);
    fclose(fid);
end

function A = get_rcs_matrix(x, maps)
    state = state_transform(x);    
    [rr, nn] = get_position_normal_at(state, maps);
    A = [cross(rr, nn), nn]';
end
