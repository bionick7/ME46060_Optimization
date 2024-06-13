% Store the results in a json file, so that blender can present it again (for sanity check)
function save_results(x, seed, maps)
    state = state_transform(x);
    [rr, nn] = get_position_normal_at(state, maps);
    A = get_rcs_matrix(x, maps);
    y = -double(det(A)^2);
    F = linsolve(A, vertcat(eye(3), zeros(3)));
    result = struct('seed', seed, 'state', state, 'y', y, 'positions', rr, 'normals', nn, 'F', F);
    encoded = jsonencode(result);
    
    fname = sprintf('../outputs/outp_%d.json', seed);
    fid = fopen(fname,'w');
    fprintf(fid,'%s',encoded);
    fclose(fid);
end

function state = state_transform(x)
    state = zeros([6, 3]);
    state(1:3,:) = mod(reshape(x, [3 3]), ones([3 3]));
    state(4,:) = state(1,:);
    state(5,:) = state(2,:);
    state(6,:) = state(3,:);
    for i = 4:6
        state(i,1) = 1 - state(i,1);
        state(i,3) = - state(i,3);
    end
end

function A = get_rcs_matrix(x, maps)
    state = state_transform(x);    
    [rr, nn] = get_position_normal_at(state, maps);
    A = [cross(rr, nn), nn]';
end
