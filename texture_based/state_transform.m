function state = state_transform(x)
    arguments 
        x (9, 1) double 
    end
    % Compleates state vector (with symmetries, modulus)

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