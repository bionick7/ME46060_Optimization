function [x_best, y_best] = powells_conjugate_directions(f, x0)
    N = length(x0);
    search_vectors = eye(N);
    x = x0;
    x_best = x;
    y_best = f(x);

    for iter = 1:100
        max_stepsize = 0;
        max_stepsize_i = 1;

        x_temp = x;
        for i = 1:N
            s_i = search_vectors(i,:)';
            [alpha_i, y_i] = fminbnd(@(alpha)f(x_temp+s_i*alpha), -0.2, 0.2);
            step_size = vecnorm(alpha_i * s_i);
            if step_size > max_stepsize
                max_stepsize = step_size;
                max_stepsize_i = i;
            end
            x_temp = x_temp + alpha_i * s_i;
            if y_i < y_best
                x_best = x_temp;
                y_best = y_i;
            end
        end
        search_vectors(max_stepsize_i,:) = (x_temp - x) / vecnorm(x_temp - x);
        %fprintf("y = %f, delta = %f\n", f(x), max_stepsize);
        x = x_temp;

        if max_stepsize < 1e-2
            return
        end
    end
end
