function [x_best, y_best, x_out, s_out, alpha_out] = powells_conjugate_directions(f, x0, search_radius, convergence_radius)
    arguments 
        f  % Function
        x0 (:,1) {mustBeNumeric}
        search_radius (1,1) {mustBePositive} = 0.25
        convergence_radius (1,1) {mustBePositive} = 0.01
    end

    % Logs the first LOG_OUTP number of line searches to x_out, s_out and a_out
    % corresponding to starting point, search direction and result (alpha)
    LOG_OUTP = 100;
    DISPLAY = false;

    N = length(x0);  % Dimensionality of the problem
    search_vectors = eye(N);  % Initialize with axis-aligned orthonormal basis

    % Initialize variables
    x = x0;
    x_best = x;
    y_best = f(x);

    total_func_count = 0;

    % Initialize detailed output
    search_index = 1;
    x_out = zeros([LOG_OUTP N]);
    s_out = zeros([LOG_OUTP N]);
    alpha_out = zeros([LOG_OUTP 1]);

    for cycle = 1:100
        max_stepsize = 0;
        max_stepsize_i = 1;

        x_temp = x;
        for i = 1:N
            s_i = search_vectors(i,:)';
            % Find optimal stepsize
            [alpha_i, y_i, ~, outp] =...
                fminbnd(@(alpha)f(x_temp+s_i*alpha), -search_radius, search_radius);
            total_func_count = total_func_count + outp.funcCount;
            step_size = vecnorm(alpha_i * s_i);

            if search_index < LOG_OUTP
                x_out(search_index,:) = x_temp;
                s_out(search_index,:) = s_i;
                alpha_out(search_index) = alpha_i;
                search_index = search_index + 1;
            end

            if step_size > max_stepsize
                % Remember the highest stepsize
                max_stepsize = step_size;
                max_stepsize_i = i;
            end
            x_temp = x_temp + alpha_i * s_i;
            if y_i < y_best  
                % Always remember the best visited splot in case it
                % jumps out of a local minimum after converging (happens sometimes)
                x_best = x_temp;
                y_best = y_i;
            end
        end
        % Replace vector of highest stepsize with normalized search direction
        search_vectors(max_stepsize_i,:) = (x_temp - x) / vecnorm(x_temp - x);

        if DISPLAY
            fprintf("y = %f, delta = %f\n", f(x), max_stepsize);
        end

        x = x_temp;

        if max_stepsize < convergence_radius
            %fprintf("%d function evals\n", total_func_count);
            return
        end
    end
end
