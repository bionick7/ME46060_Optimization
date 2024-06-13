function [x_min, y_min, step_info] = simulated_annealing(f, x0, constraint, iterations)
    N = size(x0);
    T0 = 100;

    PRINT_ITER = false;

    x = x0;
    y = f(x0);

    x_min = x;
    y_min = y;

    step_info = zeros([3, iterations]);   % y, T
    stagnation_counter = 0;

    T = T0;

    for i = 1:iterations
        % Generate new point
        rand_dir = rand([N 1]);
        rand_dir = rand_dir / norm(rand_dir);
        x_new = mod(x + rand_dir * T * 0.1, 1);
        gen_iter_count = 0;
        while constraint(x_new) > 0
            rand_dir = rand([N 1]);
            rand_dir = rand_dir / norm(rand_dir);
            x_new = mod(x + rand_dir * T * 0.1, 1);
            gen_iter_count = gen_iter_count + 1;
            if gen_iter_count > 100
                x_new = x_min;
            end
        end

        y_new = f(x_new);

        if y_new < y_min
            y_min = y_new;
            x_min = x_new;
        end
        
        accept_proabability = exp((y - y_new) / T * 400);
        if rand([1 1]) < accept_proabability && y_new < 0
            if abs(y_new - y) > 0.01  % Only count significant change
                stagnation_counter = 0;
            else
                stagnation_counter = stagnation_counter + 1;
            end
            x = x_new;
            y = y_new;
        else
            stagnation_counter = stagnation_counter + 1;
        end


        if PRINT_ITER && mod(i, 10) == 1
            fprintf("i: %d |x = (%f, %f); y = %f; T = %f\n", i, x_new(1), x_new(2), y_new, T);
        end

        step_info(:,i) = [y T stagnation_counter];

        % Decrease temperature
        T = T * 0.95;

        if stagnation_counter > 30 && (T < 1 || y > y_min)
            %T0 = T0 * 0.95;
            T = T0/2;
        elseif stagnation_counter > 100
            return
        end
    end
end