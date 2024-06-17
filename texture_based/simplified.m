echo off;

PLOT_ANNEAHLING_OUTP = false;
PCD_OPTIMIZATION = false;
SHOW_PLOT = true;
SHOW_GRID = true;
SHOW_PROFILE = true;

% Model definition
x_half = [0.50 0.48 0.44 0.38 0.38 0.38 0.32 0.32 0.44 0.47 0.50];
y_half = [0.00 0.00 0.11 0.35 0.52 0.70 0.87 0.94 0.94 1.00 1.00];
xs = [x_half(1:10) flip(1.0 - x_half)];
ys = [y_half(1:10) flip(y_half)];
t = (0:1:length(xs)-1) / (length(xs)-1);
model = struct('Xs', xs, 'Ys', ys, 'Ts', t);

tt = 0:0.001:1;
x0 = rand([2 1]);

if PLOT_ANNEAHLING_OUTP
    [x, y, step_info] = simulated_annealing(@(x)objective_func(x, model), x0, @(x) -1, 1000);
    figure(1)
    plot(1:length(step_info), step_info(1,:));
    hold on;
    plot(1:length(step_info), step_info(2,:)/100);
    %plot(1:length(step_info), step_info(3,:)/100);
    hold off;
end

if PCD_OPTIMIZATION
    x = x0;
    y = objective_func(x, model);
    for i = 1:10
        x0 = rand([2 1]);
        [xloc, yloc] = powells_conjugate_directions(@(x)objective_func(x, model), x0);
        if yloc < y
            x = xloc;
            y = yloc;
        end
        fprintf("y = %f\n", yloc);
    end
else
    x = [0.5543 -0.0436]';
    y = objective_func(x, model);
end

x_mod = mod([x' 1-x(1)]', ones([3 1]));
A = get_rcs_matrix(x_mod, model);

xx = pchip(t, xs, tt);
yy = pchip(t, ys, tt);

if SHOW_PLOT
    figure(2);
    hold on;
    plot(xs,ys,'rx');
    plot(xx,yy,'r-');

    % Retrieve state positions/directions
    px = pchip(t, xs, x_mod);
    py = pchip(t, ys, x_mod);
    A = get_rcs_matrix(x_mod, model);
    dx = A(1,:);
    dy = A(2,:);
    plot(px, py,'bo');  % Positions
    for i=1:3
        plot([px(i)+dx(i)*0.1, px(i)-dx(i)*0.1], ...
             [py(i)+dy(i)*0.1, py(i)-dy(i)*0.1], 'b-')
    end

    xlabel("x");
    ylabel("y");
    ylim([-0.1 1.1]);
    axis equal;
    hold off;
end

Y_RANGE = [-1.1 0.1];

if SHOW_GRID
    N = 300;
    grid = zeros(N);
    for x_p = 1:N/2
        for y_p = 1:N/2
            f = objective_func([x_p y_p]' / N, model);
            % Symmetry
            grid(x_p, y_p) = f;
            grid(N-x_p, y_p) = f;
            grid(x_p, N-y_p) = f;
            grid(N-x_p, N-y_p) = f;
        end
    end
    figure(3);
    %imshow(grid, 'XData', [0 1], 'YData', [0 1], 'DisplayRange', Y_RANGE);
    [mesh_xx, mesh_yy] = meshgrid((1:N)/N);
    surf(mesh_xx, mesh_yy, grid, "EdgeColor", "none");
    axis equal;
    xlabel("x_2");
    ylabel("x_1");
    %hold on;
    %scatter(mod(x(2), 1), mod(x(1), 1), 'rx');
    %hold off;
end

if SHOW_PROFILE
    xx = (1:1000) / 1000;
    ff = zeros(1000);
    for i = 1:1000
        x_p = xx(i);
        ff(i) = objective_func([x_p 0.5]', model);
    end

    figure(4);
    plot(xx, ff);
    ylim(Y_RANGE);
    xlabel("x_1");
    ylabel("f(x_1, 0.5)");
end


function y = objective_func(state, model)
    state = mod([state' 1-state(1)]', ones([3 1]));
    A = get_rcs_matrix(state, model);
    %if abs(det(A)) < 1e-5
    %    y = 1e5;
    %else
    %    F = linsolve(A, [0 0 1]');
    %    y = double(norm(reshape(F, [], 1), 20));
    %end
    %y = -norm([1/y 0.1]);
    %y = dot(F, F.*[1 1 -1]');
    y = -norm([det(A) 0.1]);
    %y = -abs(det(A));
end

function A = get_rcs_matrix(state, model)
    px = pchip(model.Ts, model.Xs, state);
    py = pchip(model.Ts, model.Ys, state);
    dx = (pchip(model.Ts, model.Xs, state+0.001) - px) / 0.001;
    dy = (pchip(model.Ts, model.Ys, state+0.001) - py) / 0.001;
    l = sqrt(dx.*dx+dy.*dy);
    dx = dx ./ l;
    dy = dy ./ l;
    torques = px.*dy - py.*dx;
    
    A = [dx dy torques]';
end