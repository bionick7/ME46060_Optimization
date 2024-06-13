clc;
%clear all;
%close all;
echo off;

PLOT_ANNEAHLING_OUTP = true;
SHOW_PLOT = true;
SHOW_GRID = false;
SHOW_PROFILE = false;

x_half = [0.50 0.48 0.44 0.38 0.38 0.38 0.32 0.32 0.44 0.47 0.50];
y_half = [0.00 0.00 0.11 0.34 0.52 0.70 0.87 0.94 0.94 1.00 1.00];
xs = [x_half(1:10) flip(1.0 - x_half)];
ys = [y_half(1:10) flip(y_half)];
t = (0:1:length(xs)-1) / (length(xs)-1);

tt = 0:0.001:1;

model = struct('Xs', xs, 'Ys', ys, 'Ts', t);
x0 = rand([1 2]);
%[x, y] = ga(@(x)objective_func(x, model), 2)

%problem = createOptimProblem('fmincon','objective', @(x)objective_func(x, model),...
    %                             'x0',x0,'lb',[0 0],'ub', [1 1]);
    %          %'nonlcon',apertureConstraint_x,'options',opts);
    %gs = GlobalSearch;
    %[x,y] = run(gs, problem)
    
%[x, y] = simulannealbnd(@(x)objective_func(x, model), x0, [], [], optimset('Display', 'iter'))
[x, y, step_info] = simulated_annealing(@(x)objective_func(x, model), x0);
[x, y] = fminunc(@(x)objective_func(x, model), x)

if PLOT_ANNEAHLING_OUTP
    figure(1)
    plot(1:length(step_info), step_info(1,:));
    hold on;
    plot(1:length(step_info), step_info(2,:)/100);
    %plot(1:length(step_info), step_info(3,:)/100);
    hold off;
end

% Works: 
% genetic algorythms
% simulated aneahling
% GlobalSearch

% Does not work:
% particleswarm

x_mod = mod([1-x(1) x]', ones([3 1]));
A = get_rcs_matrix(x_mod, model);

xx = interp1(t, xs, tt);
yy = interp1(t, ys, tt);

if SHOW_PLOT
    figure(2);
    hold on;
    plot(xs,ys,'rx');
    plot(xx,yy,'r-');
    plot(interp1(t, xs, x_mod), interp1(t, ys, x_mod),'bo');
    axis equal;
    hold off;
end

Y_RANGE = [-1 1];

if SHOW_GRID
    grid = zeros(2);
    N = 100;
    for x_p = 1:N
        for y_p = 1:N
            grid(x_p, y_p) = objective_func([x_p y_p] / N, model);
        end
    end
    figure(3);
    %imshow(grid, 'XData', [0 1], 'YData', [0 1], 'DisplayRange', Y_RANGE);
    %hold on; 
    [mesh_xx, mesh_yy] = meshgrid((1:N)/N);
    contour(mesh_xx, mesh_yy, grid);
    %scatter(mod(x(1), 1), mod(x(2), 1), 'bo');
    %hold off;
end

%figure();
%fcontour(@(x,y)domain([x, y], model), [0 1]);
%fcontour(@(x,y)x.^2+y.^2, [0 1]);

if SHOW_PROFILE
    xx = (1:1000) / 1000;
    ff = zeros(1000);
    for i = 1:1000
        x_p = xx(i);
        ff(i) = objective_func([x_p 0.5], model);
    end

    figure(4);
    plot(xx, ff);
    ylim(Y_RANGE);
end


function y = objective_func(state, model)
    state = mod([1-state(1) state]', ones([3 1]));
    A = get_rcs_matrix(state, model);
    %y = double(norm(reshape(F, [], 1), 20));
    %y = 1/y;
    %y = dot(F, F.*[1 1 -1]');
    y = -double(det(A)^2);
end

function A = get_rcs_matrix(state, model)
    px = interp1(model.Ts, model.Xs, state);
    py = interp1(model.Ts, model.Ys, state);
    dx = (interp1(model.Ts, model.Xs, state+0.001) - px) / 0.001;
    dy = (interp1(model.Ts, model.Ys, state+0.001) - py) / 0.001;
    l = sqrt(dx.*dx+dy.*dy);
    dx = dx ./ l;
    dy = dy ./ l;
    torques = px.*dy - py.*dx;
    
    A = [dx dy torques]';
end