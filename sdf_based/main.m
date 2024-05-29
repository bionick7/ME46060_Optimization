clc;
clear all;
echo off;

% Initialize the state randomly
x0 = rand([3 2]);

fun = @(x)objective_fun(x);
[x, y] = fminsearch(fun, x0)
draw_sdf_and_state(x);

%fun = @(x)objective_fun(x,maps);
%[x, y, exit, output] = fminsearch(fun, state, 'Display', 'iter')

%figure();
%imshow(position_map*.5 + .5);
%hold on;
%normalized_state = (x(:,[1 2])) .* repmat(maps.size([1 2]), [6 1]);
%scatter(normalized_state(:,1), normalized_state(:,2), 30, 'white', 'filled');
%hold off;

%save_results(x, maps);

function draw_sdf_and_state(state)
    im_size = 200;
    [X,Y] = meshgrid(-1:2/(im_size-1):1);
    Z = sdf(X,Y);
    img = 1. - (single(Z > 0) * 0.4);  % step function-ish
    figure();
    imshow(img, 'XData', [-1, 1], 'YData', [-1, 1]);
    hold on;
    scatter(state(:,1), state(:,2), 30, eye(3), 'filled');

    dir = get_sdf_gradient(state(:,1), state(:,2));
    dir = [dir(:,2) -dir(:,1)];  % Perpendicular to the gradient
    dir = dir ./ vecnorm(dir, 2, 2);
    F = thrusts(state)
    pos2 = state + dir * .05 .* F;
    scatter(pos2(:,1), pos2(:,2), 30, eye(3)*0.5, 'filled');

    hold off;
end

function y = objective_fun(state)
    F = thrusts(state);
    constraint = sum(power(sdf(state(:,1), state(:,2)),2)) * 10;
    y = norm(F, 10) + constraint;
end

function y = smooth_max(inps, p)
    y = log(sum(exp(inps*p)));
end

function F = thrusts(state)
    %test_outputs = vertcat([eye(3) -eye(3)], zeros([3 6]));
    %test_outputs = vertcat(eye(3), zeros(3));
    test_outputs = [0 0 1]';
    pos = state;
    dir = get_sdf_gradient(state(:,1), state(:,2));
    dir = [dir(:,2) -dir(:,1)];  % Perpendicular to the gradient
    dir = dir ./ vecnorm(dir, 2, 2);
    A = [pos, pos(:,1) .* dir(:,2) - pos(:,2) .* dir(:,1)]';
    F = linsolve(A, test_outputs);
end

function grad = get_sdf_gradient(x, y)
    epsilon = 0.0001;
    dFdx = (sdf(x + epsilon, y) - sdf(x - epsilon, y)) / (2*epsilon);
    dFdy = (sdf(x, y + epsilon) - sdf(x, y - epsilon)) / (2*epsilon);
    grad = [dFdx dFdy];
end
