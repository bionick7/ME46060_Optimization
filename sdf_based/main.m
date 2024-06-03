clc;
clear all;
echo off;

COM = [0.5 0.7];

% Initialize the state randomly
x0 = rand([6 1]);
%x0 = 0.1 + [0.5 0; 0 1; 1 1] * 0.8;  % Initial guess (p = 100)
%x0 = [0.4531 0.1224; 0.1194 0.7660; 0.6730 0.8767];

sdf_image = single(imread("../models/sdfs/spaceshuttle_contour.png.sdf.png")) / 128.0 - 1.0;
img_size = size(sdf_image);
img_size_xx = (0:img_size(1)-1) / (img_size(1) - 1);
img_size_yy = (0:img_size(2)-1) / (img_size(2) - 1);
global interpolant
interpolant = griddedInterpolant({img_size_xx img_size_yy}, sdf_image);

% [0.1 0.1; 0.9 0.1; 0.1 0.9]

% -x1 < -0.1
% -y1 < -0.1
% x1 + y1 < 1

% A = [-1  0;
%       0 -1;
%       1  1; 
%        ...]

% b = [-0.1 -0.1 1]'

options = optimoptions('fmincon', 'Algorithm', 'interior-point',...
                       'Display','iter-detailed', 'EnableFeasibilityMode', false...
                        );
                        %'HessianApproximation', 'finite-difference',...
                        %'SubproblemAlgorithm', 'cg'...
                        %'FiniteDifferenceStepSize', 0.01,...
                        %'FiniteDifferenceType', 'central',...

%[x, y, exitflag, opts] = fminsearch(@objective_fun, x0)
%[x, y] = ga(@objective_fun, 6)

A_sub = [-1 0; 0 -1; 1 1];
A = [
    A_sub zeros([3 2]) zeros([3 2]);
    zeros([3 2]) A_sub zeros([3 2]);
    zeros([3 2]) zeros([3 2]) A_sub;
];
b = repmat([-.1 -.1 1]', [1 3]);

[x, y, exitflag, out, lambda, grad, hessian] = ...
    fmincon(@objective_fun, x0, A, b, [], [],...
            ones([6 1]) * 0.01, ones([6 1]) * 0.99)

[c, ceq] = constraint_function(x)
x = reshape(x, [3 2]);

sd = sdf(x(:,1), x(:,2), interpolant)
draw_sdf_and_state(sdf_image, x);

%fun = @(x)objective_fun(x,maps);
%[x, y, exit, output] = fminsearch(fun, state, 'Display', 'iter')

%figure();
%imshow(position_map*.5 + .5);
%hold on;
%normalized_state = (x(:,[1 2])) .* repmat(maps.size([1 2]), [6 1]);
%scatter(normalized_state(:,1), normalized_state(:,2), 30, 'white', 'filled');
%hold off;

%save_results(x, maps);

function res = sdf_repr(sd)
    res = sin(sd*100)*.2+.5;
    res = res + single(sd < 0) * 0.3;
end

function draw_sdf_and_state(sd, state)
    %global interpolant;
    %im_size = 500;
    %[X,Y] = meshgrid(0:2/(im_size-1):1);
    %sd = sdf(X,Y, interpolant);

    img = sdf_repr(sd);  % step function-ish
    figure();
    imshow(img, 'XData', [0, 1], 'YData', [0, 1]);
    hold on;
    scatter(state(:,1), state(:,2), 30, eye(3), 'filled');

    dir = get_sdf_gradient(state(:,1), state(:,2));
    dir = [dir(:,2) -dir(:,1)];  % Perpendicular to the gradient
    dir = dir ./ vecnorm(dir, 2, 2);
    F = thrusts(state)
    pos_x = state(:,1) - 0.5;
    pos_y = state(:,2) - 0.5;
    A = [dir, pos_x .* dir(:,2) - pos_y .* dir(:,1)]';
    A * F
    pos2 = state + dir * .05 .* F;
    scatter(pos2(:,1), pos2(:,2), 30, eye(3)*0.5, 'filled');
    scatter(0.5, 0.7, 40, 'yellow', 'filled');

    % Verify Alignment
    %xx = rand([10000 2]);
    %sd2 = sdf(xx(:,1), xx(:,2), interpolant);
    %scatter(xx(:,1), xx(:,2), 30, sdf_repr(sd2), 'filled');

    hold off;
end

function [c, ceq] = constraint_function(state)
    state = reshape(state, [3 2]);
    global interpolant;
    sd = sdf(state(:,1), state(:,2), interpolant);
    ceq = double(sum(sd.*sd));
    c = [];
end

function y = objective_fun(state)
    state = reshape(state, [3 2]);
    F = thrusts(state);
    global interpolant;
    sd = sdf(state(:,1), state(:,2), interpolant);
    %y = double(norm(reshape(F, 1, []), 20) + sum(sd.*sd) * 10000);
    %y = double(norm(reshape(F, 1, []), 20));
    y = double(norm(state));
end

function F = thrusts(state)
    test_outputs = [0 0 1]';
    %test_outputs = eye(3);
    dir = get_sdf_gradient(state(:,1), state(:,2));
    dir = [dir(:,2) -dir(:,1)];  % Perpendicular to the gradient
    dir = dir ./ vecnorm(dir, 2, 2);
    pos_x = state(:,1) - 0.5;
    pos_y = state(:,2) - 0.7;
    A = [dir, pos_x .* dir(:,2) - pos_y .* dir(:,1)]';
    F = linsolve(A, test_outputs);
end

function grad = get_sdf_gradient(x, y)
    global interpolant;
    epsilon = 4.0/500.0;
    dFdx = (sdf(x + epsilon, y, interpolant) - sdf(x - epsilon, y, interpolant)) / (2*epsilon);
    dFdy = (sdf(x, y + epsilon, interpolant) - sdf(x, y - epsilon, interpolant)) / (2*epsilon);
    grad = [dFdx dFdy];
end
