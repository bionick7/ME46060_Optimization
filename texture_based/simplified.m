clc;
clear all;
close all;
echo off;

x_half = [0.50 0.48 0.44 0.38 0.38 0.38 0.32 0.32 0.44 0.47 0.50];
y_half = [0.00 0.00 0.11 0.34 0.52 0.70 0.87 0.94 0.94 1.00 1.00];
xs = [x_half(1:10) flip(1.0 - x_half)];
ys = [y_half(1:10) flip(y_half)];
t = (0:1:length(xs)-1) / (length(xs)-1);

tt = 0:0.001:1;

model = struct('Xs', xs, 'Ys', ys, 'Ts', t);
x0 = rand([1 2]);
[x, y] = fminunc(@(x)objective_func(x, model), x0)
%[x, y] = ga(@(x)objective_func(x, model), 2)
x_mod = mod([1-x(1) x]', ones([3 1]));

xx = spline(t, xs, tt);
yy = spline(t, ys, tt);

figure(1);
hold on;
plot(xs,ys,'rx');
plot(xx,yy,'r-');
plot(interp1(t, xs, x_mod), interp1(t, ys, x_mod),'bo');
axis equal;
hold off;

%grid = zeros(2);
%N = 100;
%for x_p = 1:N
%    for y_p = 1:N
%        grid(x_p, y_p) = objective_func([x_p y_p] / N, model);
%    end
%end
%figure(2);
%imshow(grid, 'XData', [0 1], 'YData', [0 1], 'DisplayRange', [-10 10]);
%hold on; 
%scatter(mod(x(1), 1), mod(x(2), 1), 'bo');
%hold off;


function y = objective_func(state, model)
    state = mod([1-state(1) state]', ones([3 1]));
    F = forces(state, model);
    y = double(norm(reshape(F, [], 1), 20));
    %y = dot(F, F.*[1 1 -1]');
end

function F = forces(state, model)
    px = spline(model.Ts, model.Xs, state);
    py = spline(model.Ts, model.Ys, state);
    dx = (spline(model.Ts, model.Xs, state+0.001) - px) / 0.001;
    dy = (spline(model.Ts, model.Ys, state+0.001) - py) / 0.001;
    dx = dx ./ sqrt(dx.*dx+dy.*dy);
    dy = dy ./ sqrt(dx.*dx+dy.*dy);

    torques = px.*dy - py.*dx;
    
    A = [dx dy torques]';
    %T = A * [0 -1 1]';
    if det(A) < 1e-5
        F = 10000;
    else
        F = linsolve(A, eye(3));
    end
    F = A * [1 1 1]';
end