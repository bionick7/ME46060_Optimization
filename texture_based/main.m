clc;
clear all;
echo off;

seed = randi(intmax("uint32"), 1, "uint32")
rng(seed);

[position_map, normal_map] = get_position_normals("test_obj");
map_size = size(position_map);
map_size_xx = (0:map_size(1)-1) / (map_size(1) - 1);
map_size_yy = (0:map_size(2)-1) / (map_size(2) - 1);
position_interpolant = griddedInterpolant({map_size_xx map_size_yy}, position_map);
normal_interpolant = griddedInterpolant({map_size_xx, map_size_yy}, normal_map);

maps = struct('size', map_size, 'position_interpolant', position_interpolant, 'normal_interpolant', normal_interpolant);

custom_x = [0 0.01 0; 0.25 0.01 0; 0 0.999 0; 0.25 0.999 0; 0.75 0.6 0.0; 0.25 0.6 0.0];

%plot_sampling_verification(position_map, normal_map, maps);

% Initialize the state randomly
x0 = rand([18 1]);
fun = @(x)objective_fun(x,maps);
%[x, y, exit, output] = fminunc(fun, state, optimset('Display', 'iter'))
%[x, y] = simulannealbnd(fun,x0, zeros(18), ones(18), optimoptions('simulannealbnd', 'InitialTemperature', 500))
%[x, y] = patternsearch(fun, x0)
%[x, y] = fminunc(fun, custom_x)
[x, y] = ga(fun, 18)

%gs = GlobalSearch;
%problem = createOptimProblem('fmincon','x0',x0,...
%   'objective',fun,'lb',zeros([18 1]),'ub',ones([18 1]));
%x = run(gs,problem)

%x = custom_x;

x = mod(reshape(x, [6 3]), ones([6 3]));

figure();
imshow(position_map*.5 + .5, 'XData', [0 1], 'YData', [0 1]);
hold on;
scatter(x(:,1), x(:,2), 30, 'white', 'filled');
hold off;

save_results(x, maps);

% Store the results in a json file, so that blender can present it again (for sanity check)
function save_results(x, maps)
    [rr, nn] = get_position_normal_at(x, maps);
    y = objective_fun(x, maps);
    F = thrust(x, maps)
    A = [cross(rr, nn), nn]';
    det(A)
    result = struct('state', x, 'y', y, 'positions', rr, 'normals', nn, 'F', F);
    encoded = jsonencode(result);
    
    fid = fopen('../models/outp.json','w');
    fprintf(fid,'%s',encoded);
    fclose(fid);
end

function F = thrust(state, maps)
    state = mod(reshape(state, [6 3]), ones([6 3]));
    %test_outputs = vertcat([eye(3) -eye(3)], zeros([3 6]));
    test_outputs = vertcat(eye(3), zeros(3));
    %test_outputs = [1 0 0 0 0 0]';
    
    [rr, nn] = get_position_normal_at(state, maps);
    A = [cross(rr, nn), nn]';
    F = linsolve(A, test_outputs);
end

function y = objective_fun(state, maps)
    F = thrust(state, maps);
    y = double(norm(reshape(F, [], 1), 20));
end

function plot_sampling_verification(position_map, normal_map, maps)
    N = 10000;
    state = rand([N 3]);
    [positions_pts, normal_pts] = get_position_normal_at(state, maps);
        
    % Positions
    figure();
    imshow(position_map*0.5 + 0.5, 'XData', [0 1], 'YData', [0 1]);
    hold on;
    scatter(state(:,1), state(:,2), 30, positions_pts*0.5 + 0.5, 'filled');
    hold off;

    % Normals
    figure();
    imshow(normal_map*0.5 + 0.5, 'XData', [0 1], 'YData', [0 1]);
    hold on;
    scatter(state(:,1), state(:,2), 30, normal_pts*0.5 + 0.5, 'filled');
    hold off;
end 

function [positions, normals] = get_position_normals(name)
    % Returns the normalized positions and normals from the files
    normal_pos = single(imread(sprintf("../models/textures/%s/normals_pos.png", name))) / 255.0;
    normal_neg = single(imread(sprintf("../models/textures/%s/normals_neg.png", name))) / 255.0;
    normals = flip(normal_pos - normal_neg);
    positions = flip(single(imread(sprintf("../models/textures/%s/position.png", name))) / 128.0 - 1.0);
end