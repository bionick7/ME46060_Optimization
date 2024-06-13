%clc;
%clear all;
echo off;

PLOT_ANNEAHLING_OUTP = true;

seed = randi(intmax("uint32"), 1, "uint32")
%seed = 1046065145;
rng(seed);

[position_map, normal_map] = get_position_normals("ih");
[~, ~, sdf_alpha] = imread(sprintf("../models/textures/%s/sdf.png", "ih"));
sdf = flip(sdf_alpha) * 2.0 - 1.0;

map_size = size(position_map);
map_size_xx = (0:map_size(1)-1) / (map_size(1) - 1);
map_size_yy = (0:map_size(2)-1) / (map_size(2) - 1);
position_interpolant = griddedInterpolant({map_size_xx map_size_yy}, position_map);
normal_interpolant = griddedInterpolant({map_size_xx, map_size_yy}, normal_map);

sdf_map_size = size(sdf);
sdf_map_size_xx = (0:sdf_map_size(1)-1) / (sdf_map_size(1) - 1);
sdf_map_size_yy = (0:sdf_map_size(2)-1) / (sdf_map_size(2) - 1);
sdf_interpolant = griddedInterpolant({sdf_map_size_xx, sdf_map_size_yy}, sdf);

maps = struct(...
    'size', map_size,...
    'sdf_size', sdf_map_size,...
    'position_interpolant', position_interpolant,...
    'normal_interpolant', normal_interpolant,...
    'sdf_interpolant', sdf_interpolant,...
    'sphere_factor', 0.0...
    );

%tst_x = [ 0.4439 1.5020 3.8477 0.3586 0.9195 0.0326 13.4726 2.6893 1.2606 0.2512 1.1282 1.6366];
%tst_state = state_transform(tst_x)
%[tst_pos, tst_norm] = get_position_normal_at(tst_state, maps)
%save_results(tst_x, -1, maps)

%custom_x = [0 0.01 0; 0.25 0.01 0; 0 0.999 0; 0.25 0.999 0; 0.75 0.6 0.0; 0.25 0.6 0.0];
%% Best spherical
%custom_x = [52.4039 -54.5439 25.8002 -56.5655 -16.4605 5.6100 17.5807 61.1205 11.4135 6.9723 14.9902 -43.1983];
%custom_x = [0.6877 0.3709 0.8006 0.4998 0.3860 0.6231 0.5913 0.5616 0.7475 0.7572 0.9968 0.4998];

%plot_sampling_verification(position_map, normal_map, maps);

% Initialize the state randomly
N_x = 9;
x0 = rand([N_x 1]);
fun = @(x)objective_fun(x,maps);
fun_con = @(x)objective_fun(x,maps) + 100 * max(constraints(x, maps), 0)^2;
%[x, y, exit, output] = fminunc(fun, state, optimset('Display', 'iter'))
%[x, y] = patternsearch(fun_con, x0)
%[x, y] = fminunc(fun, x0)

%[x, y] = ga(fun_con, N_x)
%[x, y] = ga(fun, N_x, [], [], [], [], [], [], @(x)constraints(x,maps))
%[x, y] = patternsearch(fun_con, x0)
%[x, y] = particleswarm(fun_con, N_x)
%[x, y] = simulannealbnd(fun_con, custom_x, [], [])

%fun(custom_x)
%[x, y, exit_flag, outp, lambda, grad] = fmincon(fun, custom_x,...
%    [], [], [], [], [], [], @(x)constraints(x,maps),...
%    optimset('Display', 'iter'))
%return

%ms = MultiStart('StartPointsToRun', 'bounds-ineqs', 'UseParallel', true, 'Display', 'iter');
%gs = GlobalSearch;
%problem = createOptimProblem('fmincon','x0',x,'objective',fun,'nonlcon', @(x)constraints(x,maps));

y = 0;
x = zeros([N_x 1]);

%[xloc, yloc] = simulannealbnd(fun_con, x0_loc, [], [], optimoptions(@simulannealbnd, 'Display', 'iter'))
[xloc, yloc, step_info] = simulated_annealing(fun_con, x, 100, 1000);

if PLOT_ANNEAHLING_OUTP
    figure(1)
    plot(1:length(step_info), step_info(1,:));
    hold on;
    plot(1:length(step_info), step_info(2,:)/100);
    plot(1:length(step_info), step_info(3,:)/100);
    ylim([-1 1])
    hold off;
end
return

for i = 1:10
    x0_loc = informed_initial_guess();
    while constraints(x0_loc, maps) > 0
        x0_loc = informed_initial_guess();
    end
    %x0_loc = mod(x0_loc, 1);
    %[xloc, yloc] = fmincon(fun, x0_loc, [], [], [], [], [], [], @(x)constraints(x, maps), ...
    %               optimoptions(@fmincon, 'Display', 'off'));
    %[xloc, yloc] = simulannealbnd(fun_con, x0_loc, [], [], optimoptions(@simulannealbnd, 'Display', 'off'));
    [xloc, yloc] = simulated_annealing(fun_con, x0_loc, 100, 10000);
    [xloc, yloc] = fminunc(fun_con, xloc, optimoptions(@fminunc, 'Display', 'off'));
    fprintf("%d => %f\n",i, yloc)  % Single line
    if yloc < y
        x = xloc;
        y = yloc;
        yloc
    end
end

%[x, y] = run(ms, problem, 200)
%problem = createOptimProblem('fmincon','x0',x,'objective',fun,'nonlcon', @(x)constraints(x,maps));
%[x, y] = run(gs,problem)


%[x, y] = simulannealbnd(fun_con, x, [], [], optimoptions(@simulannealbnd, 'InitialTemperature', 500))

%constraints(x, maps)
    
% Best spherical solution found
%x = [-0.2706 -0.2345 0.1676 0.3751 0.0037 -0.4793 0.0520 0.3944 -0.3436 0.3853 -0.5253 -0.4422 -0.9877 -0.5198 -0.2573 -0.7053 0.0218 0.8265];
%y = -7.9841;
%x = custom_x;
%for sphere_factor = (10:-1:0) / 10
%    maps.sphere_factor = sphere_factor;
%    fun = @(x)objective_fun(x,maps);
%    fun_con = @(x)objective_fun(x,maps) + 100 * max(constraints(x, maps), 0)^2;
%    [x, y] = simulannealbnd(fun_con,x);
%    %[x, y] = fmincon(fun, custom_x, [], [], [], [], [], [], @(x)constraints(x,maps));
%    [sphere_factor, y]
%end

state = state_transform(x);

position_map_norm = position_map ./ repmat(vecnorm(position_map, 2, 3), [1 1 3]); % Normalization

figure();
imshow(position_map*.5 + .5, 'XData', [0 1], 'YData', [0 1]);
%imshow(single(sdf > 0), 'XData', [0 1], 'YData', [0 1]);
hold on;
scatter(state(:,1), state(:,2), 30, 'white', 'filled');
hold off;

%save_results(x, seed, maps);

function x = informed_initial_guess()
    x = zeros([9 1]);
    x(1:3:7) = rand([3 1]) * 0.5;
    x(2:3:8) = rand([3 1]) * 0.5;
    x([5 8]) = x([5 8]) + 0.5;
    x(3:3:9) = rand([3 1]);
end

function state = state_transform(x)
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

function A = get_rcs_matrix(x, maps)
    state = state_transform(x);
    %test_outputs = vertcat([eye(3) -eye(3)], zeros([3 6]));
    %test_outputs = vertcat(eye(3), zeros(3));
    %test_outputs = [0 0 1 0 0 0]';
    
    [rr, nn] = get_position_normal_at(state, maps);
    A = [cross(rr, nn), nn]';
end

function y = objective_fun(state, maps)
    A = get_rcs_matrix(state, maps);
    y = -double(det(A)^2);
    %test_inps = vertcat(eye(3), zeros(3));
    %y = -(1/double(norm(linsolve(A, test_inps), 20)));
end

function [cin, ceq] = constraints(x, maps)
    state = state_transform(x);
    pt_adj = state(:, [2 1]);  % Flips x and y coordinates
    sdfs = squeeze(maps.sdf_interpolant(pt_adj));
    cin = max(-sdfs);
    ceq = [];
end

function plot_sampling_verification(position_map, normal_map, maps)
    N = 10000;
    state = rand([N 3]);
    state(:, 3) = 0;
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
    normal_pos = single(imread(sprintf("../models/textures/%s/normals_pos.png", name))) / 65535.0;
    normal_neg = single(imread(sprintf("../models/textures/%s/normals_neg.png", name))) / 65535.0;
    normals = flip(normal_pos - normal_neg);
    positions = flip(single(imread(sprintf("../models/textures/%s/positions.png", name))) / 32767.5 - 1.0);
end