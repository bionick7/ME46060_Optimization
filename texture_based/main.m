echo off;

% Display seed for reproducability
seed = randi(intmax("uint32"), 1, "uint32")
seed = 4237184565;
rng(seed);

PLOT_SEARCH_DIRECTIONS = false;

% Get position, normal and sdf maps
[position_map, normal_map] = get_position_normals("ih");
[~, ~, sdf_alpha] = imread(sprintf("../models/textures/%s/sdf.png", "ih"));
sdf = flip(sdf_alpha) * 2.0 - 1.0;

% Creates interpolants for position, normal ...
map_size = size(position_map);
map_size_xx = (0:map_size(1)-1) / (map_size(1) - 1);
map_size_yy = (0:map_size(2)-1) / (map_size(2) - 1);
position_interpolant = griddedInterpolant({map_size_xx map_size_yy}, position_map);
normal_interpolant = griddedInterpolant({map_size_xx, map_size_yy}, normal_map);

% ... and sdf maps
sdf_map_size = size(sdf);
sdf_map_size_xx = (0:sdf_map_size(1)-1) / (sdf_map_size(1) - 1);
sdf_map_size_yy = (0:sdf_map_size(2)-1) / (sdf_map_size(2) - 1);
sdf_interpolant = griddedInterpolant({sdf_map_size_xx, sdf_map_size_yy}, sdf);

% Maps structure bundles all the model parameters into one object
maps = struct(...
    'size', map_size,...
    'sdf_size', sdf_map_size,...
    'position_interpolant', position_interpolant,...
    'normal_interpolant', normal_interpolant,...
    'sdf_interpolant', sdf_interpolant...
    );

N_x = 9;  % Design variables

fun = @(x)objective_fun(x,maps);                                              % Pure objective function
fun_con = @(x)objective_fun(x,maps) + 1000 * max(constraints(x, maps), 0)^2;  % Objective function + constraints
% Initialize x (input) and y (output)
x = informed_initial_guess(maps);
y = 0;  % Set y to maximum possible value

[x, y, x_out, s_out, a_out] = powells_conjugate_directions(fun_con, x);

% Cross-sections plot
if PLOT_SEARCH_DIRECTIONS
    figure();
    hold on;
    evals = 1:9:9*5;
    optima_alpha = a_out(evals);
    optima_y = zeros(size(optima_alpha));
    counter = 1;
    for j = evals
        x0 = x_out(j,:)';
        s_i = s_out(j,:)';
        alphas = -0.25:0.001:0.25;
        yy = zeros([length(alphas) 1]);
        for i=1:length(alphas)
            yy(i) = fun_con(x0 + s_i*alphas(i));
        end
        plot(alphas, yy, 'DisplayName', sprintf("cycle %d", counter));
        ylim([-2, 1]);
        optima_y(counter) = fun_con(x0 + s_i * optima_alpha(counter));
        counter = counter + 1;
    end
    scatter(optima_alpha, optima_y, 'DisplayName', "Optima");
    xlabel("Profile along line")
    ylabel("f(x)")
    hold off;
    grid on;
end

for i = 1:100
    % Find starting point
    x0_loc = informed_initial_guess(maps);

    % Run optimiyation
    [xloc, yloc] = powells_conjugate_directions(fun_con, x0_loc);
    fprintf("%d => %f\n",i, yloc)
    if yloc < y
        % Print best x, y
        x = xloc
        y = yloc
    end
end

% Run fminsearch to refine solution
[x_precise, y_precise] = powells_conjugate_directions(fun, x, 0.05, 1e-5)

state = state_transform(x_precise);

% Display output
figure();
imshow(normal_map*.5 + .5, 'XData', [0 1], 'YData', [0 1]);
%imshow(single(sdf > 0), 'XData', [0 1], 'YData', [0 1]);
hold on;
scatter(state(:,1), state(:,2), 30, 'white', 'filled');
hold off;

%save_results(x, seed, maps);

function x = informed_initial_guess(maps)
    % Restricts domain in which x can appear to sample relevant 
    % domain more densly
    x = zeros([9 1]);
    x(1:3:7) = rand([3 1]) * 0.5;
    x(2:3:8) = rand([3 1]) * 0.5;
    x([5 8]) = x([5 8]) + 0.5;
    x(3:3:9) = rand([3 1]);

    % Regenerate if within constraints
    if constraints(x, maps) > 0
        x = informed_initial_guess(maps);
    end
end

function y = objective_fun(x, maps)
    % Finds the thruster matrix from function imput
    state = state_transform(x);    
    [rr, nn] = get_position_normal_at(state, maps);
    A = [cross(rr, nn), nn]';

    % Final objective function
    y = -double(norm([det(A), 0.1]));
end

function [cin, ceq] = constraints(x, maps)
    % Single inequality: samples sdf and returns inverse.
    state = state_transform(x);
    pt_adj = state(:, [2 1]);  % Flips x and y coordinates
    sdfs = squeeze(maps.sdf_interpolant(pt_adj));
    cin = max(-sdfs);
    ceq = [];
end

function plot_sampling_verification(position_map, normal_map, maps)
    % Verifies that the texture sampling is working correctly

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
