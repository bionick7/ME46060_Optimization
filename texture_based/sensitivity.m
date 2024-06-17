% Relies on 'maps' to be set (e.g. by rudding main.m beforehand)

x = [ 0.2282 0.4817 0.2579 0.0944 0.9819 0.0958 0.7669 1.0255 0.1008]';
y = objective_fun(x, maps);

EPSILON = 5e-4;
unit_steps_9 = eye(9);
f_sensitivities = zeros([9 1]);
for i = 1:9
    unit_step = unit_steps_9(:,i);
    y_plus = objective_fun(x + unit_step * EPSILON, maps);
    y_minus = objective_fun(x - unit_step * EPSILON, maps);
    lin_sens = (y_plus - y_minus) / (2*EPSILON);
    f_sensitivities(i) = sum(x .* unit_step) / y * lin_sens;
end

%reshape(f_sensitivities, 3, 3)'
%objective_fun(x - f_sensitivities * 1e-6, maps)

EPSILON = 3e-4;
state = state_transform(x);
[rr, dd] = get_position_normal_at(state, maps);
A = [cross(rr, dd), dd]';
dA_drd = zeros(36);
for i=1:6
    r = rr(i,:);
    d = dd(i,:);
    Omega_d = [0 d(3) -d(2)
               -d(3) 0 d(1)
               d(2) -d(1) 0];
    Omega_r = [0 -r(3) r(2)
               r(3) 0 -r(1)
              -r(2) r(1) 0];
    dA_drd((i-1)*6+1:i*6-3, (i-1)*6+1:(i-1)*6+6) = [Omega_d Omega_r];
end
for i=1:6
    i_start = (i-1)*6+4;
    dA_drd(i_start:i_start+2, i_start:i_start+2) = eye(3);
end

rd_vec = reshape([rr dd]', [], 1);
%reshape(dA_drd * rd_vec, 6, 6)

ddet_dA = reshape(inv(A)' * det(A), 1, []);
%reshape(ddet_dA * dA_drd, 6, 6)

unit_steps = reshape(eye(18), 18, 6, 3);
f_sensitivities2 = -reshape(ddet_dA * dA_drd, 6, 6);

f_sensitivities2
% Account for normalization
for i = 1:6
    dif_d = f_sensitivities2(4:6,i);
    d = dd(i,:);
    transform = eye(3)*2 - (d' * 1./d);
    f_sensitivities2(4:6,i) = transform * dif_d;
end
f_sensitivities2

% Transform to logarithmic derivatives
f_sensitivities2 = f_sensitivities2 .* [rr dd]' / y

function y = objective_fun(x, maps)
    % Finds the thruster matrix from function imput
    state = state_transform(x);
    [rr, dd] = get_position_normal_at(state, maps);
    A = [cross(rr, dd), dd]';

    % Final objective function
    y = -det(A);
end

%rr2 = [
%    -0.35  0.8  0
%     0      -1  0
%    -0.35  0.8  0
%     0.35  0.8  0
%     0      -1  0
%     0.35  0.8  0
%];
%dd2 = [
%    0 0 1
%    1 0 0
%    0 0 1
%    0 1 0
%    0 0 1
%    0 1 0
%];