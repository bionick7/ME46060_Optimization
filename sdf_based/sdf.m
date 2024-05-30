function sd = sdf(x, y, interpolant)
    N = size(x);
    x = x*2.0 - 1.0;
    y = y*2.0 - 1.0;

    % circle
    %l = sqrt(x.*x + y.*y);
    %sd = l - 0.6;

    % rect
    dx = abs(x) - 0.6;
    dy = abs(y) - 0.6;
    max_d = elem_max(dx, dy);
    dx = elem_max(dx, zeros(N));
    dy = elem_max(dy, zeros(N));
    sd = sqrt(dx.*dx + dy.*dy) + elem_min(max_d, zeros(N));

    % recle
    %x = abs(x);
    %y = abs(y);
    %if y > x
    %    tmp = x;
    %    x = y;
    %    y = tmp;
    %end
    %a = x-y;
    %b = x+y;
    %c = (2.0*b-1.0)/3.0;
    %h = sqrt(elem_max(a.*a+c.*c.*c,zeros(N)));
    %u = power(elem_max(h-a,0.0),ones(N)/3.0);
    %v = power(h+a, 1.0/3.0);
    %t = (u-v)*0.5;
    %w_x = -t + 0.75 - t.*t - x;
    %w_y = t + 0.75 - t.*t - y;
    %sd = sqrt(w_x.*w_x + w_y.*w_y) .* sign(a.*a * 0.5+b - 1.5);
end

function sd = sdf_a(x, y, interpolant)
    p = [y x];
    sd = interpolant(p);
end

function M = elem_max(A, B)
    % A and B are row-vectors
    log = single(A > B);
    M = log .* A + (1 - log) .* B;
end

function M = elem_min(A, B)
    % A and B are row-vectors
    log = single(A < B);
    M = log .* A + (1 - log) .* B;
end