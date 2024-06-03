
function sd = sdf(x, y, interpolant)
    N = size(x);
    x = x*2.0 - 1.0;
    y = y*2.0 - 1.0;

    % circle
    %l = sqrt(x.*x + y.*y);
    %sd = l - 0.6;

    % rect
    %dx = abs(x) - 0.6;
    %dy = abs(y) - 0.6;
    %max_d = elem_max(dx, dy);
    %dx = elem_max(dx, zeros(N));
    %dy = elem_max(dy, zeros(N));
    %sd = sqrt(dx.*dx + dy.*dy) + elem_min(max_d, zeros(N));

    % recle
    x = abs(x);
    y = abs(y);
    if y > x
        tmp = x;
        x = y;
        y = tmp;
    end
    a = x-y;
    b = x+y;
    c = (2.0*b-1.0)/3.0;
    h = sqrt(elem_max(a.*a+c.*c.*c,zeros(N)));
    u = power(elem_max(h-a,0.0),ones(N)/3.0);
    v = power(h+a, 1.0/3.0);
    t = (u-v)*0.5;
    w_x = -t + 0.75 - t.*t - x;
    w_y = t + 0.75 - t.*t - y;
    sd = sqrt(w_x.*w_x + w_y.*w_y) .* sign(a.*a * 0.5+b - 1.5);
end

function sd = sdf1(x, y, interpolant)
    p = [y x];
    sd = interpolant(p);
end

function sd = sdf2(x, y, UNUSED)
    assert(all(size(x) == size(y)));
    v = [0.1 0.1; 0.9 0.1; 0.1 0.9];
    M = length(v);
    dx = x-v(1,1);
    dy = y-v(1,2);
    d = dx.*dx+dy.*dy;
    s = ones(size(x));
    for i = 1:1:M
        j = mod(i-2, M) + 1;
        ex = v(j,1) - v(i,1);
        ey = v(j,2) - v(i,2);
        wx = x - v(i,1);
        wy = y - v(i,2);
        b_mult = elem_clamp((wx.*ex+wy.*ey)/(ex.*ex+ey.*ey), 0,1);
        bx = wx - ex .* b_mult;
        by = wy - ey .* b_mult;
        d = elem_min(d, bx.*bx+by.*by);

        c1 = y >= v(i,2);
        c2 = y < v(j, 2);
        c3 = ex.*wy > ey.*wx;
        s = s - 2 .* s .* double((c1 & c2 & c3) | (~c1 & ~c2 & ~c3));
    end
    sd = s.*sqrt(d) - 0.05;
end


function M = elem_clamp(A, lower, upper)
    M = elem_min(elem_max(A, lower), upper);
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