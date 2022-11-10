function [x_t, y_t, x_ref, y_ref, finished, i_updated] = target_wp(x, y, WP, R, i)

finished    = 0;
% Circle of acceptance
i_updated = i;

[~, max_i] = size(WP.pos.x);

x_ref   = WP.pos.x(i-1);
y_ref   = WP.pos.y(i-1);
x_t     = WP.pos.x(i);
y_t     = WP.pos.y(i);

% disp(x_t);
% disp(x);
% disp("-----------")
% pause

if (x_t-x)^2 + (y_t-y)^2 <= R^2
    disp(i)
    if max_i == i
        finished = 1;
    else
        x_ref = x_t;
        y_ref = y_t;

        i_updated = i_updated+1;
        x_t = WP.pos.x(i_updated);
        y_t = WP.pos.y(i_updated);  
    end
end

end