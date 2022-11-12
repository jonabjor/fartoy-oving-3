function chi_d = guidance(x_t, y_t, x_ref, y_ref, x, y, delta)

[~, ~, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y);

pi_p = atan2(y_t-y_ref, x_t-x_ref);


K_p = 1/delta;

% Formula 12.78 in Fossen
chi_d = ssa(pi_p - atan(K_p*y_e));

end

