
function [Ad, Bd, Cd, Dd, Ed]   = KF_setup(h)
    % From part 2, Nomoto calculation
    K       = 0.008229598;
    T       = (192.307692308 + 6.877579092 - 13.961532797)/10;

    % Using Nomoto model, with subtracted rudder bias at the input
    A       = [0, 1, 0; 0, -1/K, -K/T; 0, 0, 0;];
    B       = [0, K/T, 0]';
    C       = [1, 0, 0]';           % Measuring yaw angle
    D       = 0;                    % Zero disturbances
    E       = [0, 0; 1, 0; 0, 1;];  % Model noise.
    
    % First order linearization, Following (13.60 - 13.65 in [Fossen 2021])
    Ad      = eye(3) + h * A;
    Bd      = h * B;
    Cd      = C;
    Dd      = D;
    Ed      = h * E;
    
    % Observability
    fprintf("Observability [1-Yes/0-No]: \t %d\n", ...
        isequal(rank(obsv(Ad, Cd')), 3));
    
end