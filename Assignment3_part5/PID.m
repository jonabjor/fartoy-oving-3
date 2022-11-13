function [delta_c, psi_int] = PID(psi, psi_d, r_d, psi_int, h, K, T)
%PID Summary of this function goes here
% psi       : yaw, state parameter
% psi_d     : output from reference model, yaw
% r_d       : output from reference model, yaw-rate
% psi_int   : integral part from previous time-step
% h         : sample time
    sat     = @(x, limit) min(max(x, -limit), limit); % saturation function
    delta_max = deg2rad(30);
    
    % control law, PID controller [Kanskje bruke eksempel 15.7 aktivt?]
    w_b     = 0.06;         % [rad/s]
%     w_b     = 0.1;
    zeta    = 1;
    
    w_n     = 1 / (sqrt(1 - 2*zeta^2 + sqrt(4*zeta^4 - 4*zeta^2 + 2)))*w_b;
    
    % Nomoto model design
%     T       = 5.81*10;
%     K       = 0.00074*10;
%     
%     K       = 0.008229598;
%     T       = (192.307692308 + 6.877579092 - 13.961532797);

    m       = T/K;
    d       = 1/K;
    k       = 0;
    
    K_p     = (m)*w_n^2 - k;
    K_d     = 2*zeta*w_n*m - d;
    K_i     = w_n * K_p / 10;
%     disp(K_i)
    
    psi_int = psi_int + h * ssa(psi - psi_d);
    
    delta_unsat = -K_p * ssa(psi - psi_d) -K_d * r_d -K_i * psi_int; % unsaturated
    
    if abs(delta_unsat) > delta_max
        delta_c = sign(delta_unsat) * delta_max;
        psi_int = psi_int - (h/K_i) * (delta_c - delta_unsat); %Anti-windup
    else
        delta_c = delta_unsat;
    end
    
        
    
    
    
end

