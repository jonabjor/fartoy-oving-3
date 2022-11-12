function [psi_r, r_r, a_r] = ref_model(psi_d, psi_r, h, r_r, a_r)
%REF_MODEL Summary of this function goes here
% phi_d     : Reference from e.g. operator
% phi_r     : Output reference from earlier step
% h         : sample time
% r_r       : yaw-rate from earlier step
% a_r       : yaw-rate-rate from earlier step :)
% NOTE: Some mess up with notation, but since the given out code had
% reference from operator as psi_d I let it be :=)

    sat     = @(x, limit) min(max(x, -limit), limit); % saturation function
    r_max   = 1;
    a_max   = 0.5;

%     w_ref       = 0.03;
    w_ref       = 0.5;
    zeta_ref    = 1;
    
    a_r         = a_r + (-(2*zeta_ref + 1)*w_ref*sat(a_r, a_max) ...
                    - (2*zeta_ref + 1)*w_ref^2 * sat(r_r, r_max) ...
                    + (psi_d - psi_r) * w_ref^3) *  h;
    r_r         = r_r + sat(a_r, a_max)*h;
    psi_r       = psi_r + sat(r_r, r_max)*h;
    
end

