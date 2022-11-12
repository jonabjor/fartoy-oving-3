function [K, T]     = nomoto(ud)
    m = 17.0677e6;          % mass (kg)
    Iz = 2.1732e10;         % yaw moment of inertia (kg m^2)
    xg = -3.7;              % CG x-ccordinate (m)
    L = 161;                % length (m)
    B = 21.8;               % beam (m)
    T = 8.9;                % draft (m)
    KT = 0.7;               % propeller coefficient (-)
    Dia = 3.3;              % propeller diameter (m)
    rho = 1025;             % density of water (m/s^3)

    Ja = 0;

    % rudder limitations
    delta_max  = 40 * pi/180;        % max rudder angle      (rad)
    Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

    % added mass matrix
    Xudot = -8.9830e5;
    Yvdot = -5.1996e6;
    Yrdot =  9.3677e5;
    Nvdot =  Yrdot;
    Nrdot = -2.4283e10;
    MA = -[ Xudot 0 0; 0 Yvdot Yrdot; 0 Nvdot Nrdot ];

    % rigid-body mass matrix
    MRB = [ m 0    0 
            0 m    m*xg
            0 m*xg Iz ];
    Minv = inv(MRB + MA);

    % linear damping
    T1 = 20; T2 = 20; T6 = 10;
    Xu = -(m-Xudot)/T1;
    Yv = -(m-Yvdot)/T2;
    Nr = -(Iz-Nrdot)/T6;
    D = -diag([Xu Yv Nr]);


    % rudder coefficients
    b = 2;
    AR = 8;
    CB = 0.8;

    lambda = b^2 / AR;
    tR = 0.45 - 0.28*CB;
    CN = 6.13*lambda / (lambda + 2.25);
    aH = 0.75;
    xH = -0.4 * L;
    xR = -0.5 * L;

    % input matrix
    t_thr = 0.05;                                        % thrust deduction number
    X_delta2 = 0.5 * (1 - tR) * rho * AR * CN;           % rudder coefficients (Section 9.5)
    Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
    N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

    % ud = 7;
    % ud = 9;

    % Current
    V_c     = 1;    % [m/s]
    beta_Vc = pi/4;   % Degrees
    psi     = 10* pi/180;

    uc      = V_c * cos(beta_Vc - psi);
    vc      = V_c * sin(beta_Vc - psi);

    CAstar  = [0,   0,  -Yvdot * vc;
               0,   0,  -Xudot*(ud - uc);
               Yvdot * vc, (-Yvdot * ud + Xudot * (ud - uc)), -Yrdot*ud]; 
    CRBstar = [0 0 0;
               0 0 m*ud;
               0 0 m*xg*ud];



    A = -Minv(2:3, 2:3)*(CRBstar(2:3, 2:3)+CAstar(2:3, 2:3)+D(2:3, 2:3));
    B = Minv(2:3, 2:3)*(2*ud*[-Y_delta; -N_delta]);
    C = [0 1];
    Du = 0;

    [b, a]  = ss2tf(A, B, C, Du);

    sys     = tf(b, a);

    % get that num --> 0.0001111 s + 8.561e-06
    % get that den --> s^2 + 0.1506 s - 0.001831

    T1T2    = -1 ./ pole(sys);
    T1      = T1T2(1);
    T2      = T1T2(2);
    T3      = -1 / zero(sys);

    T       = T1 + T2 - T3;
    K       = b(end)/a(end);

    fprintf("K: \t%d\nT: \t%d\n", K, T);
end





% Assignment3_part3
% 
% 
% 
% K   = 0.008229598;
% 
% roots([1330.140994945, 200.319233839, 1])
% 
% T   = 192.307692308 + 6.877579092 - 13.961532797;