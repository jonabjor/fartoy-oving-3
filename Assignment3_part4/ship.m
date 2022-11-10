function [xdot,u, Qm] = ship(x,u,nu_c,tau_ext, Qm)
% [xdot,u] = ship(x,u,nu_c,tau_ext) returns the time derivative of the state vector: 
% x = [ u v r x y psi delta n ]' for a ship with L = 161 m where:
%
% u     = surge velocity, must be positive  (m/s)    
% v     = sway velocity                     (m/s)
% r     = yaw velocity                      (rad/s)
% x     = position in x-direction           (m)
% y     = position in y-direction           (m)
% psi   = yaw angle                         (rad)
% delta = actual rudder angle               (rad)
% n     = actual shaft velocity             (rpm)
% 
% The input vector is :
%
% u       = [ delta_c  n_c ]'  where
%
% delta_c = commanded rudder angle          (rad)
% n_c     = commanded shaft velocity        (rpm)
%
% The current vector is : 
%
% nu_c    = [ u_c v_c 0 ]'  where
%
% uc = current velocity in surge            (m/s)
% vc = current velocity in sway             (m/s)
%
% The external environmental force vector is denoted by tau_ext
%
% Author:    name
% Date:      date

% Check of input and state dimensions
if (length(x)~= 8),error('x-vector must have dimension 8 !');end
if (length(u)~= 2),error('u-vector must have dimension 2 !');end

% Dimensional states and input
delta_c = u(1); 
n_c     = u(2);

nu    = x(1:3);
eta   = x(4:6);
delta = x(7);
n     = x(8); 
h     = 0.1;

nu_r = nu - nu_c;
uc = nu_c(1);
vc = nu_c(2);

% ship parameters 
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

Bi = @(u_r,delta) [ (1-t_thr)  -u_r^2 * X_delta2 * delta
                        0      -u_r^2 * Y_delta
                        0      -u_r^2 * N_delta            ];
    
% state-dependent time-varying matrices
CRB = m * nu_r(3) * [ 0 -1 -xg 
                      1  0  0 
                      xg 0  0  ];

CA = [ 0 0 Yvdot*nu_r(2) + Yrdot*nu_r(3); 0 0 -Xudot*nu_r(1); -Yvdot*nu_r(2)-Yrdot*nu_r(3) Xudot*nu_r(1) 0];

% linear damping
T1 = 20; T2 = 20; T6 = 10;
Xu = -(m-Xudot)/T1;
Yv = -(m-Yvdot)/T2;
Nr = -(Iz-Nrdot)/T6;
D = -diag([Xu Yv Nr]);

% nonlinear surge damping
eps = 0.001;
CR = 0;
k = 0.1;
S = B*L + 2*T*(B+L);
v = 1e-6;
Rn = L / v * abs(nu_r(1));
Cf = 0.075 / (log10(Rn) - 2 + eps)^2 + CR;
Xns = -0.5*rho*S*(1+k)*Cf*abs(nu_r(1))*nu_r(1);

% nonlinear cross-flow drag
Cd_2d = Hoerner(B,T);
dx = L/10;
Ycf = 0; Ncf = 0;
for xL = -L/2:dx:L/2
Ucf = abs(nu_r(2) + xL * nu_r(3)) * (nu_r(2) + xL * nu_r(3));
Ycf = Ycf - 0.5 * rho * T * Cd_2d * Ucf * dx;
Ncf = Ncf - 0.5 * rho * T * Cd_2d * xL * Ucf * dx;
end
d = -[Xns Ycf Ncf]';


% propeller dynamics
PD = 1.5; AEAO = 0.65; z = 4;
[KT, KQ] = wageningen(Ja, PD, AEAO, z);
Qd = rho * Dia^5 * KQ * abs(n_c) * n_c;
[thr, Qm, n_dot] = propeller(n, Qd, Qm, Dia, rho, KT, KQ, h);


% ship dynamics
R = Rzyx(0,0,eta(3));
u = [ thr delta ]';
tau = Bi(nu_r(1),delta) * u;
nu_dot = [nu(3)*vc -nu(3)*uc 0]' + Minv * (tau_ext + tau - (CRB + CA + D) * nu_r - d);
eta_dot = R * nu;    

% Rudder saturation and dynamics (Sections 9.5.2)
if abs(delta_c) >= delta_max
    delta_c = sign(delta_c)*delta_max;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max
    delta_dot = sign(delta_dot)*Ddelta_max;
end    

xdot = [nu_dot' eta_dot' delta_dot n_dot]';
u = [delta_c n_c]';
end