% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.1;    % sampling time [s]
Ns = 64000; %10000;    % no. of samples

psi_ref = -110 * pi/180;  % desired yaw angle (rad)
U_ref   = 9; %7;            % desired surge speed (m/s)

% initial states
eta = [0 0 deg2rad(-0)]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;
x = [nu' eta' delta n]';
Qm = 0;

% Initial disturbance
% Current
V_c     = 1;   %1  % [m/s]  Problem 2: Turn off the current
beta_Vc = pi/4;   % Degrees
% Wind
V_w     = 10;
beta_Vw = 3*pi/4;
rho_a   = 1.247;
c_y     = 0.95;
c_n     = 0.15;
L_oa    = 161;
A_Lw    = 10*L_oa;

% Ref. model
a_d     = 0;
r_d     = 0;
psi_d   = psi_ref; 
psi_r   = 0;
psi_int = 0;
u_d     = U_ref;

% Propeller Revolution and Speed Control
Ja = 0; PD = 1.5; AEAO = 0.65; z = 4;
Dia = 3.3; m = 17.0677e6; T1 = 20; rho = 1025;
[KT, KQ] = wageningen(Ja, PD, AEAO, z);

Xudot   = -8.9830e5;
Xu      = -(m-Xudot)/T1;

% Guidance
j           = 2;      % next waypoint number
finished    = 0;
kappa       = 0.7;
y_int       = 0;
Delta       = 1000;
ship_length = 161;
R_switch    = 5 * ship_length;
addpath('WP_and_pathplotter')
WP          = load('WP.mat').WP;
wpt.pos.x   = WP(1, :);
wpt.pos.y   = WP(2, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,17);       % table of simulation data

for i=1:Ns+1
%     fprintf('(x, y) = (%.2f, %.2f) \n',x(4), x(5));
    t       = (i-1) * h;    % time (s)
    psi     = x(6);
    
    % current disturbance
    uc      = V_c * cos(beta_Vc - psi);
    vc      = V_c * sin(beta_Vc - psi);
    nu_c    = [ uc vc 0 ]';
    
    % wind disturbance, after 200 s
    if t >= 200
        u_w         = V_w * cos(beta_Vw - psi);
        v_w         = V_w * sin(beta_Vw - psi);
        V_rw        = sqrt((x(1) - u_w)^2 + (x(2) - v_w)^2);
        gamma_rw    = -atan2((x(2) - v_w), (x(1) - u_w));
        Ywind       = 0.5 * rho_a * V_rw^2 * c_y * sin(gamma_rw) * A_Lw;
        Nwind       = 0.5 * rho_a * V_rw^2 * c_n * sin(2*gamma_rw)...
                        * A_Lw * L_oa;
%         % New set-point, -20 degrees
%         psi_int     = 0;
%         psi_ref     = -20*pi/180;
%         psi_d       = psi_ref;
    else
        Ywind       = 0;
        Nwind       = 0;
    end
    tau_wind        = [0 Ywind Nwind]';
    
    % Sideslip and crab angle
    u_r         = x(1) - uc;
    v_r         = x(2) - vc;
    beta        = asin(v_r/sqrt(u_r^2 + v_r^2));
    beta_c      = asin(x(2)/ x(1));
    
    % guidance law
    % Computing chi_d
    x_now = x(4);                 
    y_now = x(5); 
    if not(finished) 
         [x_t, y_t, x_ref, y_ref, finished, j_updated] = target_wp(x_now, y_now, wpt, R_switch, j);
         j = j_updated;
    end
    if finished
        U_ref = 0;
    end

    
    % chi_d             = guidance(x_t, y_t, x_ref, y_ref, x_now, y_now, Delta);  
    [chi_d,y_int]       = guidance_ILOS(x_t, y_t, x_ref, y_ref, x_now, y_now, Delta, y_int, kappa, h);
    psi_d               = chi_d;
        
    % reference models
    [psi_r, r_d, a_d]   = ref_model(chi_d, psi_r, h, r_d, a_d); %ref_model(psi_d, psi_r, h, r_d, a_d); % Pre guidance law.
    %[psi_r, r_d, a_d]   = ref_model((chi_d - beta_c), psi_r, h, r_d, a_d);    % Compensating for crab angle
    
    % control law
    [delta_c, psi_int]  = PID(psi, psi_r, r_d, psi_int, h); % rudder angle command (rad)

%   When putting t as something else than time... t \in [0.05-0.2], put it 
%   to 0.05
    Td              = U_ref * Xu / (0.05 - 1);
    nnd             = Td / (rho * Dia^4 * KT);
    n_c             = sqrt(abs(nnd)) * sign(nnd);  % propeller speed (rps)
    
    % ship dynamics
    u               = [delta_c n_c]';
    [xdot,u,Qm]     = ship(x,u,nu_c,tau_wind,Qm);
    
    
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d beta beta_c, chi_d];
    
 
    % Euler integration
    x = euler2(xdot,x,h);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
u       = simdata(:,2);                 % m/s
v       = simdata(:,3);                 % m/s
r       = (180/pi) * simdata(:,4);      % deg/s
x       = simdata(:,5);                 % m
y       = simdata(:,6);                 % m
psi     = (180/pi) * simdata(:,7);      % deg
delta   = (180/pi) * simdata(:,8);      % deg
n       = 60 * simdata(:,9);            % rpm
delta_c = (180/pi) * simdata(:,10);     % deg
n_c     = 60 * simdata(:,11);           % rpm
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
beta    = (180/pi) * simdata(:,15);     % deg, sideslip
beta_c  = (180/pi) * simdata(:,16);     % deg, crab angle
chi_d   = (180/pi) * simdata(:,17)';    % deg, desired course angle

chi     = (psi + beta_c)';

figure(1)
figure(gcf)
subplot(211)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)'); xlabel('(m)'); ylabel('(m)'); 
subplot(212)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
% subplot(313)
% plot(t,r,t,r_d,'linewidth',2);
% title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

figure(3)
figure(gcf)
plot(t, beta, t, beta_c, 'linewidth', 2);
legend('Sideslip, \beta', 'Crab angle, \beta_c')
title('sideslip and crab angle');

figure(4)
figure(gcf)
plot(t, chi, t, chi_d, 'linewidth', 2);
legend('Course angle, \chi', 'Desired course angle, \chi_d')
title('Course angle and desired course angle');


% pathplotter
pathplotter(x, y)
