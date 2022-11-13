% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.1;    % sampling time [s]
Ns = 10000;%64000; %10000;    % no. of samples

psi_ref = -110 * pi/180;  % desired yaw angle (rad)
U_ref   = 9; %7;            % desired surge speed (m/s)
fprintf("Getting K and T, using Nomoto:\n");
[K_nomoto, T_nomoto]  = nomoto(U_ref);

% initial states
eta = [0 0 deg2rad(-110)]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;
x = [nu' eta' delta n]'; % x_0
Qm = 0;

% Initial disturbance
% Current
%     Problem 2: Turn off the current
V_c     = 1;    %1;    % [m/s]
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
clear LOSchi
j           = 2;      % next waypoint number
finished    = 0;
Delta       = 1000;
ship_length = 161;
R_switch    = 5 * ship_length;
WP          = load('WP.mat').WP;
wpt.pos.x   = WP(1, :);
wpt.pos.y   = WP(2, :);


% Kalman filter setup
fprintf("Setting up for kalman filter, checking obervability (1c):\n");
[Ad, Bd, Cd, Dd, Ed]    = KF_setup(h);
x_hat                   = [0, 0, 0]';           % x_0
P_hat                   = eye(3);               %P0
Rd                      = deg2rad(0.5)^2;
q11                     = 0.01;
q22                     = 0.00001;
Qd                      = diag([q11, q22]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,22);       % table of simulation data

for i=1:Ns+1
%     fprintf('(x, y) = (%.2f, %.2f) \n',x(4), x(5));
    t       = (i-1) * h;    % time (s)
    psi     = x(6);
    psi_n   = x(6) + deg2rad(normrnd(0, 0.5));  % 2a) used to plot
    r_n     = x(3) + deg2rad(normrnd(0, 0.1));  % 2a) used to plot  
    
    
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

    
    chi_d               = guidance(x_t, y_t, x_ref, y_ref, x_now, y_now, Delta);
    psi_d               = chi_d;
        
    % reference models
    [psi_r, r_d, a_d]   = ref_model(chi_d, psi_r, h, r_d, a_d); %ref_model(psi_d, psi_r, h, r_d, a_d); % Pre guidance law.
%     [psi_r, r_d, a_d]   = ref_model((chi_d - beta_c), psi_r, h, r_d, a_d);    % Compensating for crab angle
    % control law
    [delta_c, psi_int]  = PID(psi, psi_r, r_d, psi_int, h, K_nomoto, T_nomoto); % rudder angle command (rad)

%   When putting t as something else than time... t \in [0.05-0.2], put it 
%   to 0.05
    Td              = U_ref * Xu / (0.05 - 1);
    nnd             = Td / (rho * Dia^4 * KT);
    n_c             = sqrt(abs(nnd)) * sign(nnd);  % propeller speed (rps)
    
    % ship dynamics
    u               = [delta_c n_c]';
    [xdot,u,Qm]     = ship(x,u,nu_c,tau_wind,Qm);
    
    
    
    
    
 
    % Euler integration
    x = euler2(xdot,x,h);    
    
    % Kalman filter
    %KF-gain
    K       = P_prd * Cd * inv(Cd' * P_prd * Cd + Rd);
    IKC     = eye(3) - K * Cd';
    
    % measuring yaw angle
    y       = x(6) + deg2rad(normrnd(0, 0.5));
    
    % state- and covariance correcting
    x_hat   = x_prd + K * (y - Cd' * x_prd);
%     x_hat
%     pause
    P_hat   = IKC * P_prd * IKC' + K * Rd * K';
    
    % state- and covariance prediction
    x_prd   = Ad * x_hat + Bd * delta_c;
    P_prd   = Ad * P_hat * Ad' + Ed * Qd * Ed';
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d beta beta_c, chi_d, psi_n, r_n, x_hat'];
    
    
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
psi_n   = (180/pi) * simdata(:,18)';    % deg, noisy yaw angle
r_n     = (180/pi) * simdata(:,19)';    % deg, noisy yaw rate
psi_est = (180/pi) * simdata(:,20)';    % deg, est. yaw angle
r_est   = (180/pi) * simdata(:,21)';    % deg, est. yaw rate
delta_e = (180/pi) * simdata(:,22)';    % deg, est. rudder bias

chi     = (psi + beta_c)';

% figure(1)
% figure(gcf)
% subplot(211)
% plot(y,x,'linewidth',2); axis('equal')
% title('North-East positions (m)'); xlabel('(m)'); ylabel('(m)'); 
% subplot(212)
% plot(t,psi,t,psi_d,'linewidth',2);
% title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
% % subplot(313)
% % plot(t,r,t,r_d,'linewidth',2);
% % title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');
% 
% figure(2)
% figure(gcf)
% subplot(311)
% plot(t,u,t,u_d,'linewidth',2);
% title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
% subplot(312)
% plot(t,n,t,n_c,'linewidth',2);
% title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
% subplot(313)
% plot(t,delta,t,delta_c,'linewidth',2);
% title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');
% 
% figure(3)
% figure(gcf)
% plot(t, beta, t, beta_c, 'linewidth', 2);
% legend('Sideslip', 'Crab angle');
% title('sideslip and crab angle');
% 
% figure(4)
% figure(gcf)
% plot(t, chi, t, chi_d, 'linewidth', 2);
% legend('Course angle', 'Desired course angle')
% title('Course angle and desired course angle');

figure(5)
figure(gcf)
subplot(211)
plot(t, psi_n, t, psi,'linewidth',2);
title('Noisy yaw angle and true yaw angel'); xlabel('time (s)'); ylabel('degree'); 
subplot(212)
plot(t,r_n,t,r,'linewidth',2);
title('Noisy yaw rate and true yaw rate'); xlabel('time (s)');

figure(6)
figure(gcf)
subplot(311)
plot(t, psi_est, t, psi,'linewidth',2);
title('Estimated yaw angle and true yaw angel'); xlabel('time (s)'); ylabel('degree'); 
subplot(312)
plot(t,r_est,t,r,'linewidth',2);
title('Estimated yaw rate and true yaw rate'); xlabel('time (s)'); ylabel('degree/second'); 
subplot(313)
plot(t,delta_e,'linewidth',2);
title('Rudder bias estimate'); xlabel('time (s)'); 


% pathplotter
pathplotter(x, y)
