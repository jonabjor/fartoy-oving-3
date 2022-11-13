%% exIntWindUp
% M-script for demonstration of integrator wind-up when the control law is
% saturated. The solution is to modify the integrator state by using
% anti-wind-up, see Section 6.5 in Beard & McLain (2012). The control
% objective is to regulate x to zero for a vanishing disturbance w.
%
% Continious-time system:  
%                            .     
%                            x  = a x + u + w
% Control law:                                      .
%                            u  = - Kp x - Ki z,    z = x
% 
% Definitions:               w = disturbance (step + white noise)
%                            u = control input
%                            x = state
%                            z = integral state
%
% Author:                    2020-09-04 Thor I. Fossen

%% USER INPUTS
h  = 0.1;                   % sample time (s)
N  = 2000;                  % number of samples

a  = -0.1;                  % model parameter
w0 =  0.15;                 % constant disturbance for t < t_switch
t_switch = 100;             

Kp  = 1;                    % PI controller gains
Ki  = 0.1;q
u_max = 0.1;                % saturation (u_max > w0 is OK, but test u_max < w0)
%u_max = 0.2;

table1 = zeros(N+1,4);      % table for simulation data
table2 = zeros(N+1,4);      % table for simulation data

%% FOR LOOP 1 (without anti-wind-up)
x = 0; z = 0;
for i = 1:N+1
    
   t = (i-1)*h;                             % time
   u = -Kp * x - Ki * z;                    % control law 
   
   if abs(u) > u_max                        % saturation
       u = sign(u) * u_max;          
   end

   w = w0 + 0.01 * randn(1);                % process disturbance
   if t > t_switch
       w = 0.01 * randn(1);
   end
   
   x_dot = a * x + u + w;                   % system model
   z_dot = x;                               % integral state
   
   table1(i,:) = [t x z u];                 % store data in table1
   
   % Euler integration   
   z = z + h * z_dot;   
   x = x + h * x_dot;               
   
end
 
%% FOR LOOP 2 (with anti-wind-up)
x = 0; z = 0;  
for i = 1:N+1
    
   t = (i-1)*h;                             % time
   u_unsat = -Kp * x - Ki * z;              % unsaturated control law 
   
   % saturated control law and integrator anti-wind-up
   if abs(u_unsat) > u_max
       u = sign(u_unsat) * u_max;           % saturation
       z = z - (h/Ki) * (u - u_unsat);      % anti-wind-up
   else
       u = u_unsat;                         % no saturation
   end

   w = w0 + 0.01 * randn(1);                % process disturbance
   if t > t_switch
       w = 0.01 * randn(1);
   end
   
   x_dot = a * x + u + w;                   % system model
   z_dot = x;                               % integral state
   
   table2(i,:) = [t x z u];                 % store data in table2
   
   % Euler integration 
   z = z + h * x;   
   x = x + h * x_dot;               
   
end

%% PLOTS 
t1   = table1(:,1);  
x1   = table1(:,2);  
z1   = table1(:,3); 
u1   = table1(:,4); 

t2   = table2(:,1);  
x2   = table2(:,2);  
z2   = table2(:,3); 
u2   = table2(:,4);  

figure(gcf)
subplot(311),plot(t1,x1,t2,x2,'linewidth',2),xlabel('time (s)'),title('state x'),grid
legend('without anti-wind-up','with anti-wind-up')
subplot(312),plot(t1,z1,t2,z2,'linewidth',2),xlabel('time (s)'),title('integral state z'),grid
legend('without anti-wind-up','with anti-wind-up')
subplot(313),plot(t1,u1,t2,u2,'linewidth',2),xlabel('time (s)'),title('control u'),grid
legend('without anti-wind-up','with anti-wind-up')

