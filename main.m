clc
clear all
close all

set(0,'language','en')
% set(0,'language','ja')
set(0,'defaultAxesFontSize',13,'defaultTextFontSize',13)
set(0,'defaultAxesFontName','Helvetica', 'defaultTextFontName','Helvetica')
set(0, 'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2) 




%% set parameters

% -- simulation time-- 
t_end = 30;
dt = 0.001;
tspan = 0:dt:t_end;

% -- initial state values --
v_0 = 3;
w_0 = 0;
x_0 = [v_0; w_0];

% -- initial estimate --
v_0_est = v_0;
w_0_est = w_0;
mi_0_est = 1;
iz_0_est = 0.1;
z_0 = [v_0_est; w_0_est; 1 / mi_0_est; 1 / iz_0_est];

% -- physical parameters
m = 20;      % mass (variable)
b = 0.05;   % damping
iz = 2;   % inertia (variable)
l = 0.3;    % right-left wheel length / 2
param = [m, b, iz, l];


% -- dimension --
dim_x = 2;
dim_z = 4;
dim_u = 2;
dim_y = 2;


% -- initial covariance matrix in EKF --
% * variance of mi_est & iz_est may be large
P_0 = zeros(dim_z, dim_z);
P_0(3,3) = 100;
P_0(4,4) = 100;

% -- process covariance --
Q = zeros(dim_z, dim_z);
Q(1,1) = 0.01;  % velocist
Q(2,2) = 0.01;  % ang-vel
Q(3,3) = 0;     % m_inv
Q(4,4) = 0;     % iz_inv

% -- measurement covariance --
R = [1, 0; 0, 1];

% -- measurement noise property --
noise_mean = 0;
noise_cov = 0.1;

% -- controller parameters --
Qc = diag([1,10000,1000]);
Rc = diag([1,1]);
K_ini = [0, 1, 1];

% -- others --
save_interval = 10;
display_interval = 1000;
lqr_recalc_interval = 1000;



%% reference velocity & ang-vel

v_ref_mat = 3 + sat(4 * sin(tspan), 1, -1);
w_ref_mat = 5 * (sin(tspan) - sat(sin(tspan), 0.5, -0.5));

%% Simulation


% for save data
l_size = fix(length(tspan) / save_interval);
X = zeros(l_size, dim_x);
Y = zeros(l_size, dim_y);
U = zeros(l_size, dim_u);
Z = zeros(l_size, dim_z);
T = zeros(l_size, 1);
P_mat = zeros(l_size, dim_z, dim_z);

% initialization
x_next = x_0;
z_pred = z_0;
P_pred = P_0;
u_prev = [0; 0];
K_lqr = K_ini;
v_e_I = 0;

j = 1;
for i = 1:length(tspan)
    
    t = tspan(i);
    x_curr = x_next;
    
    % linearize at predicted state
    [F, H] = linearized_dynamics_func(z_pred, u_prev, dt, param);
    
    % output
    y = output_func(x_curr) + (noise_mean + noise_cov * randn);
    
    % == EKF: update process ==
    error = y - output_func(z_pred);
    K_opt = P_pred * H' / (H * P_pred * H' + R);
    z_curr = z_pred + K_opt * error;
    P_curr = (eye(dim_z) - K_opt * H) * P_pred;

    % == Adaptive LQR controller ==
    v_ref = v_ref_mat(i);
    w_ref = w_ref_mat(i);
    v_e = z_curr(1) - v_ref;
    w_e = z_curr(2) - w_ref;
    v_e_I = v_e_I + v_e * dt; % integral term
    servo_error = [v_e_I; v_e; w_e];
    if rem(i, lqr_recalc_interval) == 0
        [A, B] = continuous_linear_dynamics_func(z_curr, param);
        K_lqr = lqr(A, B, Qc, Rc);
    end
    % feedfoward input
    u_ff = [b * v_ref; b * v_ref];
    % feedback input
    u_fb = -K_lqr * servo_error;
    u = u_ff + u_fb;
    % zero-order hold
    u_ZOH = u;
    
    % == Continuous dynamics model ==
    % calc next state variables by 4th-order Runge Kutta
    k1 = continuous_dynamics_func(x_curr, u_ZOH, param);   
    k2 = continuous_dynamics_func(x_curr + k1 * dt / 2, u_ZOH, param);   
    k3 = continuous_dynamics_func(x_curr + k2 * dt / 2, u_ZOH, param);   
    k4 = continuous_dynamics_func(x_curr + dt * k3, u_ZOH, param);   
    x_next = x_curr + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6 ;
      
    % == EKF: discritized model ==
    % prediction
    z_pred = discritized_dynamics_func(z_curr, u, dt, param);
    P_pred = F * P_curr * F' + Q;
    
    
    u_prev = u;
    
    % == save data ==
    if rem(i, save_interval) == 0
        X(j,:) = x_curr';   
        Y(j,:) = y';    
        Z(j,:) = z_curr';
        U(j,:) = u';
        P_mat(j,:,:) = P_curr;
        T(j) = t;
        j = j + 1;
    end
    
    % == display ==
    if rem(i, display_interval) == 0
        fprintf('calc time = %3.3f, (%2.2f %%)\n', t, 100 * t / t_end);
    end
end
    

%% plot result
    
figure(1);
subplot(2,1,1)
plot(tspan, v_ref_mat, T, X(:,1), T, Y(:,1), T, Z(:,1));
legend('ref vel', 'real vel', 'measured vel','estimated vel');
grid on;
subplot(2,1,2);
plot(tspan, w_ref_mat, T, X(:,2), T, Y(:,2), T, Z(:,2))
legend('ref ang-vel','real ang-vel', 'measured ang-vel','estimated ang-vel');
grid on;
    

figure(2)
subplot(2,1,1)
plot(T, Z(:,3).^(1));grid on
legend('estimated inv mass'); hold on
plot(T, (Z(:,3) + sqrt(P_mat(:,3,3))).^(1), 'r--');
plot(T, (Z(:,3) - sqrt(P_mat(:,3,3))).^(1), 'r--');
subplot(2,1,2)
plot(T, Z(:,4).^(1));grid on
legend('estimated inv inertia'); hold on
plot(T, (Z(:,4) + sqrt(P_mat(:,4,4))).^(1), 'r--');
plot(T, (Z(:,4) - sqrt(P_mat(:,4,4))).^(1), 'r--');

    
    
    
    
   

    
    
    
    
    
    
    
    
    