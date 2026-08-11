%% 8.6 Observer-Based SMC Comparison
% Time histories, regulation error, and SMC parameter sweep

clear; clc; close all;

%% =========================================================
% Model Parameters
% ==========================================================

M = 1.0;
Cv = 0.3;
Kv = 2.25;

phi_p_0 = 1.0;
T_0     = 0.0;
Ve_0    = 1.0;

alpha_B = 0.20;
alpha_E = 0.12;

gamma_Phi = 0.15;
k_Phi    = 0.04;
k_E      = 0.03;

C_th  = 1.0;
h     = 0.25;
tau_e = 0.10;

%% =========================================================
% Simulation Parameters
% ==========================================================

tspan = 0:0.01:10;

%% =========================================================
% Linearized Perturbation Model
% ==========================================================

A = [
    0, 1, 0, 0, 0;
   -Kv/M, -Cv/M, 2*alpha_B*phi_p_0/M, 0, 2*alpha_E*Ve_0/M;
    0, 0, -gamma_Phi, 0, 0;
    0, 0, 2*k_Phi*phi_p_0/C_th, -h/C_th, 2*k_E*Ve_0/C_th;
    0, 0, 0, 0, -1/tau_e
];

B = [0; 0; 0; 0; 1/tau_e];

C = [
    1 0 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1
];

%% =========================================================
% LQG Controller and Common LQE
% ==========================================================

Q = diag([100 10 5 5 1]);
R = 1;

K_lqr = lqr(A,B,Q,R);

W = 0.01*eye(5);
V = 0.01*eye(4);

% Common observer used by both LQG and SMC
L = lqe(A,eye(5),C,W,V);

%% =========================================================
% Representative Initial Condition
% ==========================================================

dx_0 = [0.20; 0.00; 0.05; 0.03; 0.10];

% Absolute nonlinear-state initial condition
x_0_nl = [
    dx_0(1);
    dx_0(2);
    phi_p_0 + dx_0(3);
    T_0     + dx_0(4);
    Ve_0    + dx_0(5)
];

% LQG starts with zero state estimate
z_0_lqg = [dx_0; zeros(5,1)];

%% =========================================================
% Uncontrolled Nonlinear Response
% ==========================================================

f_nl = @(t,x) nonlinear_dynamics( ...
    t,x,M,Cv,Kv,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,Ve_0);

[t_nl,x_nl] = ode45(f_nl,tspan,x_0_nl);

%% =========================================================
% LQG Controlled Response
% ==========================================================

f_lqg = @(t,z) lqg_dynamics( ...
    t,z,A,B,C,K_lqr,L);

[t_lqg,z_lqg] = ode45(f_lqg,tspan,z_0_lqg);

x_lqg    = z_lqg(:,1:5);
xhat_lqg = z_lqg(:,6:10);

%% =========================================================
% Observer-Based SMC Parameters
% ==========================================================

eta_s     = 2.0;
Delta_max = 248.7;
k_s       = 1.05*Delta_max;
phi_s     = 0.05;

u_min = 0.0;
u_max = 3.0;

sat = @(s) max(-1,min(1,s));

%% =========================================================
% Observer-Based SMC Response
% ==========================================================

% Initial observer estimate in perturbation coordinates
xhat_0_smc = [
    dx_0(1);
    0.0;
    dx_0(3);
    dx_0(4);
    dx_0(5)
];

z_0_smc = [x_0_nl; xhat_0_smc];

f_smc = @(t,z) smc_observer_dynamics( ...
    t,z,M,Cv,Kv,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,Ve_0, ...
    eta_s,k_s,phi_s,u_min,u_max, ...
    sat,A,B,C,L);

[t_smc,z_smc] = ode45(f_smc,tspan,z_0_smc);

x_smc    = z_smc(:,1:5);
xhat_smc = z_smc(:,6:10);

%% =========================================================
% Figure 13: Stress-Coordinate Time History
% ==========================================================

figure;

plot(t_nl,x_nl(:,1),'r','LineWidth',2);
hold on;

plot(t_lqg,x_lqg(:,1),'b','LineWidth',2);
plot(t_smc,x_smc(:,1),'k--','LineWidth',2);

yline(1e-2,'--','LineWidth',1);
yline(-1e-2,'--','LineWidth',1);

grid on;

xlabel('$t$ (s)','Interpreter','latex');
ylabel('$q(t)$','Interpreter','latex');

title('Dominant Stress Coordinate Regulation');

legend( ...
    'Nonlinear Uncontrolled', ...
    'LQG Controlled', ...
    'Observer-Based SMC Controlled', ...
    '$\pm 10^{-2}$ tolerance', ...
    'Location','best', ...
    'Interpreter','latex');

%% =========================================================
% Figure 14: Absolute Regulation Error
% ==========================================================

q_ref = 0;

e_lqg = abs(x_lqg(:,1) - q_ref);
e_smc = abs(x_smc(:,1) - q_ref);

figure;

semilogy(t_lqg,e_lqg,'b','LineWidth',2);
hold on;

semilogy(t_smc,e_smc,'k--','LineWidth',2);
yline(1e-2,'r--','LineWidth',0.5);

grid on;

xlabel('$t$ (s)','Interpreter','latex');
ylabel('$|q(t)-q_{ref}|$','Interpreter','latex');

title('Dominant Stress Coordinate Regulation Error');

legend( ...
    'LQG Controlled', ...
    'Observer-Based SMC Controlled', ...
    '$10^{-2}$ tolerance', ...
    'Location','best', ...
    'Interpreter','latex');

%% =========================================================
% Figure 15: SMC Observer Convergence
% ==========================================================

figure;

plot(t_smc,x_smc(:,1),'k','LineWidth',2);
hold on;

plot(t_smc,xhat_smc(:,1),'m--','LineWidth',2);

grid on;

xlabel('$t$ (s)','Interpreter','latex');
ylabel('$q(t),\ \hat{q}(t)$','Interpreter','latex');

title('SMC Observer Convergence');

legend( ...
    '$q(t)$', ...
    '$\hat{q}(t)$', ...
    'Location','best', ...
    'Interpreter','latex');

%% =========================================================
% Local Functions
% ==========================================================

function dx = nonlinear_dynamics(~,x,M,Cv,Kv,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e,phi_p_0,T_0,Ve_0)

    q   = x(1);
    qdot = x(2);
    Phi = x(3);
    T   = x(4);
    Ve  = x(5);

    dx = zeros(5,1);

    dx(1) = qdot;

    dx(2) = ( ...
        -Cv*qdot ...
        -Kv*q ...
        +alpha_B*Phi^2 ...
        +alpha_E*Ve^2 ...
        -alpha_B*phi_p_0^2 ...
        -alpha_E*Ve_0^2) / M;

    dx(3) = -gamma_Phi*(Phi - phi_p_0);

    dx(4) = ( ...
        k_E*Ve^2 ...
        +k_Phi*Phi^2 ...
        -h*(T - T_0) ...
        -k_E*Ve_0^2 ...
        -k_Phi*phi_p_0^2) / C_th;

    dx(5) = -(Ve - Ve_0)/tau_e;

end


function dz = lqg_dynamics(~,z,A,B,C,K,L)

    x = z(1:5);
    xhat = z(6:10);

    u = -K*xhat;

    dx = A*x + B*u;

    dxhat = ...
        A*xhat ...
        +B*u ...
        +L*(C*x - C*xhat);

    dz = [dx; dxhat];

end


function dz = smc_observer_dynamics(~,z,M,Cv,Kv,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,Ve_0,eta,k_s,phi_s,u_min,u_max, ...
    sat,A,B,C,L)

    %% -----------------------------------------------------
    % Plant states
    % ------------------------------------------------------

    x = z(1:5);

    q   = x(1);
    qdot = x(2);
    Phi = x(3);
    T   = x(4);
    Ve  = x(5);

    %% -----------------------------------------------------
    % Observer states
    % ------------------------------------------------------

    xhat = z(6:10);

    qhat    = xhat(1);
    qdothat = xhat(2);
    dPhihat = xhat(3);

    %% -----------------------------------------------------
    % Measured outputs in perturbation coordinates
    % ------------------------------------------------------

    y = [
        q;
        Phi - phi_p_0;
        T   - T_0;
        Ve  - Ve_0
    ];

    %% -----------------------------------------------------
    % Observer-Based Sliding Surface
    % ------------------------------------------------------

    s_hat = qdothat + eta*qhat;

    %% -----------------------------------------------------
    % Desired acceleration
    % ------------------------------------------------------

    qddot_des = ...
        -eta*qdothat ...
        -k_s*sat(s_hat/phi_s);

    %% -----------------------------------------------------
    % Equivalent electrostatic control
    % ------------------------------------------------------

    dVe_cmd = ...
        (M/(2*alpha_E*Ve_0)) * ( ...
            qddot_des ...
            +(Cv/M)*qdothat ...
            +(Kv/M)*qhat ...
            -(2*alpha_B*phi_p_0/M)*dPhihat);

    %% -----------------------------------------------------
    % Absolute voltage command
    % ------------------------------------------------------

    Ve_cmd = Ve_0 + dVe_cmd;

    %% Actuator saturation

    Ve_cmd = min(max(Ve_cmd,u_min),u_max);

    %% -----------------------------------------------------
    % Nonlinear plant dynamics
    % ------------------------------------------------------

    dx = zeros(5,1);

    dx(1) = qdot;

    dx(2) = ( ...
        -Cv*qdot ...
        -Kv*q ...
        +alpha_B*Phi^2 ...
        +alpha_E*Ve^2 ...
        -alpha_B*phi_p_0^2 ...
        -alpha_E*Ve_0^2) / M;

    dx(3) = -gamma_Phi*(Phi - phi_p_0);

    dx(4) = ( ...
        k_E*Ve^2 ...
        +k_Phi*Phi^2 ...
        -h*(T - T_0) ...
        -k_E*Ve_0^2 ...
        -k_Phi*phi_p_0^2) / C_th;

    dx(5) = (Ve_cmd - Ve)/tau_e;

    %% -----------------------------------------------------
    % Common LQE
    % ------------------------------------------------------

    delta_u = Ve_cmd - Ve_0;

    dxhat = ...
        A*xhat ...
        +B*delta_u ...
        +L*(y - C*xhat);

    dz = [dx; dxhat];

end