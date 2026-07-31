%% Script 3: PID

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; clc; close all;

%% Model Parameters

M = 1.0;
C_v = 0.3;
K_v = 2.25;

phi_p_0 = 1.0;
T_0 = 0.0;
V_e_0 = 1.0;

alpha_B = 0.20;
alpha_E = 0.12;

gamma_Phi = 0.15;
k_Phi = 0.04;
k_E = 0.03;
C_th = 1.0;
h = 0.25;
tau_e = 0.10;

%% PID Gains

Kp = 8.0;
Ki = 2.0;
Kd = 1.0;

u_min = 0.0;
u_max = 3.0;

%% Initial Condition Study

rng(7);

num_IC = 1000;
t_final = 10;
tspan = linspace(0,t_final,1000);

q_max     = 0.30;
q_dot_max = 0.15;
dphi_max  = 0.10;
dT_max    = 0.08;
dV_e_max  = 0.15;

eps_q     = 1e-2;
eps_q_dot = 1e-2;
eps_norm  = 2e-2;

converged_pid = false(num_IC,1);
converged_nl  = false(num_IC,1);

settle_times_pid = nan(num_IC,1);
control_energy_pid = nan(num_IC,1);
max_Ve_pid = nan(num_IC,1);

figure;
hold on; grid on;

for i = 1:num_IC

    dx_0 = [
        q_max     *(2*rand - 1);
        q_dot_max *(2*rand - 1);
        dphi_max  *(2*rand - 1);
        dT_max    *(2*rand - 1);
        dV_e_max  *(2*rand - 1)
    ];

    x_0_nl = [
        dx_0(1);
        dx_0(2);
        phi_p_0 + dx_0(3);
        T_0     + dx_0(4);
        V_e_0   + dx_0(5)
    ];

    % PID augmented state:
    % z = [q; qdot; Phi_p; T; Ve; I_e]
    % I_e = integral of tracking error q - q_ref
    z_0_pid = [x_0_nl; 0];

    %% Uncontrolled nonlinear plant

    f_nl = @(t,x) [
        x(2);
        (-C_v*x(2) ...
        - K_v*x(1) ...
        + alpha_B*x(3)^2 ...
        + alpha_E*x(5)^2 ...
        - alpha_B*phi_p_0^2 ...
        - alpha_E*V_e_0^2)/M;
        -gamma_Phi*(x(3) - phi_p_0);
        (k_E*x(5)^2 ...
        + k_Phi*x(3)^2 ...
        - h*(x(4) - T_0) ...
        - k_E*V_e_0^2 ...
        - k_Phi*phi_p_0^2)/C_th;
        -(x(5) - V_e_0)/tau_e
    ];

    [t_nl,x_nl] = ode45(f_nl,tspan,x_0_nl);

    q_nl = x_nl(:,1);
    qdot_nl = x_nl(:,2);

    q_nl_final = q_nl(end);
    qdot_nl_final = qdot_nl(end);

    converged_nl(i) = ...
        abs(q_nl_final) < eps_q && ...
        abs(qdot_nl_final) < eps_q_dot;

    %% PID controlled nonlinear plant

    f_pid = @(t,z) pid_dynamics(t,z,M,C_v,K_v,alpha_B,alpha_E, ...
        gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
        phi_p_0,T_0,V_e_0,Kp,Ki,Kd,u_min,u_max);

    [t_pid,z_pid] = ode45(f_pid,tspan,z_0_pid);

    q_pid = z_pid(:,1);
    qdot_pid = z_pid(:,2);
    Ve_pid = z_pid(:,5);

    q_pid_final = q_pid(end);
    qdot_pid_final = qdot_pid(end);

    terminal_norm = norm([q_pid_final; qdot_pid_final],2);

    converged_pid(i) = ...
        abs(q_pid_final) < eps_q && ...
        abs(qdot_pid_final) < eps_q_dot && ...
        terminal_norm < eps_norm;

    err_norm = vecnorm([q_pid qdot_pid],2,2);
    idx = find(err_norm < eps_norm,1,'first');

    if ~isempty(idx)
        settle_times_pid(i) = t_pid(idx);
    end

    control_energy_pid(i) = trapz(t_pid,(Ve_pid - V_e_0).^2);
    max_Ve_pid(i) = max(abs(Ve_pid));

    %% Plot representative trajectories only

    if i <= 10
        plot(q_nl,qdot_nl,'r'); hold on;
        plot(q_pid,qdot_pid,'g');
    end
end

xlabel('$q$');
ylabel('$\dot{q}$');
title('PID Phase Response Under Multiple Initial Conditions');

legend('Nonlinear Uncontrolled','PID Controlled','Location','best');

%% Report Results

fprintf('\nPID initial-condition study results:\n');
fprintf('Number of ICs tested: %d\n',num_IC);
fprintf('Simulation time: %.2f s\n',t_final);

fprintf('\nPID gains:\n');
fprintf('Kp = %.3f\n',Kp);
fprintf('Ki = %.3f\n',Ki);
fprintf('Kd = %.3f\n',Kd);

fprintf('\nConvergence criterion:\n');
fprintf('|q(tf)| < %.3e\n',eps_q);
fprintf('|qdot(tf)| < %.3e\n',eps_q_dot);
fprintf('norm([q(tf); qdot(tf)]) < %.3e\n',eps_norm);

fprintf('\nNonlinear uncontrolled convergence rate: %.1f%%\n', ...
    100*mean(converged_nl));

fprintf('PID controlled convergence rate: %.1f%%\n', ...
    100*mean(converged_pid));

fprintf('\nPID performance metrics:\n');
fprintf('Mean settling time: %.3f s\n',mean(settle_times_pid,'omitnan'));
fprintf('Mean control energy: %.3f\n',mean(control_energy_pid,'omitnan'));
fprintf('Mean max |Ve|: %.3f\n',mean(max_Ve_pid,'omitnan'));

%% Representative Single-Trajectory Time Response

dx_0 = [0.20; 0.00; 0.05; 0.03; 0.10];

x_0_nl = [
    dx_0(1);
    dx_0(2);
    phi_p_0 + dx_0(3);
    T_0     + dx_0(4);
    V_e_0   + dx_0(5)
];

z_0_pid = [x_0_nl; 0];

[t_nl,x_nl] = ode45(f_nl,tspan,x_0_nl);
[t_pid,z_pid] = ode45(f_pid,tspan,z_0_pid);

figure;
plot(t_nl,x_nl(:,1),'r','LineWidth',2); hold on;
plot(t_pid,z_pid(:,1),'g','LineWidth',2);
yline(1e-2,'k--','LineWidth',1);
yline(-1e-2,'k--','LineWidth',1);

grid on;
xlabel('$t$ (s)');
ylabel('$q(t)$');
title('PID Dominant Stress Coordinate Regulation');
legend('Nonlinear Uncontrolled','PID Controlled','$\pm 10^{-2}$ tolerance', ...
    'Location','best');

figure;
semilogy(t_pid,abs(z_pid(:,1)),'g','LineWidth',2); hold on;
yline(1e-2,'r--','LineWidth',0.5);

grid on;
xlabel('$t$ (s)');
ylabel('$|q(t)-q_{ref}|$');
title('PID Baseline Regulation Error');
legend('PID Controlled','$10^{-2}$ tolerance','Location','best');

%% Local PID Dynamics Function

function dz = pid_dynamics(~,z,M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,V_e_0,Kp,Ki,Kd,u_min,u_max)

    q    = z(1);
    qdot = z(2);
    Phi  = z(3);
    T    = z(4);
    Ve   = z(5);
    I_e  = z(6);

    q_ref = 0;
    e = q - q_ref;

    dPhi = Phi - phi_p_0;

    %% PID acceleration command
    qddot_des = -Kp*e - Ki*I_e - Kd*qdot;

    %% Linearized inverse electrostatic command
    dVe_cmd = (M/(2*alpha_E*V_e_0))*( ...
        qddot_des ...
        + (C_v/M)*qdot ...
        + (K_v/M)*q ...
        - (2*alpha_B*phi_p_0/M)*dPhi );

    u = V_e_0 + dVe_cmd;

    %% Actuator saturation
    u = min(max(u,u_min),u_max);

    dz = zeros(6,1);

    dz(1) = qdot;

    dz(2) = ( -C_v*qdot ...
              - K_v*q ...
              + alpha_B*Phi^2 ...
              + alpha_E*Ve^2 ...
              - alpha_B*phi_p_0^2 ...
              - alpha_E*V_e_0^2 ) / M;

    dz(3) = -gamma_Phi*(Phi - phi_p_0);

    dz(4) = ( k_E*Ve^2 ...
              + k_Phi*Phi^2 ...
              - h*(T - T_0) ...
              - k_E*V_e_0^2 ...
              - k_Phi*phi_p_0^2 ) / C_th;

    dz(5) = (u - Ve)/tau_e;

    %% Integral error dynamics
    dz(6) = e;

end