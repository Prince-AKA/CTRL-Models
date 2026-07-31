%% 8.6 Observer-Based SMC Comparison: Time Histories and Error Response
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

tspan = 0:0.01:10;

%% Linearized Model

A = [
    0,        1,        0,                         0,             0;
   -K_v/M,   -C_v/M,   2*alpha_B*phi_p_0/M,        0,             2*alpha_E*V_e_0/M;
    0,        0,       -gamma_Phi,                 0,             0;
    0,        0,        2*k_Phi*phi_p_0/C_th,     -h/C_th,        2*k_E*V_e_0/C_th;
    0,        0,        0,                         0,            -1/tau_e
];

B = [0; 0; 0; 0; 1/tau_e];

C = [
    1 0 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1
];

D = zeros(4,1);

%% LQR and Observer Gains

Q_lqr = diag([100 10 5 5 1]);
R_lqr = 1;

K_lqr = lqr(A,B,Q_lqr,R_lqr);

observer_poles_lqg = [-10 -9 -8 -7 -6];
L = place(A',C',observer_poles_lqg)';

observer_poles_smc = [-8 -9 -10 -11 -12];
L_o = place(A',C',observer_poles_smc)';

%% Representative Initial Condition

dx_0 = [0.20; 0.00; 0.05; 0.03; 0.10];

x_0_nl = [
    dx_0(1);
    dx_0(2);
    phi_p_0 + dx_0(3);
    T_0     + dx_0(4);
    V_e_0   + dx_0(5)
];

x_0_lin = dx_0;
xhat_0  = zeros(5,1);

z_0_lqg = [x_0_lin; xhat_0];

%% Uncontrolled Nonlinear Response

f_nl = @(t,x) nonlinear_dynamics(t,x,M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e,phi_p_0,T_0,V_e_0);

[t_nl,x_nl] = ode45(f_nl,tspan,x_0_nl);

%% LQG Controlled Response

f_lqg = @(t,z) lqg_dynamics(t,z,A,B,C,K_lqr,L);

[t_lqg,z_lqg] = ode45(f_lqg,tspan,z_0_lqg);

x_lqg     = z_lqg(:,1:5);
xhat_lqg  = z_lqg(:,6:10);

%% Observer-Based SMC Parameters

eta_s = 2.0;
Delta_max = 248.7;
k_s = 1.05*Delta_max;
eps_s = 0.05;

u_min = 0.0;
u_max = 3.0;

sat = @(z) max(-1,min(1,z));

%% Observer-Based SMC Response

xhat_0_smc = [dx_0(1); 0.0; dx_0(3); dx_0(4); dx_0(5)];    
z_0_smc = [x_0_nl; xhat_0_smc];

f_smc = @(t,z) smc_observer_dynamics(t,z,M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,V_e_0,eta_s,k_s,eps_s,u_min,u_max,sat,A,B,C,L_o);

[t_smc,z_smc] = ode45(f_smc,tspan,z_0_smc);

x_smc    = z_smc(:,1:5);
xhat_smc = z_smc(:,6:10);

%% Figure 13: Stress-Coordinate Time History

figure;
plot(t_nl,x_nl(:,1),'r','LineWidth',2); hold on;
plot(t_lqg,x_lqg(:,1),'b','LineWidth',2);
plot(t_smc,x_smc(:,1),'k--','LineWidth',2);

yline(1e-2,'--','LineWidth',1);
yline(-1e-2,'--','LineWidth',1);

grid on;
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$q(t)$','Interpreter','latex');
title('Dominant Stress Coordinate Regulation');

legend('Nonlinear Uncontrolled','LQG Controlled','Observer-Based SMC Controlled', ...
       '$\pm 10^{-2}$ tolerance','Location','best','Interpreter','latex');

%% Figure 14: Absolute Regulation Error

q_ref = 0;

e_lqg = abs(x_lqg(:,1) - q_ref);
e_smc = abs(x_smc(:,1) - q_ref);

figure;
semilogy(t_lqg,e_lqg,'b','LineWidth',2); hold on;
semilogy(t_smc,e_smc,'k--','LineWidth',2);
yline(1e-2,'r--','LineWidth',0.5);

grid on;
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$|q(t)-q_{ref}|$','Interpreter','latex');
title('Dominant Stress Coordinate Regulation Error');

legend('LQG Controlled','Observer-Based SMC Controlled','$10^{-2}$ tolerance', ...
       'Location','best','Interpreter','latex');

%% Optional Figure: SMC Observer Convergence

figure;
plot(t_smc,x_smc(:,1),'k','LineWidth',2); hold on;
plot(t_smc,xhat_smc(:,1),'m--','LineWidth',2);

grid on;
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$q(t),\ \hat{q}(t)$','Interpreter','latex');
title('SMC Observer Convergence');

legend('$q(t)$','$\hat{q}(t)$','Location','best','Interpreter','latex');

%% Optional SMC eta/eps Sweep

run_smc_sweep = 1;

if run_smc_sweep

    eta_list = [0.5 1.0 1.5 2.0 3.0 4.0 5.0];
    eps_list = [0.01 0.025 0.05 0.10 0.20];

    Nmc = 1000;

    q_tol = 1e-2;
    qdot_tol = 1e-2;

    results_smc = [];
    idx_result = 0;

    for ieta = 1:length(eta_list)

        eta = eta_list(ieta);

        for ieps = 1:length(eps_list)

            eps_s = eps_list(ieps);

            conv_count = 0;
            settle_times = nan(Nmc,1);
            control_energy = nan(Nmc,1);
            max_u = nan(Nmc,1);
            max_s = nan(Nmc,1);
            mean_obs_error = nan(Nmc,1);

            parfor imc = 1:Nmc

                dx0 = [
                    0.30*(2*rand-1);
                    0.15*(2*rand-1);
                    0.10*(2*rand-1);
                    0.08*(2*rand-1);
                    0.15*(2*rand-1)
                ];

                x0_nl = [
                    dx0(1);
                    dx0(2);
                    phi_p_0 + dx0(3);
                    T_0     + dx0(4);
                    V_e_0   + dx0(5)
                ];

                z0_smc = [x0_nl; zeros(5,1)];

                f_smc = @(t,z) smc_observer_dynamics(t,z,M,C_v,K_v,alpha_B,alpha_E, ...
                    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
                    phi_p_0,T_0,V_e_0,eta,k_s,eps_s,u_min,u_max,sat,A,B,C,L_o);

                [t,z] = ode23(f_smc,tspan,z0_smc);

                x = z(:,1:5);
                xhat = z(:,6:10);

                q = x(:,1);
                qdot = x(:,2);
                Ve = x(:,5);

                qhat = xhat(:,1);
                qdothat = xhat(:,2);

                s_hat = qdothat + eta*qhat;

                final_ok = abs(q(end)) < q_tol && abs(qdot(end)) < qdot_tol;

                if final_ok
                    conv_count = conv_count + 1;
                end

                err_norm = vecnorm([q qdot],2,2);
                idt = find(err_norm < 2e-2,1,'first');

                if ~isempty(idt)
                    settle_times(imc) = t(idt);
                end

                control_energy(imc) = trapz(t,(Ve - V_e_0).^2);
                max_u(imc) = max(abs(Ve));
                max_s(imc) = max(abs(s_hat));

                true_delta_x = [
                    x(:,1), ...
                    x(:,2), ...
                    x(:,3)-phi_p_0, ...
                    x(:,4)-T_0, ...
                    x(:,5)-V_e_0
                ];

                obs_err = true_delta_x - xhat;
                mean_obs_error(imc) = mean(vecnorm(obs_err,2,2));

            end

            idx_result = idx_result + 1;

            results_smc(idx_result).eta = eta;
            results_smc(idx_result).eps_s = eps_s;
            results_smc(idx_result).k_s = k_s;
            results_smc(idx_result).convergence_rate = conv_count/Nmc;
            results_smc(idx_result).mean_settling_time = mean(settle_times,'omitnan');
            results_smc(idx_result).mean_control_energy = mean(control_energy,'omitnan');
            results_smc(idx_result).mean_max_u = mean(max_u,'omitnan');
            results_smc(idx_result).mean_max_s = mean(max_s,'omitnan');
            results_smc(idx_result).mean_obs_error = mean(mean_obs_error,'omitnan');

        end
    end

    fprintf('\nObserver-Based SMC eta/eps Sweep Results:\n');

    for i = 1:length(results_smc)
        fprintf('eta = %.2f, eps = %.3f: conv = %.1f%%, Ts = %.3f s, Energy = %.3f, max|Ve| = %.3f, obs err = %.3e\n', ...
            results_smc(i).eta, ...
            results_smc(i).eps_s, ...
            100*results_smc(i).convergence_rate, ...
            results_smc(i).mean_settling_time, ...
            results_smc(i).mean_control_energy, ...
            results_smc(i).mean_max_u, ...
            results_smc(i).mean_obs_error);
    end

end

%% Local Functions

function dx = nonlinear_dynamics(~,x,M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e,phi_p_0,T_0,V_e_0)

    q    = x(1);
    qdot = x(2);
    Phi  = x(3);
    T    = x(4);
    Ve   = x(5);

    dx = zeros(5,1);

    dx(1) = qdot;

    dx(2) = ( -C_v*qdot ...
              - K_v*q ...
              + alpha_B*Phi^2 ...
              + alpha_E*Ve^2 ...
              - alpha_B*phi_p_0^2 ...
              - alpha_E*V_e_0^2 ) / M;

    dx(3) = -gamma_Phi*(Phi - phi_p_0);

    dx(4) = ( k_E*Ve^2 ...
              + k_Phi*Phi^2 ...
              - h*(T - T_0) ...
              - k_E*V_e_0^2 ...
              - k_Phi*phi_p_0^2 ) / C_th;

    dx(5) = -(Ve - V_e_0)/tau_e;

end

function dz = lqg_dynamics(~,z,A,B,C,K_lqr,L)

    x = z(1:5);
    xhat = z(6:10);

    u = -K_lqr*xhat;

    dx = A*x + B*u;
    dxhat = A*xhat + B*u + L*(C*x - C*xhat);

    dz = [dx; dxhat];

end

function dz = smc_observer_dynamics(~,z,M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,V_e_0,eta,k_s,eps_s,u_min,u_max,sat,A,B,C,L_o)

    %% Nonlinear plant states, absolute coordinates
    x = z(1:5);

    q    = x(1);
    qdot = x(2);
    Phi  = x(3);
    T    = x(4);
    Ve   = x(5);

    %% Observer states, perturbation coordinates
    xhat = z(6:10);

    qhat    = xhat(1);
    qdothat = xhat(2);
    dPhihat = xhat(3);

    %% Measured output in perturbation coordinates
    y = [q; Phi - phi_p_0; T - T_0; Ve  - V_e_0];

    %% Observer-based sliding surface
    s_hat = qdothat + eta*qhat;

    %% Desired acceleration from estimated sliding condition
    qddot_des = -eta*qdothat - k_s*sat(s_hat/eps_s);

    %% Electrostatic command computed from estimated state
    dVe_cmd = (M/(2*alpha_E*V_e_0))*( ...
        qddot_des ...
        + (C_v/M)*qdothat ...
        + (K_v/M)*qhat ...
        - (2*alpha_B*phi_p_0/M)*dPhihat );

    u = V_e_0 + dVe_cmd;

    %% Actuator saturation
    u = min(max(u,u_min),u_max);

    %% Nonlinear plant dynamics
    dx = zeros(5,1);

    dx(1) = qdot;

    dx(2) = ( -C_v*qdot ...
              - K_v*q ...
              + alpha_B*Phi^2 ...
              + alpha_E*Ve^2 ...
              - alpha_B*phi_p_0^2 ...
              - alpha_E*V_e_0^2 ) / M;

    dx(3) = -gamma_Phi*(Phi - phi_p_0);

    dx(4) = ( k_E*Ve^2 ...
              + k_Phi*Phi^2 ...
              - h*(T - T_0) ...
              - k_E*V_e_0^2 ...
              - k_Phi*phi_p_0^2 ) / C_th;

    dx(5) = (u - Ve)/tau_e;

    %% Luenberger observer dynamics
    delta_u = u - V_e_0;

    dxhat = A*xhat + B*delta_u + L_o*(y - C*xhat);

    dz = [dx; dxhat];

end