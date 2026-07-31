%% LQG Formulation

set(groot,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

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

%% Linearized Model
A = [0,        1,        0,                         0,             0;
    -K_v/M,   -C_v/M,   2*alpha_B*phi_p_0/M,        0,             2*alpha_E*V_e_0/M;
     0,        0,       -gamma_Phi,                  0,             0;
     0,        0,        2*k_Phi*phi_p_0/C_th,      -h/C_th,        2*k_E*V_e_0/C_th;
     0,        0,        0,                         0,            -1/tau_e];

B = [0; 0; 0; 0; 1/tau_e];

C = [1 0 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];

D = zeros(4,1);

%% LQR Gain
Q_lqr = diag([100, 10, 5, 5, 1]);
R_lqr = 1;

K_lqr = lqr(A,B,Q_lqr,R_lqr);

A_cl = A - B*K_lqr;
eigs_A_cl = eig(A_cl);

disp('K_lqr = ');
disp(K_lqr);

disp('Closed-Loop Eigenvalues = ');
disp(eigs_A_cl);

%% Kalman Filter Gain
W = diag([1e-3, 1e-2, 1e-4, 1e-4, 1e-3]);   % Process noise covariance
V = diag([1e-3, 1e-4, 1e-3, 1e-3]);         % Measurement noise covariance

L = lqe(A, eye(5), C, W, V);

A_est = A - L*C;
eigs_A_est = eig(A_est);

disp('Kalman Gain L = ');
disp(L);

disp('Estimator Eigenvalues = ');
disp(eigs_A_est);

%% Multiple Initial Condition Study

rng(7);                     % Reproducibility
num_IC = 1000;                 % Number of initial conditions
t_final = 10;               % Simulation time, seconds
tspan = linspace(0,t_final,1000);

% Bounded perturbation set around nominal equilibrium
q_max      = 0.30;
q_dot_max   = 0.15;
dphi_max   = 0.10;
dT_max     = 0.08;
dV_e_max    = 0.15;

% Convergence tolerances
eps_q    = 1e-2;
eps_q_dot = 1e-2;
eps_norm = 2e-2;

terminalData = zeros(num_IC,6);
converged_lqg = false(num_IC,1);
converged_nl  = false(num_IC,1);

figure;
hold on; grid on;

for i = 1:num_IC

    % Random perturbation initial condition
    dx_0 = [
        q_max    *(2*rand - 1);
        q_dot_max *(2*rand - 1);
        dphi_max *(2*rand - 1);
        dT_max   *(2*rand - 1);
        dV_e_max  *(2*rand - 1)
    ];

    % Linear perturbation-state IC
    x_0_lin = dx_0;
    x_hat_0  = zeros(5,1);
    z_0     = [x_0_lin; x_hat_0];

    % Nonlinear physical-coordinate IC
    x_0_nl = [
        dx_0(1);
        dx_0(2);
        phi_p_0 + dx_0(3);
        T_0     + dx_0(4);
        V_e_0   + dx_0(5)
    ];

    %% Nonlinear Plant Model
    % x = [q; qdot; Phi_p; T; V_e]
    f_nl = @(t,x) [x(2); ...
            (-C_v*x(2) -K_v*x(1) +alpha_B*x(3)^2 +alpha_E*x(5)^2 -alpha_B*phi_p_0^2 -alpha_E*V_e_0^2)/M; ...
            -gamma_Phi*(x(3) -phi_p_0); ...
            (k_E*x(5)^2 +k_Phi*x(3)^2 -h*(x(4) -T_0) -k_E*V_e_0^2 -k_Phi*phi_p_0^2)/C_th; ...
            -(x(5) - V_e_0)/tau_e];
    
    [~,x_nl] = ode45(f_nl,tspan,x_0_nl);
    q_nl_final    = x_nl(end,1);
    qdot_nl_final = x_nl(end,2);

    converged_nl(i) = ...
        abs(q_nl_final) < eps_q && ...
        abs(qdot_nl_final) < eps_q_dot;


%% Combined Controlled Plant + Estimator
% z = [x; x_hat]
f_lqg = @(t,z) [A*z(1:5) + B*(-K_lqr*z(6:10)); A*z(6:10) + B*(-K_lqr*z(6:10)) + L*(C*z(1:5) - C*z(6:10))];

%% LQG controlled simulation
[~,z_lqg] = ode45(f_lqg,tspan,z_0);

x_cl    = z_lqg(:,1:5);
x_hat_cl = z_lqg(:,6:10);

q_nl    = x_nl(:,1);
q_dot_nl = x_nl(:,2);

q_lqg_final    = x_cl(end,1);
qdot_lqg_final = x_cl(end,2);

terminal_norm = norm([q_lqg_final; qdot_lqg_final],2);

converged_lqg(i) = ...
    abs(q_lqg_final) < eps_q && ...
    abs(qdot_lqg_final) < eps_q_dot && ...
    terminal_norm < eps_norm;

terminalData(i,:) = [
    dx_0(1), dx_0(2), q_lqg_final, qdot_lqg_final, ...
    q_nl_final, qdot_nl_final
];

    %% Plot representative trajectories only
    if i <= 10
        plot(q_nl,q_dot_nl,'r'); hold on;
        plot(x_cl(:,1),x_cl(:,2),'b');
        plot(x_hat_cl(:,1),x_hat_cl(:,2),'--','Color', '[1 0.65 0.12]');
    end
end

xlabel('$q$');
ylabel('$\dot{q}$');
title('Phase Response Under Multiple Initial Conditions');

legend('Nonlinear Uncontrolled', ...
       'LQG Controlled True', ...
       'LQG Estimated', ...
       'Location','best');

%% Report Results

fprintf('\nInitial-condition study results:\n');
fprintf('Number of ICs tested: %d\n',num_IC);
fprintf('Simulation time: %.2f s\n',t_final);
fprintf('Perturbation bounds:\n');
fprintf('|q(0)| <= %.3f\n',q_max);
fprintf('|qdot(0)| <= %.3f\n',q_dot_max);
fprintf('|dPhi(0)| <= %.3f\n',dphi_max);
fprintf('|dT(0)| <= %.3f\n',dT_max);
fprintf('|dVe(0)| <= %.3f\n',dV_e_max);

fprintf('\nConvergence criterion:\n');
fprintf('|q(tf)| < %.3e\n',eps_q);
fprintf('|qdot(tf)| < %.3e\n',eps_q_dot);
fprintf('norm([q(tf); qdot(tf)]) < %.3e\n',eps_norm);

fprintf('\nNonlinear uncontrolled convergence rate: %.1f%%\n', ...
    100*mean(converged_nl));

fprintf('LQG controlled convergence rate: %.1f%%\n', ...
    100*mean(converged_lqg));

x_0_lin = dx_0;
xhat_0 = zeros(5,1);
z_0 = [x_0_lin; xhat_0];
x_0_nl = [
    dx_0(1);
    dx_0(2);
    phi_p_0 + dx_0(3);
    T_0     + dx_0(4);
    V_e_0   + dx_0(5)
];

%% Phase Plot: q vs q_dot

figure;
plot(q_nl,q_dot_nl,'r'); hold on;
plot(x_cl(:,1),x_cl(:,2),'b');
plot(x_hat_cl(:,1),x_hat_cl(:,2),'--','Color', '[1 0.65 0.12]');

grid on;
xlabel('$q$');
ylabel('$\dot{q}$');
legend('Uncontrolled Nonlinear','LQG True','LQG Estimated', ...
       'Location','best');

title('LQG Phase Response');

Nmc = 5000;
tspan = 0:0.01:10;

Q_list = {
    diag([10 1 1 1 1])
    diag([50 5 1 1 1])
    diag([100 10 1 1 1])
    diag([250 25 1 1 1])
    diag([500 50 1 1 1])
    diag([1000 100 1 1 1])
};

R = 1;

q_tol = 1e-2;
qdot_tol = 1e-2;

results_lqr = [];

for iq = 1:length(Q_list)

    Q = Q_list{iq};
    K_test = lqr(A,B,Q,R);

    conv_count = 0;
    settle_times = nan(Nmc,1);
    control_energy = nan(Nmc,1);
    max_u = nan(Nmc,1);

    for imc = 1:Nmc

        dx0 = [
            0.30*(2*rand-1);
            0.15*(2*rand-1);
            0.10*(2*rand-1);
            0.08*(2*rand-1);
            0.15*(2*rand-1)
        ];

        f_cl = @(t,x) A*x + B*(-K_test*x);

        [t,x] = ode45(f_cl,tspan,dx0);

        u = -(K_test*x')';

        q = x(:,1);
        qdot = x(:,2);

        final_ok = abs(q(end)) < q_tol && abs(qdot(end)) < qdot_tol;

        if final_ok
            conv_count = conv_count + 1;
        end

        err_norm = vecnorm([q qdot],2,2);
        idx = find(err_norm < 2e-2,1,'first');

        if ~isempty(idx)
            settle_times(imc) = t(idx);
        end

        control_energy(imc) = trapz(t,u.^2);
        max_u(imc) = max(abs(u));

    end

    results_lqr(iq).Q = Q;
    results_lqr(iq).K = K_test;
    results_lqr(iq).convergence_rate = conv_count/Nmc;
    results_lqr(iq).mean_settling_time = mean(settle_times,'omitnan');
    results_lqr(iq).mean_control_energy = mean(control_energy,'omitnan');
    results_lqr(iq).mean_max_u = mean(max_u,'omitnan');

end

fprintf('\nLQR Q Sweep Results:\n');
for iq = 1:length(results_lqr)
    fprintf('Q case %d: conv = %.1f%%, Ts = %.3f s, Energy = %.3f, max|u| = %.3f\n', ...
        iq, ...
        100*results_lqr(iq).convergence_rate, ...
        results_lqr(iq).mean_settling_time, ...
        results_lqr(iq).mean_control_energy, ...
        results_lqr(iq).mean_max_u);
end

