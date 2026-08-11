%% ==========================================================
% LQG Nominal Performance Study
%
% Nested Monte Carlo structure:
%   N_param = number of parameter realizations
%   N_IC    = number of initial-condition realizations
%
% Every parameter realization is simulated over ALL loaded
% initial conditions.
%
% Total closed-loop simulations = N_param * N_IC
% ==========================================================

clear;
clc;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ==========================================================
% Nominal Model Parameters
% ==========================================================

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

%% ==========================================================
% Monte Carlo Configuration
% ==========================================================

rng(7);

% Number of parameter realizations
%
% For the current nominal study, this is one realization.
% Increase this when parameter uncertainty distributions
% are introduced.
N_param = 10;

% Number of initial-condition realizations
N_IC = 1000;

Tsim = 10;
tspan = 0:0.01:Tsim;

% ----------------------------------------------------------
% Initial-condition perturbation bounds
% ----------------------------------------------------------

q_max       = 0.30;
q_dot_max   = 0.15;
dphi_max    = 0.10;
dT_max      = 0.08;
dV_e_max    = 0.15;

% ----------------------------------------------------------
% Convergence criterion
% ----------------------------------------------------------

eps_q     = 1e-2;
eps_q_dot = 1e-2;

% ----------------------------------------------------------
% Settling criterion
% ----------------------------------------------------------

settle_tol = 2e-2;

% ----------------------------------------------------------
% Physical actuator-voltage limits
% ----------------------------------------------------------

V_min = 0;
V_max = 10;

%% ==========================================================
% Generate Initial-Condition Ensemble ONCE
%
% Every parameter realization is evaluated over exactly the
% same set of initial conditions.
% ==========================================================

ICs = zeros(5,N_IC);

for j = 1:N_IC

    ICs(:,j) = [
        q_max      * (2*rand - 1);
        q_dot_max  * (2*rand - 1);
        dphi_max   * (2*rand - 1);
        dT_max     * (2*rand - 1);
        dV_e_max   * (2*rand - 1)
    ];

end

fprintf('\n');
fprintf('=========================================================\n');
fprintf('Monte Carlo Configuration\n');
fprintf('=========================================================\n');
fprintf('Parameter realizations: %d\n',N_param);
fprintf('Initial conditions per realization: %d\n',N_IC);
fprintf('Total closed-loop simulations: %d\n',N_param*N_IC);
fprintf('Simulation duration: %.2f s\n',Tsim);
fprintf('Physical voltage limits: [%.2f, %.2f] V\n',V_min,V_max);

%% ==========================================================
% Parameter Realization Storage
%
% For the current nominal study, the parameters remain at
% their nominal values. This structure is ready for adding
% parameter uncertainty later.
% ==========================================================

M_param          = zeros(N_param,1);
Cv_param         = zeros(N_param,1);
Kv_param         = zeros(N_param,1);
alpha_B_param    = zeros(N_param,1);
alpha_E_param    = zeros(N_param,1);
gamma_Phi_param  = zeros(N_param,1);
k_Phi_param      = zeros(N_param,1);
k_E_param        = zeros(N_param,1);
C_th_param       = zeros(N_param,1);
h_param          = zeros(N_param,1);
tau_e_param      = zeros(N_param,1);

for p = 1:N_param

    % ------------------------------------------------------
    % Current nominal realization
    %
    % Replace these assignments with random draws when
    % parameter uncertainty is introduced.
    % ------------------------------------------------------

    M_param(p)         = M;
    Cv_param(p)        = C_v;
    Kv_param(p)        = K_v;
    alpha_B_param(p)   = alpha_B;
    alpha_E_param(p)   = alpha_E;
    gamma_Phi_param(p) = gamma_Phi;
    k_Phi_param(p)     = k_Phi;
    k_E_param(p)       = k_E;
    C_th_param(p)      = C_th;
    h_param(p)         = h;
    tau_e_param(p)     = tau_e;

end

%% ==========================================================
% Preallocate Results
%
% One row per closed-loop simulation:
%
%   parameter realization x initial condition
% ==========================================================

N_total = N_param * N_IC;

LQG_RMS = zeros(N_total,1);
LQG_energy = zeros(N_total,1);
LQG_settle = zeros(N_total,1);
LQG_maxVoltage = zeros(N_total,1);

LQG_saturated = false(N_total,1);
LQG_saturationDuration = zeros(N_total,1);

LQG_converged = false(N_total,1);

% Store parameter and IC indices
Result_paramIndex = zeros(N_total,1);
Result_ICIndex = zeros(N_total,1);

%% ==========================================================
% Representative Trajectory Storage
% ==========================================================

representativeSaved = false;

%% ==========================================================
% Nested Monte Carlo Simulation
% ==========================================================

resultIndex = 0;

for p = 1:N_param

    %% ------------------------------------------------------
    % Parameter realization
    % ------------------------------------------------------

    M_i         = M_param(p);
    C_v_i       = Cv_param(p);
    K_v_i       = Kv_param(p);
    alpha_B_i   = alpha_B_param(p);
    alpha_E_i   = alpha_E_param(p);
    gamma_Phi_i = gamma_Phi_param(p);
    k_Phi_i     = k_Phi_param(p);
    k_E_i       = k_E_param(p);
    C_th_i      = C_th_param(p);
    h_i         = h_param(p);
    tau_e_i     = tau_e_param(p);

    %% ------------------------------------------------------
    % Linearized Reduced-Order Model
    % x = [q; qdot; Phi_p; T; V_e]
    % ------------------------------------------------------

    A = [
        0,        1,        0,                         0,             0;
        -K_v_i/M_i, -C_v_i/M_i, 2*alpha_B_i*phi_p_0/M_i, 0, ...
            2*alpha_E_i*V_e_0/M_i;
        0,        0,       -gamma_Phi_i,              0,             0;
        0,        0,        2*k_Phi_i*phi_p_0/C_th_i, ...
            -h_i/C_th_i, 2*k_E_i*V_e_0/C_th_i;
        0,        0,        0,                        0, ...
            -1/tau_e_i
    ];

    B = [
        0;
        0;
        0;
        0;
        1/tau_e_i
    ];

    C = [
        1 0 0 0 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 0 1
    ];

    D = zeros(4,1);

    %% ------------------------------------------------------
    % LQR State Feedback
    % ------------------------------------------------------

    Q_lqr = diag([100, 10, 5, 5, 1]);
    R_lqr = 1;

    K_lqr = lqr(A,B,Q_lqr,R_lqr);

    %% ------------------------------------------------------
    % Kalman State Estimator
    % ------------------------------------------------------

    W = diag([1e-3, 1e-2, 1e-4, 1e-4, 1e-3]);
    V = diag([1e-3, 1e-4, 1e-3, 1e-3]);

    L = lqe(A,eye(5),C,W,V);

    %% ------------------------------------------------------
    % Closed-loop / estimator diagnostics
    % ------------------------------------------------------

    A_cl = A - B*K_lqr;
    A_est = A - L*C;

    if p == 1

        disp(' ');
        disp('=========================================================');
        disp('LQG Controller Diagnostics');
        disp('=========================================================');

        disp('K_lqr = ');
        disp(K_lqr);

        disp('Closed-Loop Eigenvalues = ');
        disp(eig(A_cl));

        disp('Kalman Gain L = ');
        disp(L);

        disp('Estimator Eigenvalues = ');
        disp(eig(A_est));

    end

    %% ======================================================
    % Initial-Condition Ensemble
    % ======================================================

    for j = 1:N_IC

        resultIndex = resultIndex + 1;

        Result_paramIndex(resultIndex) = p;
        Result_ICIndex(resultIndex) = j;

        %% --------------------------------------------------
        % Retrieve shared initial condition
        % --------------------------------------------------

        dx_0 = ICs(:,j);

        %% --------------------------------------------------
        % Physical nonlinear plant initial condition
        % --------------------------------------------------

        x_0 = [
            dx_0(1);
            dx_0(2);
            phi_p_0 + dx_0(3);
            T_0     + dx_0(4);
            V_e_0   + dx_0(5)
        ];

        %% --------------------------------------------------
        % Observer initialized at nominal equilibrium
        %
        % Observer states are perturbation coordinates.
        % --------------------------------------------------

        xhat_0 = zeros(5,1);

        z_0 = [x_0; xhat_0];

        %% --------------------------------------------------
        % Nonlinear Plant + Kalman Observer
        % --------------------------------------------------

        f_lqg = @(t,z) lqg_nonlinear_dynamics( ...
            t,z,A,B,C,L,K_lqr,...
            M_i,C_v_i,K_v_i,...
            phi_p_0,T_0,V_e_0,...
            alpha_B_i,alpha_E_i,...
            gamma_Phi_i,k_Phi_i,k_E_i,...
            C_th_i,h_i,tau_e_i,...
            V_min,V_max);

        [t,z] = ode45(f_lqg,tspan,z_0);

        %% --------------------------------------------------
        % Extract States
        % --------------------------------------------------

        x = z(:,1:5);
        x_hat = z(:,6:10);

        q = x(:,1);
        q_dot = x(:,2);

        %% --------------------------------------------------
        % Reconstruct Actual Control Input
        %
        % x_hat is expressed in perturbation coordinates:
        %
        %     u_cmd = -K*x_hat
        %
        % This is a perturbation voltage command.
        %
        % Convert to physical voltage, saturate in [0,10],
        % then convert the saturated voltage back into
        % perturbation coordinates.
        % --------------------------------------------------

        u = zeros(length(t),1);
        V_e_cmd = zeros(length(t),1);

        for k = 1:length(t)

            % Perturbation-coordinate control command
            u_cmd = -K_lqr*x_hat(k,:)';

            % Convert to physical actuator voltage
            V_cmd = V_e_0 + u_cmd;

            % Apply physical actuator constraint
            V_cmd = max(V_min,min(V_max,V_cmd));

            % Convert actual applied voltage back to
            % perturbation coordinates
            u(k) = V_cmd - V_e_0;

            % Store physical voltage
            V_e_cmd(k) = V_cmd;

        end

        %% --------------------------------------------------
        % Performance Metrics
        % --------------------------------------------------

        % RMS regulation error
        LQG_RMS(resultIndex) = sqrt(mean(q.^2));

        % Control energy in perturbation voltage
        LQG_energy(resultIndex) = trapz(t,u.^2);

        % Maximum physical actuator voltage
        LQG_maxVoltage(resultIndex) = max(V_e_cmd);

        %% --------------------------------------------------
        % Physical actuator saturation
        % --------------------------------------------------

        saturated = ...
            (V_e_cmd <= V_min + 1e-10) | ...
            (V_e_cmd >= V_max - 1e-10);

        LQG_saturated(resultIndex) = any(saturated);

        LQG_saturationDuration(resultIndex) = ...
            100 * trapz(t,double(saturated)) / Tsim;

        %% --------------------------------------------------
        % Convergence
        % --------------------------------------------------

        LQG_converged(resultIndex) = ...
            abs(q(end)) < eps_q && ...
            abs(q_dot(end)) < eps_q_dot;

        %% --------------------------------------------------
        % Settling Time
        % --------------------------------------------------

        err_norm = sqrt(q.^2 + q_dot.^2);

        idx_settle = find( ...
            err_norm < settle_tol, ...
            1,'first');

        if isempty(idx_settle)

            LQG_settle(resultIndex) = Tsim;

        else

            remaining_err = err_norm(idx_settle:end);

            if all(remaining_err < settle_tol)

                LQG_settle(resultIndex) = t(idx_settle);

            else

                LQG_settle(resultIndex) = Tsim;

            end

        end

        %% --------------------------------------------------
        % Save One Representative Trajectory
        % --------------------------------------------------

        if ~representativeSaved

            representative_t = t;
            representative_q = q;
            representative_qhat = x_hat(:,1);
            representative_u = u;
            representative_Ve = V_e_cmd;

            representativeSaved = true;

        end

    end

    %% ------------------------------------------------------
    % Progress indicator
    % ------------------------------------------------------

    fprintf( ...
        'Parameter realization %d/%d complete (%d simulations)\n', ...
        p,N_param,p*N_IC);

end

%% ==========================================================
% Summary Statistics
% ==========================================================

SuccessRate = ...
    100 * mean(LQG_converged);

MeanEnergy = ...
    mean(LQG_energy);

MeanRMS = ...
    mean(LQG_RMS);

MeanSettlingTime = ...
    mean(LQG_settle);

MeanMaxVoltage = ...
    mean(LQG_maxVoltage);

SaturationIncidence = ...
    100 * mean(LQG_saturated);

MeanSaturationDuration = ...
    mean(LQG_saturationDuration);

%% ==========================================================
% Results Display
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('LQG Nominal Performance Results\n');
fprintf('=========================================================\n');

fprintf('Parameter realizations: %d\n',N_param);
fprintf('Initial conditions per realization: %d\n',N_IC);
fprintf('Total closed-loop simulations: %d\n',N_total);
fprintf('Simulation duration: %.2f s\n',Tsim);

fprintf('\nConvergence rate: %.2f %%\n',SuccessRate);
fprintf('Mean RMS regulation error: %.4f\n',MeanRMS);
fprintf('Mean settling time: %.4f s\n',MeanSettlingTime);
fprintf('Mean control energy: %.4f\n',MeanEnergy);
fprintf('Mean maximum physical voltage: %.4f V\n',MeanMaxVoltage);
fprintf('Saturation incidence: %.2f %%\n',SaturationIncidence);
fprintf('Mean saturation duration: %.4f %%\n',MeanSaturationDuration);

%% ==========================================================
% Results Table
% ==========================================================

LQG_results_table = table( ...
    SuccessRate,...
    MeanRMS,...
    MeanSettlingTime,...
    MeanEnergy,...
    MeanMaxVoltage,...
    SaturationIncidence,...
    MeanSaturationDuration);

disp(LQG_results_table);

%% ==========================================================
% Representative Phase Plot
% ==========================================================

figure;

plot( ...
    representative_q,...
    gradient(representative_q,representative_t),...
    'LineWidth',1.2);

grid off;

xlabel('$q$');
ylabel('$\dot{q}$');

title('LQG Nominal Phase Response');

set(gca,'FontSize',9);

%% ==========================================================
% Local Nonlinear LQG Dynamics Function
% ==========================================================

function dz = lqg_nonlinear_dynamics( ...
    t,z,A,B,C,L,K,...
    M,C_v,K_v,...
    phi_p_0,T_0,V_e_0,...
    alpha_B,alpha_E,...
    gamma_Phi,k_Phi,k_E,...
    C_th,h,tau_e,...
    V_min,V_max)

%% ----------------------------------------------------------
% Plant state and observer state
% ----------------------------------------------------------

x = z(1:5);
x_hat = z(6:10);

%% ----------------------------------------------------------
% Estimated-state feedback
%
% x_hat is in perturbation coordinates.
% ----------------------------------------------------------

u_cmd = -K*x_hat;

%% ----------------------------------------------------------
% Convert perturbation command to physical voltage
% ----------------------------------------------------------

V_e_cmd = V_e_0 + u_cmd;

%% ----------------------------------------------------------
% Physical actuator saturation
% ----------------------------------------------------------

V_e_cmd = max(V_min,min(V_max,V_e_cmd));

%% ----------------------------------------------------------
% Convert saturated physical voltage back to perturbation
% coordinates
% ----------------------------------------------------------

u = V_e_cmd - V_e_0;

%% ----------------------------------------------------------
% Nonlinear plant
% ----------------------------------------------------------

q     = x(1);
q_dot = x(2);
Phi_p = x(3);
T     = x(4);
V_e   = x(5);

x_dot = [
    q_dot;

    (-C_v*q_dot ...
    -K_v*q ...
    +alpha_B*Phi_p^2 ...
    +alpha_E*V_e^2 ...
    -alpha_B*phi_p_0^2 ...
    -alpha_E*V_e_0^2)/M;

    -gamma_Phi*(Phi_p - phi_p_0);

    (k_E*V_e^2 ...
    +k_Phi*Phi_p^2 ...
    -h*(T-T_0) ...
    -k_E*V_e_0^2 ...
    -k_Phi*phi_p_0^2)/C_th;

    -(V_e-V_e_0)/tau_e ...
    + u/tau_e
];

%% ----------------------------------------------------------
% Convert physical plant state to perturbation coordinates
% ----------------------------------------------------------

x_pert = [
    x(1);
    x(2);
    x(3) - phi_p_0;
    x(4) - T_0;
    x(5) - V_e_0
];

%% ----------------------------------------------------------
% Measurements in perturbation coordinates
% ----------------------------------------------------------

y = C*x_pert;

y_hat = C*x_hat;

%% ----------------------------------------------------------
% Linear Kalman observer
% ----------------------------------------------------------

x_hat_dot = ...
    A*x_hat ...
    + B*u ...
    + L*(y-y_hat);

%% ----------------------------------------------------------
% Combined dynamics
% ----------------------------------------------------------

dz = [
    x_dot;
    x_hat_dot
];

end