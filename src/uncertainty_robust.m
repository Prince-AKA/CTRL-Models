%% =========================================================
% LQG vs Observer-Based SMC Robustness Study
%
% Monte Carlo robustness comparison using ALL loaded initial
% conditions from MC_initial_conditions_1000.mat.
%
% Each Monte Carlo realization uses:
%   - one uncertain plant realization
%   - one loaded initial condition
%   - one disturbance realization
%
% The SAME uncertain plant, initial condition, and disturbance
% realization are used for both LQG and SMC, providing a paired
% controller comparison.
%
% Physical actuator constraint:
%       Ve_min <= Ve <= Ve_max
%
% Controller operates in perturbation coordinates:
%       dVe = Ve - Ve_0
%
% Actuator convention:
%
%   Ve_des = SMC desired/virtual physical voltage.
%
%   Ve_cmd = physical actuator command BEFORE hard saturation.
%            Saturation statistics are based exclusively on
%            whether Ve_cmd violates [Ve_min, Ve_max].
%
%   Ve     = actual physical actuator state:
%            Ve = Ve_0 + x(5)
%
% Robust convergence criterion:
%
%   RMS(q - q_ref) over final 0.5 s < 0.05
%   RMS(qdot)       over final 0.5 s < 0.05
%
% IMPORTANT:
%   qdot is taken directly from plant state x(2).
% =========================================================

clear;
clc;
close all;

%% =========================================================
% Initial Conditions
% ==========================================================

IC_file = 'MC_initial_conditions_1000.mat';

if ~isfile(IC_file)
    error('Initial-condition file not found: %s',IC_file);
end

load(IC_file,'ICs');

if ~isnumeric(ICs) || ndims(ICs) ~= 2
    error('ICs must be a numeric two-dimensional matrix.');
end

if size(ICs,2) ~= 5
    error('ICs must be an N_IC x 5 matrix.');
end

N_IC = size(ICs,1);

if N_IC < 1
    error('ICs must contain at least one initial condition.');
end

fprintf('\nLoaded %d initial conditions from %s.\n', ...
    N_IC,IC_file);

%% =========================================================
% Reproducibility
% =========================================================
%
% Keep a fixed seed during development so that changes to the
% code can be compared using exactly the same random ensemble.
%
% Change to rng('shuffle') for an independent ensemble after
% the implementation has been validated.
% ==========================================================

rng(1,'twister');

%% =========================================================
% Nominal Parameters
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
k_Phi     = 0.04;
k_E       = 0.03;

C_th  = 1.0;
h     = 0.25;
tau_e = 0.10;

%% =========================================================
% Simulation Parameters
% ==========================================================

dt   = 1e-3;
Tsim = 10;

t  = 0:dt:Tsim;
Nt = length(t);

% IMPORTANT:
% Use every loaded initial condition.
N = N_IC;

% Perturbation-coordinate reference
q_ref = 0;

% Physical actuator limits
Ve_min = 0.0;
Ve_max = 10.0;

% Actuator command tracking rate
lambda_e = 1000;

%% =========================================================
% Convergence Criteria
% ==========================================================

% Pointwise tolerance used only for settling-time metric
tol = 2e-2;

% Robust final-window RMS tolerances
tol_q    = 1e-2;
tol_qdot = 1e-2;

% Final window used to assess robust convergence
convergence_window = 0.5;

N_convergence = round(convergence_window/dt);

if N_convergence < 1 || N_convergence > Nt
    error('Invalid convergence window.');
end

%% =========================================================
% Disturbance Parameters
% =========================================================
%
% d(t) =
%     A1*sin(2*pi*f1*t)
%   + A2*sin(2*pi*f2*t)
%   + d_colored(t)
%
% The colored component is generated as a first-order
% Ornstein-Uhlenbeck process.
% ==========================================================

d_A1    = 0.05;
d_f1    = 0.5;

d_A2    = 0.025;
d_f2    = 1.2;

d_sigma = 0.02;
d_tau   = 0.15;

%% =========================================================
% Nominal Linearized Model
% ==========================================================

A = [
    0, 1, 0, 0, 0;
   -Kv/M, -Cv/M, ...
    2*alpha_B*phi_p_0/M, 0, ...
    2*alpha_E*Ve_0/M;
    0, 0, -gamma_Phi, 0, 0;
    0, 0, ...
    2*k_Phi*phi_p_0/C_th, ...
    -h/C_th, ...
    2*k_E*Ve_0/C_th;
    0, 0, 0, 0, -1/tau_e
];

B = [
    0;
    0;
    0;
    0;
    1/tau_e
];

B_d = [
    0;
    1/M;
    0;
    0;
    0
];

C = [
    1 0 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1
];

%% =========================================================
% LQG Controller and Common LQE
% ==========================================================

Q = diag([100 40 10 10 1]);
R = 0.25;

K = lqr(A,B,Q,R);

% Same observer/noise assumptions for both controllers
W = 0.1*eye(5);
V = 0.1*eye(4);

% Common nominal LQE observer
L = [1.7923    0.0490    0.0005    0.0001;
    1.1063    0.1373    0.0010    0.0011;
    0.0049    0.8610    0.0049   -0.0000;
    0.0005    0.0487    0.1538    0.0003;
    0.0001   -0.0000    0.0003    0.0499];

%= lqe(A,eye(5),C,W,V);

%% =========================================================
% SMC Parameters
% ==========================================================

Ks  = 2.3;
eta = 5;
phi = 0.16;

%% =========================================================
% Storage
% ==========================================================

LQG_RMS        = zeros(N,1);
LQG_energy     = zeros(N,1);
LQG_Vmax       = zeros(N,1);
LQG_saturation = zeros(N,1);
LQG_success    = false(N,1);
LQG_settle     = nan(N,1);

LQG_qRMS_final    = zeros(N,1);
LQG_qdotRMS_final = zeros(N,1);

SMC_RMS        = zeros(N,1);
SMC_energy     = zeros(N,1);
SMC_Vmax       = zeros(N,1);
SMC_saturation = zeros(N,1);
SMC_success    = false(N,1);
SMC_settle     = nan(N,1);

SMC_qRMS_final    = zeros(N,1);
SMC_qdotRMS_final = zeros(N,1);

%% ---------------------------------------------------------
% Phase-plane trajectory storage
% ---------------------------------------------------------

N_phase = 20;

% Select trials distributed across the entire Monte Carlo ensemble
phase_trials = round(linspace(1,N,N_phase));

% Uncontrolled
UNC_q    = zeros(Nt,N_phase);
UNC_qdot = zeros(Nt,N_phase);

% LQG using true state
LQG_true_q    = zeros(Nt,N_phase);
LQG_true_qdot = zeros(Nt,N_phase);

% LQG using estimated state
LQG_est_q    = zeros(Nt,N_phase);
LQG_est_qdot = zeros(Nt,N_phase);

%% =========================================================
% Monte Carlo Robustness Simulation
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('LQG vs Observer-Based SMC Robustness Study\n');
fprintf('=========================================================\n');
fprintf('Monte Carlo realizations: %d\n',N);
fprintf('Simulation duration: %.2f s\n',Tsim);
fprintf('Sampling interval: %.4f s\n',dt);
fprintf('Total controller simulations: %d\n',2*N);
fprintf('Final-window convergence interval: %.2f s\n', ...
    convergence_window);
fprintf('q RMS tolerance: %.4f\n',tol_q);
fprintf('qdot RMS tolerance: %.4f\n',tol_qdot);
fprintf('=========================================================\n\n');

%% =========================================================
% Monte Carlo Loop
% ==========================================================

for i = 1:N

    %% -----------------------------------------------------
    % Random uncertain physical parameters
    % ------------------------------------------------------

    M_i = M*(1 + 0.80*(2*rand-1));

    Cv_i = Cv*(1 + 1.00*(2*rand-1));

    Kv_i = Kv*(1 + 0.80*(2*rand-1));

    alpha_B_i = ...
        alpha_B*(1 + 0.70*(2*rand-1));

    alpha_E_i = ...
        alpha_E*(1 + 0.70*(2*rand-1));

    gamma_i = ...
        gamma_Phi*(1 + 0.70*(2*rand-1));

    kPhi_i = ...
        k_Phi*(1 + 0.70*(2*rand-1));

    kE_i = ...
        k_E*(1 + 0.70*(2*rand-1));

    h_i = ...
        h*(1 + 0.70*(2*rand-1));

    tau_i = ...
        tau_e*(1 + 0.70*(2*rand-1));

    %% -----------------------------------------------------
    % Uncertain plant
    % ------------------------------------------------------

    A_i = [
        0, 1, 0, 0, 0;
       -Kv_i/M_i, -Cv_i/M_i, ...
        2*alpha_B_i*phi_p_0/M_i, 0, ...
        2*alpha_E_i*Ve_0/M_i;
        0, 0, -gamma_i, 0, 0;
        0, 0, ...
        2*kPhi_i*phi_p_0/C_th, ...
        -h_i/C_th, ...
        2*kE_i*Ve_0/C_th;
        0, 0, 0, 0, -1/tau_i
    ];

    B_i = [
        0;
        0;
        0;
        0;
        1/tau_i
    ];

    %% -----------------------------------------------------
    % Initial condition
    % ------------------------------------------------------

    dx0 = ICs(i,:).';

    %% -----------------------------------------------------
    % Generate ONE disturbance realization.
    %
    % The exact same disturbance is used by LQG and SMC.
    % ------------------------------------------------------

    d_profile = generate_disturbance( ...
        t, ...
        dt, ...
        d_A1, ...
        d_f1, ...
        d_A2, ...
        d_f2, ...
        d_sigma, ...
        d_tau);

    %% =====================================================
    % LQG
    % ======================================================

    x    = dx0;
    xhat = zeros(5,1);

    X       = zeros(Nt,1);
    Xdot    = zeros(Nt,1);
    U       = zeros(Nt,1);
    Ve_hist = zeros(Nt,1);

    sat_count = 0;

    for k = 1:Nt

        %% Measurement

        y = C*x;

        %% External disturbance

        d = d_profile(k);

        %% LQG perturbation-voltage command

        dVe_cmd = -K*xhat;

        %% Convert perturbation command to physical voltage

        Ve_cmd = Ve_0 + dVe_cmd;

        %% Saturation statistics BEFORE clipping

        if Ve_cmd < Ve_min || Ve_cmd > Ve_max
            sat_count = sat_count + 1;
        end

        %% Apply physical actuator saturation

        Ve_cmd = max(Ve_min,min(Ve_max,Ve_cmd));

        %% Convert saturated physical command
        % to perturbation coordinate

        u = Ve_cmd - Ve_0;

        %% Disturbance input matrix

        B_d_i = [
            0;
            1/M_i;
            0;
            0;
            0
        ];

        %% Uncertain plant

        xdot = ...
            A_i*x + ...
            B_i*u + ...
            B_d_i*d;

        x = x + dt*xdot;

        %% Common nominal LQE observer

        xhat_dot = ...
            A*xhat + ...
            B*u + ...
            L*(y - C*xhat);

        xhat = xhat + dt*xhat_dot;

        %% Actual physical actuator state

        Ve = Ve_0 + x(5);

        %% Store

        X(k)       = x(1);
        Xdot(k)    = x(2);
        U(k)       = u;
        Ve_hist(k) = Ve;

    end

    %% -----------------------------------------------------
    % LQG metrics
    % ------------------------------------------------------

    LQG_RMS(i) = rms(X - q_ref);

    LQG_energy(i) = trapz(t,U.^2);

    LQG_Vmax(i) = max(Ve_hist);

    LQG_saturation(i) = sat_count/Nt;

    %% -----------------------------------------------------
    % LQG final-window convergence
    % ------------------------------------------------------

    idx_start = max(1,Nt-N_convergence+1);

    q_tail = X(idx_start:end);

    qdot_tail = Xdot(idx_start:end);

    q_rms_tail = rms(q_tail - q_ref);

    qdot_rms_tail = rms(qdot_tail);

    LQG_qRMS_final(i) = q_rms_tail;

    LQG_qdotRMS_final(i) = qdot_rms_tail;

    LQG_success(i) = ...
        (q_rms_tail < tol_q) && ...
        (qdot_rms_tail < tol_qdot);

    %% -----------------------------------------------------
    % LQG settling-time metric
    % ------------------------------------------------------

    err = abs(X - q_ref);

    idx = find(arrayfun( ...
        @(k) all(err(k:end) < tol), ...
        1:Nt),1);

    if ~isempty(idx)
        LQG_settle(i) = t(idx);
    end

    %% =====================================================
    % Observer-Based SMC
    % ======================================================

    x    = dx0;
    xhat = zeros(5,1);

    X       = zeros(Nt,1);
    Xdot    = zeros(Nt,1);
    U       = zeros(Nt,1);
    Ve_hist = zeros(Nt,1);

    %% -----------------------------------------------------
    % Temporarily save LQG trajectory for phase-plane plot
    % ------------------------------------------------------

    LQG_phase_q_trial    = X;
    LQG_phase_qdot_trial = Xdot;

    sat_count = 0;

    for k = 1:Nt

        %% Measurement

        y = C*x;

        %% External disturbance
        %
        % SAME disturbance realization used by LQG.
        %

        d = d_profile(k);

        %% Estimated states

        qhat    = xhat(1);
        qdothat = xhat(2);
        Phihat  = xhat(3);

        %% Sliding surface

        error_hat = qhat - q_ref;

        s = qdothat + eta*error_hat;

        %% Boundary-layer saturation

        sat_s = max(-1,min(1,s/phi));

        %% Switching term

        usw = Ks*sat_s;

        %% Desired acceleration

        qddot_des = ...
            -eta*qdothat ...
            -usw;

        %% Desired physical voltage
        %
        % Controller uses nominal plant parameters here.
        %

        dVe_des = ...
            (M/(2*alpha_E*Ve_0)) * ( ...
                  qddot_des ...
                + (Cv/M)*qdothat ...
                + (Kv/M)*qhat ...
                - (2*alpha_B*phi_p_0/M)*Phihat );

        Ve_des = Ve_0 + dVe_des;

        %% Current physical actuator state

        Ve = Ve_0 + x(5);

        %% Physical actuator command BEFORE saturation

        Ve_cmd = ...
            Ve + tau_i*lambda_e*(Ve_des - Ve);

        %% Saturation statistics BEFORE clipping

        if Ve_cmd < Ve_min || Ve_cmd > Ve_max
            sat_count = sat_count + 1;
        end

        %% Apply physical actuator saturation

        Ve_cmd = max(Ve_min,min(Ve_max,Ve_cmd));

        %% Convert physical command to perturbation coordinate

        u = Ve_cmd - Ve_0;

        %% Uncertain plant

        xdot = ...
            A_i*x + ...
            B_i*u + ...
            B_d_i*d;

        x = x + dt*xdot;

        %% Common nominal LQE observer

        xhat_dot = ...
            A*xhat + ...
            B*u + ...
            L*(y - C*xhat);

        xhat = xhat + dt*xhat_dot;

        %% Actual physical actuator state

        Ve = Ve_0 + x(5);

        %% Store

        X(k)       = x(1);
        Xdot(k)    = x(2);
        U(k)       = u;
        Ve_hist(k) = Ve;

    end

    %% -----------------------------------------------------
    % SMC metrics
    % ------------------------------------------------------

    SMC_RMS(i) = rms(X - q_ref);

    SMC_energy(i) = trapz(t,U.^2);

    SMC_Vmax(i) = max(Ve_hist);

    SMC_saturation(i) = sat_count/Nt;

    %% -----------------------------------------------------
    % SMC final-window convergence
    % ------------------------------------------------------

    q_tail = X(idx_start:end);

    qdot_tail = Xdot(idx_start:end);

    q_rms_tail = rms(q_tail - q_ref);

    qdot_rms_tail = rms(qdot_tail);

    SMC_qRMS_final(i) = q_rms_tail;

    SMC_qdotRMS_final(i) = qdot_rms_tail;

    SMC_success(i) = ...
        (q_rms_tail < tol_q) && ...
        (qdot_rms_tail < tol_qdot);

    %% -----------------------------------------------------
    % SMC settling-time metric
    % ------------------------------------------------------

    err = abs(X - q_ref);

    idx = find(arrayfun( ...
        @(k) all(err(k:end) < tol), ...
        1:Nt),1);

    if ~isempty(idx)
        SMC_settle(i) = t(idx);
    end

    %% =====================================================
    % Phase-plane data for selected Monte Carlo trials
    % ======================================================
    
    phase_idx = find(phase_trials == i,1);
    
    if ~isempty(phase_idx)
    
        %% -------------------------------------------------
        % LQG — Estimated State
        % -------------------------------------------------
        
        LQG_est_q(:,phase_idx)    = X;
        LQG_est_qdot(:,phase_idx) = Xdot;
    
    
        %% =================================================
        % LQG — True State
        % =================================================
        
        x_true = dx0;
    
        X_true    = zeros(Nt,1);
        Xdot_true = zeros(Nt,1);
    
        for k_true = 1:Nt
    
            % Same disturbance realization
            d_true = d_profile(k_true);
    
            % Ideal full-state feedback
            dVe_true = -K*x_true;
    
            % Physical voltage command
            Ve_cmd_true = Ve_0 + dVe_true;
    
            % Apply same physical actuator saturation
            Ve_cmd_true = ...
                max(Ve_min,min(Ve_max,Ve_cmd_true));
    
            % Convert to perturbation coordinate
            u_true = Ve_cmd_true - Ve_0;
    
            % Plant dynamics
            xdot_true = ...
                A_i*x_true + ...
                B_i*u_true + ...
                B_d_i*d_true;
    
            % Integrate
            x_true = x_true + dt*xdot_true;
    
            % Store
            X_true(k_true)    = x_true(1);
            Xdot_true(k_true) = x_true(2);
    
        end
    
        LQG_true_q(:,phase_idx)    = X_true;
        LQG_true_qdot(:,phase_idx) = Xdot_true;
    
    
        %% =================================================
        % Uncontrolled Plant
        % =================================================
        
        x_unc = dx0;
    
        X_unc    = zeros(Nt,1);
        Xdot_unc = zeros(Nt,1);
    
        for k_unc = 1:Nt
    
            % Same disturbance realization
            d_unc = d_profile(k_unc);
    
            % No control input
            u_unc = 0;
    
            % Plant dynamics
            xdot_unc = ...
                A_i*x_unc + ...
                B_i*u_unc + ...
                B_d_i*d_unc;
    
            % Integrate
            x_unc = x_unc + dt*xdot_unc;
    
            % Store
            X_unc(k_unc)    = x_unc(1);
            Xdot_unc(k_unc) = x_unc(2);
    
        end
    
        UNC_q(:,phase_idx)    = X_unc;
        UNC_qdot(:,phase_idx) = Xdot_unc;
    
    end

    %% =====================================================
    % Trial diagnostic
    % ======================================================

    fprintf(['Trial %4d/%4d | ', ...
             'LQG: %d  SMC: %d | ', ...
             'LQG qRMS = %.4f, qdotRMS = %.4f | ', ...
             'SMC qRMS = %.4f, qdotRMS = %.4f\n'], ...
             i,N, ...
             LQG_success(i), ...
             SMC_success(i), ...
             LQG_qRMS_final(i), ...
             LQG_qdotRMS_final(i), ...
             SMC_qRMS_final(i), ...
             SMC_qdotRMS_final(i));

end

%% =========================================================
% LQG vs SMC Comparison
% ==========================================================

comparison_table = table( ...
    [mean(LQG_success)*100; ...
     mean(SMC_success)*100], ...
    [mean(LQG_energy); ...
     mean(SMC_energy)], ...
    [mean(LQG_Vmax); ...
     mean(SMC_Vmax)], ...
    [mean(LQG_saturation)*100; ...
     mean(SMC_saturation)*100], ...
    'VariableNames', ...
    {'SuccessRate', ...
     'MeanEnergy', ...
     'MeanMaxVoltage', ...
     'MeanSaturationTimePercent'}, ...
    'RowNames', ...
    {'LQG','SMC'});

%% =========================================================
% Print Results
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('LQG vs Observer-Based SMC\n');
fprintf('=========================================================\n\n');

disp(comparison_table);

%% =========================================================
% Detailed Success Statistics
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('SUCCESS STATISTICS\n');
fprintf('=========================================================\n');

fprintf('\nLQG:\n');
fprintf('  Successful trials: %d / %d\n', ...
    sum(LQG_success),N);

fprintf('  Success rate: %.2f %%\n', ...
    100*mean(LQG_success));

fprintf('  Mean final q RMS: %.6f\n', ...
    mean(LQG_qRMS_final));

fprintf('  Mean final qdot RMS: %.6f\n', ...
    mean(LQG_qdotRMS_final));

fprintf('\nSMC:\n');
fprintf('  Successful trials: %d / %d\n', ...
    sum(SMC_success),N);

fprintf('  Success rate: %.2f %%\n', ...
    100*mean(SMC_success));

fprintf('  Mean final q RMS: %.6f\n', ...
    mean(SMC_qRMS_final));

fprintf('  Mean final qdot RMS: %.6f\n', ...
    mean(SMC_qdotRMS_final));

%% =========================================================
% Settling-Time Statistics
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('SETTLING-TIME STATISTICS\n');
fprintf('=========================================================\n');

fprintf('\nMean LQG settling time: %.3f s\n', ...
    mean(LQG_settle,'omitnan'));

fprintf('Mean SMC settling time: %.3f s\n', ...
    mean(SMC_settle,'omitnan'));

fprintf('\nLQG trials with finite settling time: %d / %d\n', ...
    sum(~isnan(LQG_settle)),N);

fprintf('SMC trials with finite settling time: %d / %d\n', ...
    sum(~isnan(SMC_settle)),N);

%% =========================================================
% Saturation Statistics
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('SATURATION STATISTICS\n');
fprintf('=========================================================\n');

fprintf('\nLQG saturation incidence: %.2f %%\n', ...
    100*mean(LQG_saturation));

fprintf('SMC saturation incidence: %.2f %%\n', ...
    100*mean(SMC_saturation));

fprintf('\nLQG maximum observed physical voltage: %.6f\n', ...
    max(LQG_Vmax));

fprintf('SMC maximum observed physical voltage: %.6f\n', ...
    max(SMC_Vmax));

%% =========================================================
% Results — LQG vs Observer-Based SMC
% ==========================================================

controller = categorical([
    repmat("LQG",N,1);
    repmat("SMC",N,1)
]);

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ---------------------------------------------------------
% RMS Regulation Error
% ---------------------------------------------------------

figure;

boxchart(controller,[LQG_RMS; SMC_RMS]);

ylabel('RMS Error');
title('');

grid off;
set(gca,'FontSize',12);

%% ---------------------------------------------------------
% Control Energy
% ---------------------------------------------------------

figure;

boxchart(controller,[LQG_energy; SMC_energy]);

set(gca,'YScale','log');

ylabel('Energy');

grid off;
set(gca,'FontSize',12);

%% ---------------------------------------------------------
% Settling Time
% ---------------------------------------------------------

valid_LQG_settle = ...
    LQG_settle(~isnan(LQG_settle));

valid_SMC_settle = ...
    SMC_settle(~isnan(SMC_settle));

controller_settle = categorical([
    repmat("LQG",length(valid_LQG_settle),1);
    repmat("SMC",length(valid_SMC_settle),1)
]);

figure;

boxchart(controller_settle, ...
    [valid_LQG_settle; valid_SMC_settle]);

ylabel('Time (s)');

ylim([0 Tsim]);

grid off;
set(gca,'FontSize',12);

%% =========================================================
% Phase-Plane Comparison — Multiple Monte Carlo Trials
% ==========================================================

figure;
hold on;

%% ---------------------------------------------------------
% Model colors
% ---------------------------------------------------------

color_unc  = [0.00 0.00 0.00];   % Black
color_true = [0.00 0.45 0.74];   % Blue
color_est  = [0.85 0.33 0.10];   % Orange


%% ---------------------------------------------------------
% Uncontrolled trajectories
% ---------------------------------------------------------

plot(UNC_q, ...
     UNC_qdot, ...
     '--', ...
     'Color',color_unc, ...
     'LineWidth',0.8);


%% ---------------------------------------------------------
% LQG — True State trajectories
% ---------------------------------------------------------

plot(LQG_true_q, ...
     LQG_true_qdot, ...
     '-', ...
     'Color',color_true, ...
     'LineWidth',1.0);


%% ---------------------------------------------------------
% LQG — Estimated State trajectories
% ---------------------------------------------------------

plot(LQG_est_q, ...
     LQG_est_qdot, ...
     '-.', ...
     'Color',color_est, ...
     'LineWidth',1.0);


%% ---------------------------------------------------------
% Create legend handles
% ---------------------------------------------------------

h_unc = plot(nan,nan, ...
    '--', ...
    'Color',color_unc, ...
    'LineWidth',1.2);

h_true = plot(nan,nan, ...
    '-', ...
    'Color',color_true, ...
    'LineWidth',1.2);

h_est = plot(nan,nan, ...
    '-.', ...
    'Color',color_est, ...
    'LineWidth',1.2);


%% ---------------------------------------------------------
% Labels
% ---------------------------------------------------------

xlabel('$q$');
ylabel('$\dot{q}$');

legend([h_unc,h_true,h_est], ...
       {'Uncontrolled', ...
        'LQG -- True State', ...
        'LQG -- Estimated State'}, ...
       'Location','best');


%% ---------------------------------------------------------
% Axes formatting
% ---------------------------------------------------------

grid off;
box on;

set(gca,'FontSize',12);

%% ---------------------------------------------------------
% IEEE figure dimensions
% ---------------------------------------------------------

set(gcf,'Units','inches');
set(gcf,'Position',[1 1 7.16 4.50]);

hold off;

%% =========================================================
% Local Function — Correlated External Disturbance
% ==========================================================

function d = generate_disturbance( ...
    t,dt,A1,f1,A2,f2,sigma_d,tau_d)

    Nt = length(t);

    %% Deterministic low-frequency components

    d_det = ...
        A1*sin(2*pi*f1*t) + ...
        A2*sin(2*pi*f2*t);

    %% First-order correlated stochastic component
    %
    % Continuous-time OU process:
    %
    %   dd_c/dt = -(1/tau_d)*d_c
    %             + sqrt(2*sigma_d^2/tau_d)*w(t)
    %
    % Euler-Maruyama discretization.

    d_colored = zeros(1,Nt);

    noise_gain = ...
        sqrt(2*sigma_d^2/tau_d);

    for k = 2:Nt

        d_colored(k) = ...
            d_colored(k-1) ...
            + dt*(-d_colored(k-1)/tau_d) ...
            + sqrt(dt)*noise_gain*randn;

    end

    %% Total disturbance

    d = d_det + d_colored;

end