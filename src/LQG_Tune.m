%% ----------------------------------------------------------
% Monte Carlo Configuration
% ----------------------------------------------------------

rng(7);

% Number of fixed initial conditions
N_IC = 1000;

Tsim = 10;
tspan = 0:0.01:Tsim;

% Initial-condition perturbation bounds
q_max       = 0.30;
q_dot_max   = 0.15;
dphi_max    = 0.10;
dT_max      = 0.08;
dV_e_max    = 0.15;

% Convergence criterion
eps_q      = 1e-2;
eps_q_dot  = 1e-2;

% Settling criterion
settle_tol = 2e-2;

% Physical actuator-voltage saturation limit
V_max = 15;

%% ----------------------------------------------------------
% Fixed Monte Carlo Initial-Condition Ensemble
% ----------------------------------------------------------
%
% IMPORTANT:
%
% Every controller and observer candidate is evaluated using
% exactly the same 1000 initial conditions.
%
% Therefore, differences between candidates are attributable
% to the controller/observer design rather than differences
% in the Monte Carlo ensemble.
%
% Total simulations:
%
%   Stage 1 = N_K_candidates * N_IC
%   Stage 2 = N_top_K * N_L_candidates * N_IC
%
% where N_IC = 1000.
%

ICs = zeros(N_IC,5);

for i = 1:N_IC

    ICs(i,:) = [
        q_max      * (2*rand - 1), ...
        q_dot_max  * (2*rand - 1), ...
        dphi_max   * (2*rand - 1), ...
        dT_max     * (2*rand - 1), ...
        dV_e_max   * (2*rand - 1)
    ];

end

%% ==========================================================
% STAGE 1: LQR Q/R Controller Sweep
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('STAGE 1: LQR Q/R Controller Sweep\n');
fprintf('=========================================================\n');

%% ----------------------------------------------------------
% Q weighting candidates
% ----------------------------------------------------------
%
% The nominal Q matrix is:
%
%   Q = diag([100 10 5 5 1])
%
% The sweep varies the relative importance of each state.
%

Qq_values   = [25 50 100 200 400];
Qqd_values  = [2.5 5 10 20 40];

Qphi_values = [1 2.5 5 10];
QT_values   = [1 2.5 5 10];
QVe_values  = [0.25 0.5 1 2.5];

% Control penalty
R_values = [0.25 0.5 1 2 4];

%% ----------------------------------------------------------
% Number of Q/R candidates
% ----------------------------------------------------------

N_K_candidates = ...
    length(Qq_values) * ...
    length(Qqd_values) * ...
    length(Qphi_values) * ...
    length(QT_values) * ...
    length(QVe_values) * ...
    length(R_values);

fprintf('LQR Q/R candidates: %d\n',N_K_candidates);
fprintf('Initial conditions per candidate: %d\n',N_IC);
fprintf('Stage 1 nonlinear simulations: %d\n', ...
    N_K_candidates*N_IC);

%% ----------------------------------------------------------
% Fixed nominal Kalman observer for Stage 1
% ----------------------------------------------------------
%
% The observer is held fixed while searching for the best
% state-feedback controller.
%
% This prevents the Q/R search from being confounded by
% simultaneous changes in estimator dynamics.
%

W_nominal = diag([ ...
    1e-3, ...
    1e-2, ...
    1e-4, ...
    1e-4, ...
    1e-3]);

V_nominal = diag([ ...
    1e-3, ...
    1e-4, ...
    1e-3, ...
    1e-3]);

L_nominal = lqe( ...
    A, ...
    eye(5), ...
    C, ...
    W_nominal, ...
    V_nominal);

%% ----------------------------------------------------------
% Stage 1 result storage
% ----------------------------------------------------------

sweep_Q = zeros(N_K_candidates,5);
sweep_R = zeros(N_K_candidates,1);

sweep_success = zeros(N_K_candidates,1);
sweep_RMS = zeros(N_K_candidates,1);
sweep_energy = zeros(N_K_candidates,1);
sweep_settle = zeros(N_K_candidates,1);
sweep_maxVoltage = zeros(N_K_candidates,1);
sweep_saturation = zeros(N_K_candidates,1);
sweep_satDuration = zeros(N_K_candidates,1);

sweep_K = zeros(N_K_candidates,5);

candidate = 0;

%% ==========================================================
% Stage 1 Q/R Sweep
% ==========================================================

for iq = 1:length(Qq_values)

    for iqd = 1:length(Qqd_values)

        for iphi = 1:length(Qphi_values)

            for iT = 1:length(QT_values)

                for iVe = 1:length(QVe_values)

                    for ir = 1:length(R_values)

                        candidate = candidate + 1;

                        %% --------------------------------------------------
                        % Construct Q and R
                        % --------------------------------------------------

                        Q_test = diag([ ...
                            Qq_values(iq), ...
                            Qqd_values(iqd), ...
                            Qphi_values(iphi), ...
                            QT_values(iT), ...
                            QVe_values(iVe)]);

                        R_test = R_values(ir);

                        %% --------------------------------------------------
                        % Compute LQR gain
                        % --------------------------------------------------

                        try

                            K_test = lqr( ...
                                A, ...
                                B, ...
                                Q_test, ...
                                R_test);

                        catch

                            warning( ...
                                'LQR failed for candidate %d.', ...
                                candidate);

                            continue;

                        end

                        %% --------------------------------------------------
                        % Candidate metric storage
                        % --------------------------------------------------

                        RMS_i = zeros(N_IC,1);
                        energy_i = zeros(N_IC,1);
                        settle_i = zeros(N_IC,1);
                        maxVoltage_i = zeros(N_IC,1);

                        saturated_i = false(N_IC,1);
                        saturationDuration_i = zeros(N_IC,1);
                        converged_i = false(N_IC,1);

                        %% ==================================================
                        % Monte Carlo simulations
                        % ==================================================

                        for i = 1:N_IC

                            %% Fixed initial condition

                            dx_0 = ICs(i,:).';

                            x_0 = [
                                dx_0(1);
                                dx_0(2);
                                phi_p_0 + dx_0(3);
                                T_0     + dx_0(4);
                                V_e_0   + dx_0(5)
                            ];

                            %% Observer initialized at nominal equilibrium

                            xhat_0 = zeros(5,1);

                            z_0 = [x_0; xhat_0];

                            %% Nonlinear LQG dynamics

                            f_lqg = @(t,z) ...
                                lqg_nonlinear_dynamics( ...
                                t,z,A,B,C,L_nominal,K_test,...
                                M,C_v,K_v,...
                                phi_p_0,T_0,V_e_0,...
                                alpha_B,alpha_E,...
                                gamma_Phi,k_Phi,k_E,...
                                C_th,h,tau_e,V_max);

                            [t,z] = ode45( ...
                                f_lqg, ...
                                tspan, ...
                                z_0);

                            %% Extract states

                            x = z(:,1:5);
                            x_hat = z(:,6:10);

                            q = x(:,1);
                            q_dot = x(:,2);

                            %% --------------------------------------------------
                            % Reconstruct physical control input
                            % --------------------------------------------------

                            u_cmd = -x_hat*K_test.';

                            u = max( ...
                                min(u_cmd,V_max), ...
                                -V_max);

                            %% --------------------------------------------------
                            % Performance metrics
                            % --------------------------------------------------

                            RMS_i(i) = sqrt(mean(q.^2));

                            energy_i(i) = ...
                                trapz(t,u.^2);

                            % Physical actuator voltage
                            maxVoltage_i(i) = ...
                                max(abs(u));

                            %% --------------------------------------------------
                            % Saturation
                            % --------------------------------------------------

                            saturated = ...
                                abs(u) >= V_max - 1e-10;

                            saturated_i(i) = ...
                                any(saturated);

                            saturationDuration_i(i) = ...
                                100 * ...
                                trapz( ...
                                t, ...
                                double(saturated)) / Tsim;

                            %% --------------------------------------------------
                            % Convergence
                            % --------------------------------------------------

                            converged_i(i) = ...
                                abs(q(end)) < eps_q && ...
                                abs(q_dot(end)) < eps_q_dot;

                            %% --------------------------------------------------
                            % Settling time
                            % --------------------------------------------------

                            err_norm = ...
                                sqrt(q.^2 + q_dot.^2);

                            idx_settle = ...
                                find( ...
                                err_norm < settle_tol, ...
                                1, ...
                                'first');

                            if isempty(idx_settle)

                                settle_i(i) = Tsim;

                            else

                                remaining_err = ...
                                    err_norm(idx_settle:end);

                                if all(remaining_err < settle_tol)

                                    settle_i(i) = ...
                                        t(idx_settle);

                                else

                                    settle_i(i) = Tsim;

                                end

                            end

                        end

                        %% ==================================================
                        % Candidate summary
                        % ==================================================

                        sweep_Q(candidate,:) = ...
                            diag(Q_test).';

                        sweep_R(candidate) = ...
                            R_test;

                        sweep_K(candidate,:) = ...
                            K_test;

                        sweep_success(candidate) = ...
                            100*mean(converged_i);

                        sweep_RMS(candidate) = ...
                            mean(RMS_i);

                        sweep_energy(candidate) = ...
                            mean(energy_i);

                        sweep_settle(candidate) = ...
                            mean(settle_i);

                        sweep_maxVoltage(candidate) = ...
                            mean(maxVoltage_i);

                        sweep_saturation(candidate) = ...
                            100*mean(saturated_i);

                        sweep_satDuration(candidate) = ...
                            mean(saturationDuration_i);

                        %% --------------------------------------------------
                        % Progress
                        % --------------------------------------------------

                        fprintf( ...
                            'K candidate %4d/%4d | ', ...
                            candidate,N_K_candidates);

                        fprintf( ...
                            'Q = [%.2g %.2g %.2g %.2g %.2g] ', ...
                            diag(Q_test));

                        fprintf( ...
                            '| R = %.2g | Success = %.1f %% | RMS = %.4g\n', ...
                            R_test,...
                            sweep_success(candidate),...
                            sweep_RMS(candidate));

                    end

                end

            end

        end

    end

end

%% ==========================================================
% Stage 1 Results
% ==========================================================

Qq_result   = sweep_Q(:,1);
Qqd_result  = sweep_Q(:,2);
Qphi_result = sweep_Q(:,3);
QT_result   = sweep_Q(:,4);
QVe_result  = sweep_Q(:,5);

K1 = sweep_K(:,1);
K2 = sweep_K(:,2);
K3 = sweep_K(:,3);
K4 = sweep_K(:,4);
K5 = sweep_K(:,5);

K_results_table = table( ...
    Qq_result, ...
    Qqd_result, ...
    Qphi_result, ...
    QT_result, ...
    QVe_result, ...
    sweep_R, ...
    K1, K2, K3, K4, K5, ...
    sweep_success, ...
    sweep_RMS, ...
    sweep_settle, ...
    sweep_energy, ...
    sweep_maxVoltage, ...
    sweep_saturation, ...
    sweep_satDuration, ...
    'VariableNames',{ ...
    'Q_q', ...
    'Q_qdot', ...
    'Q_Phi', ...
    'Q_T', ...
    'Q_Ve', ...
    'R', ...
    'K_q', ...
    'K_qdot', ...
    'K_Phi', ...
    'K_T', ...
    'K_Ve', ...
    'SuccessRate', ...
    'MeanRMS', ...
    'MeanSettlingTime', ...
    'MeanEnergy', ...
    'MeanMaxVoltage', ...
    'SaturationIncidence', ...
    'MeanSaturationDuration'});

%% ----------------------------------------------------------
% Sort Stage 1 results
% ----------------------------------------------------------

K_results_table = sortrows( ...
    K_results_table, ...
    {'SuccessRate','MeanRMS','MeanEnergy'}, ...
    {'descend','ascend','ascend'});

fprintf('\n');
fprintf('=========================================================\n');
fprintf('TOP LQR CONTROLLERS\n');
fprintf('=========================================================\n');

disp(K_results_table( ...
    1:min(20,height(K_results_table)),:));

%% ==========================================================
% Select Top LQR Controllers for LQE Sweep
% ==========================================================
%
% Rather than performing a prohibitively large joint sweep of
% every Q/R/W/V combination, retain the strongest LQR
% controllers and optimize the observer around them.
%

N_top_K = min(10,height(K_results_table));

top_K_table = K_results_table(1:N_top_K,:);

fprintf('\n');
fprintf('Retaining top %d LQR controllers for LQE sweep.\n', ...
    N_top_K);

%% ==========================================================
% STAGE 2: LQE W/V Observer Sweep
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('STAGE 2: LQE W/V Observer Sweep\n');
fprintf('=========================================================\n');

%% ----------------------------------------------------------
% Process-noise covariance scaling
% ----------------------------------------------------------
%
% W = W_nominal * W_scale
%
% Larger W_scale tells the Kalman filter that the plant model
% is less trustworthy and that process disturbances are larger.
%

W_scale_values = [0.1 1 10];

%% ----------------------------------------------------------
% Measurement-noise covariance scaling
% ----------------------------------------------------------
%
% V = V_nominal * V_scale
%
% Larger V_scale tells the Kalman filter that measurements
% are noisier and should therefore be trusted less.
%

V_scale_values = [0.1 1 10];

%% ----------------------------------------------------------
% Number of observer candidates
% ----------------------------------------------------------

N_L_candidates = ...
    length(W_scale_values) * ...
    length(V_scale_values);

fprintf('LQE W/V candidates per controller: %d\n', ...
    N_L_candidates);

fprintf('Top LQR controllers evaluated: %d\n', ...
    N_top_K);

fprintf('Initial conditions per LQG candidate: %d\n', ...
    N_IC);

fprintf('Stage 2 nonlinear simulations: %d\n', ...
    N_top_K*N_L_candidates*N_IC);

%% ----------------------------------------------------------
% Result storage
% ----------------------------------------------------------

N_LQG_candidates = ...
    N_top_K * N_L_candidates;

lqg_K_index = zeros(N_LQG_candidates,1);

lqg_W_scale = zeros(N_LQG_candidates,1);
lqg_V_scale = zeros(N_LQG_candidates,1);

lqg_success = zeros(N_LQG_candidates,1);
lqg_RMS = zeros(N_LQG_candidates,1);
lqg_energy = zeros(N_LQG_candidates,1);
lqg_settle = zeros(N_LQG_candidates,1);
lqg_maxVoltage = zeros(N_LQG_candidates,1);
lqg_saturation = zeros(N_LQG_candidates,1);
lqg_satDuration = zeros(N_LQG_candidates,1);

lqg_L = zeros(N_LQG_candidates,25);

candidate = 0;

%% ==========================================================
% Stage 2 W/V Sweep
% ==========================================================

for k_index = 1:N_top_K

    %% ------------------------------------------------------
    % Recover LQR controller
    % ------------------------------------------------------

    Q_test = diag([ ...
        top_K_table.Q_q(k_index), ...
        top_K_table.Q_qdot(k_index), ...
        top_K_table.Q_Phi(k_index), ...
        top_K_table.Q_T(k_index), ...
        top_K_table.Q_Ve(k_index)]);

    R_test = top_K_table.R(k_index);

    K_test = lqr(A,B,Q_test,R_test);

    for iW = 1:length(W_scale_values)

        for iV = 1:length(V_scale_values)

            candidate = candidate + 1;

            W_scale = W_scale_values(iW);
            V_scale = V_scale_values(iV);

            %% --------------------------------------------------
            % Construct W and V
            % --------------------------------------------------

            W_test = W_scale * W_nominal;
            V_test = V_scale * V_nominal;

            %% --------------------------------------------------
            % Compute Kalman gain
            % --------------------------------------------------

            try

                L_test = lqe( ...
                    A, ...
                    eye(5), ...
                    C, ...
                    W_test, ...
                    V_test);

            catch

                warning( ...
                    'LQE failed for candidate %d.', ...
                    candidate);

                continue;

            end

            %% --------------------------------------------------
            % Candidate metric storage
            % --------------------------------------------------

            RMS_i = zeros(N_IC,1);
            energy_i = zeros(N_IC,1);
            settle_i = zeros(N_IC,1);
            maxVoltage_i = zeros(N_IC,1);

            saturated_i = false(N_IC,1);
            saturationDuration_i = zeros(N_IC,1);
            converged_i = false(N_IC,1);

            %% ==================================================
            % Monte Carlo simulations
            % ==================================================

            for i = 1:N_IC

                %% Fixed initial condition

                dx_0 = ICs(i,:).';

                x_0 = [
                    dx_0(1);
                    dx_0(2);
                    phi_p_0 + dx_0(3);
                    T_0     + dx_0(4);
                    V_e_0   + dx_0(5)
                ];

                %% Observer initialized at nominal equilibrium

                xhat_0 = zeros(5,1);

                z_0 = [x_0; xhat_0];

                %% Nonlinear LQG dynamics

                f_lqg = @(t,z) ...
                    lqg_nonlinear_dynamics( ...
                    t,z,A,B,C,L_test,K_test,...
                    M,C_v,K_v,...
                    phi_p_0,T_0,V_e_0,...
                    alpha_B,alpha_E,...
                    gamma_Phi,k_Phi,k_E,...
                    C_th,h,tau_e,V_max);

                [t,z] = ode45( ...
                    f_lqg, ...
                    tspan, ...
                    z_0);

                %% Extract states

                x = z(:,1:5);
                x_hat = z(:,6:10);

                q = x(:,1);
                q_dot = x(:,2);

                %% --------------------------------------------------
                % Reconstruct control
                % --------------------------------------------------

                u_cmd = -x_hat*K_test.';

                u = max( ...
                    min(u_cmd,V_max), ...
                    -V_max);

                %% --------------------------------------------------
                % Performance metrics
                % --------------------------------------------------

                RMS_i(i) = ...
                    sqrt(mean(q.^2));

                energy_i(i) = ...
                    trapz(t,u.^2);

                maxVoltage_i(i) = ...
                    max(abs(u));

                %% --------------------------------------------------
                % Saturation
                % --------------------------------------------------

                saturated = ...
                    abs(u) >= V_max - 1e-10;

                saturated_i(i) = ...
                    any(saturated);

                saturationDuration_i(i) = ...
                    100 * ...
                    trapz( ...
                    t, ...
                    double(saturated))/Tsim;

                %% --------------------------------------------------
                % Convergence
                % --------------------------------------------------

                converged_i(i) = ...
                    abs(q(end)) < eps_q && ...
                    abs(q_dot(end)) < eps_q_dot;

                %% --------------------------------------------------
                % Settling time
                % --------------------------------------------------

                err_norm = ...
                    sqrt(q.^2 + q_dot.^2);

                idx_settle = ...
                    find( ...
                    err_norm < settle_tol, ...
                    1, ...
                    'first');

                if isempty(idx_settle)

                    settle_i(i) = Tsim;

                else

                    remaining_err = ...
                        err_norm(idx_settle:end);

                    if all(remaining_err < settle_tol)

                        settle_i(i) = ...
                            t(idx_settle);

                    else

                        settle_i(i) = Tsim;

                    end

                end

            end

            %% ==================================================
            % Candidate summary
            % ==================================================

            lqg_K_index(candidate) = k_index;

            lqg_W_scale(candidate) = W_scale;

            lqg_V_scale(candidate) = V_scale;

            lqg_L(candidate,:) = ...
                L_test(:).';

            lqg_success(candidate) = ...
                100*mean(converged_i);

            lqg_RMS(candidate) = ...
                mean(RMS_i);

            lqg_energy(candidate) = ...
                mean(energy_i);

            lqg_settle(candidate) = ...
                mean(settle_i);

            lqg_maxVoltage(candidate) = ...
                mean(maxVoltage_i);

            lqg_saturation(candidate) = ...
                100*mean(saturated_i);

            lqg_satDuration(candidate) = ...
                mean(saturationDuration_i);

            %% --------------------------------------------------
            % Progress
            % --------------------------------------------------

            fprintf( ...
                'LQG candidate %3d/%3d | K #%2d | ', ...
                candidate, ...
                N_LQG_candidates, ...
                k_index);

            fprintf( ...
                'W = %.2g | V = %.2g | ', ...
                W_scale,V_scale);

            fprintf( ...
                'Success = %.1f %% | RMS = %.4g\n', ...
                lqg_success(candidate), ...
                lqg_RMS(candidate));

        end

    end

end

%% ==========================================================
% Stage 2 LQG Results Table
% ==========================================================

Qq_result = zeros(N_LQG_candidates,1);
Qqd_result = zeros(N_LQG_candidates,1);
Qphi_result = zeros(N_LQG_candidates,1);
QT_result = zeros(N_LQG_candidates,1);
QVe_result = zeros(N_LQG_candidates,1);
R_result = zeros(N_LQG_candidates,1);

for i = 1:N_LQG_candidates

    k_index = lqg_K_index(i);

    Qq_result(i) = top_K_table.Q_q(k_index);
    Qqd_result(i) = top_K_table.Q_qdot(k_index);
    Qphi_result(i) = top_K_table.Q_Phi(k_index);
    QT_result(i) = top_K_table.Q_T(k_index);
    QVe_result(i) = top_K_table.Q_Ve(k_index);
    R_result(i) = top_K_table.R(k_index);

end

%% ----------------------------------------------------------
% Create LQG results table
% ----------------------------------------------------------

LQG_results_table = table( ...
    lqg_K_index, ...
    Qq_result, ...
    Qqd_result, ...
    Qphi_result, ...
    QT_result, ...
    QVe_result, ...
    R_result, ...
    lqg_W_scale, ...
    lqg_V_scale, ...
    lqg_success, ...
    lqg_RMS, ...
    lqg_settle, ...
    lqg_energy, ...
    lqg_maxVoltage, ...
    lqg_saturation, ...
    lqg_satDuration, ...
    'VariableNames',{ ...
    'KIndex', ...
    'Q_q', ...
    'Q_qdot', ...
    'Q_Phi', ...
    'Q_T', ...
    'Q_Ve', ...
    'R', ...
    'W_Scale', ...
    'V_Scale', ...
    'SuccessRate', ...
    'MeanRMS', ...
    'MeanSettlingTime', ...
    'MeanEnergy', ...
    'MeanMaxVoltage', ...
    'SaturationIncidence', ...
    'MeanSaturationDuration'});

%% ==========================================================
% Sort LQG Results
% ==========================================================
%
% Primary criterion:
%   Maximum convergence rate
%
% Secondary criteria:
%   RMS regulation error
%   Control energy
%   Saturation incidence
%

LQG_results_table = sortrows( ...
    LQG_results_table, ...
    {'SuccessRate', ...
     'MeanRMS', ...
     'MeanEnergy', ...
     'SaturationIncidence'}, ...
    {'descend', ...
     'ascend', ...
     'ascend', ...
     'ascend'});

fprintf('\n');
fprintf('=========================================================\n');
fprintf('TOP LQG CONTROLLER/OBSERVER COMBINATIONS\n');
fprintf('=========================================================\n');

disp(LQG_results_table( ...
    1:min(20,height(LQG_results_table)),:));

%% ==========================================================
% Select Best LQG Configuration
% ==========================================================

% Maximum observed convergence rate
maxSuccess = ...
    max(LQG_results_table.SuccessRate);

% Controllers within 1 percentage point of maximum
success_margin = 1.0;

eligible = ...
    LQG_results_table.SuccessRate >= ...
    maxSuccess - success_margin;

eligible_table = ...
    LQG_results_table(eligible,:);

% Among essentially equivalent convergence rates,
% prioritize regulation and then control effort.
eligible_table = sortrows( ...
    eligible_table, ...
    {'MeanRMS', ...
     'MeanEnergy', ...
     'SaturationIncidence'}, ...
    {'ascend', ...
     'ascend', ...
     'ascend'});

best_LQG = eligible_table(1,:);

%% ==========================================================
% Recover Best Q, R, W, V, K, and L
% ==========================================================

best_K_index = best_LQG.KIndex;

Q_best = diag([ ...
    best_LQG.Q_q, ...
    best_LQG.Q_qdot, ...
    best_LQG.Q_Phi, ...
    best_LQG.Q_T, ...
    best_LQG.Q_Ve]);

R_best = best_LQG.R;

W_best = ...
    best_LQG.W_Scale * W_nominal;

V_best = ...
    best_LQG.V_Scale * V_nominal;

K_best = lqr( ...
    A, ...
    B, ...
    Q_best, ...
    R_best);

L_best = lqe( ...
    A, ...
    eye(5), ...
    C, ...
    W_best, ...
    V_best);

%% ==========================================================
% Final Results
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('BEST LQG CONFIGURATION\n');
fprintf('=========================================================\n');

disp(best_LQG);

fprintf('\nBest Q matrix:\n');
disp(Q_best);

fprintf('Best R:\n');
disp(R_best);

fprintf('Best W matrix:\n');
disp(W_best);

fprintf('Best V matrix:\n');
disp(V_best);

fprintf('Best LQR gain K:\n');
disp(K_best);

fprintf('Best Kalman gain L:\n');
disp(L_best);

fprintf('Closed-loop controller eigenvalues:\n');
disp(eig(A-B*K_best));

fprintf('Estimator eigenvalues:\n');
disp(eig(A-L_best*C));

%% ==========================================================
% Final Summary
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('FINAL LQG MONTE CARLO SUMMARY\n');
fprintf('=========================================================\n');

fprintf('Fixed initial conditions: %d\n',N_IC);
fprintf('Simulation duration: %.2f s\n',Tsim);

fprintf('LQR candidates evaluated: %d\n', ...
    N_K_candidates);

fprintf('Top LQR controllers retained: %d\n', ...
    N_top_K);

fprintf('LQE candidates per controller: %d\n', ...
    N_L_candidates);

fprintf('Total Stage 1 simulations: %d\n', ...
    N_K_candidates*N_IC);

fprintf('Total Stage 2 simulations: %d\n', ...
    N_top_K*N_L_candidates*N_IC);

fprintf('Total nonlinear simulations: %d\n', ...
    N_K_candidates*N_IC + ...
    N_top_K*N_L_candidates*N_IC);

fprintf('\nBest convergence rate: %.2f %%\n', ...
    best_LQG.SuccessRate);

fprintf('Best mean RMS error: %.6f\n', ...
    best_LQG.MeanRMS);

fprintf('Best mean settling time: %.6f s\n', ...
    best_LQG.MeanSettlingTime);

fprintf('Best mean control energy: %.6f\n', ...
    best_LQG.MeanEnergy);

fprintf('Best mean maximum voltage: %.6f\n', ...
    best_LQG.MeanMaxVoltage);

fprintf('Best saturation incidence: %.2f %%\n', ...
    best_LQG.SaturationIncidence);

fprintf('Best mean saturation duration: %.6f %%\n', ...
    best_LQG.MeanSaturationDuration);

%% ==========================================================
% Local Nonlinear LQG Dynamics Function
% ==========================================================

function dz = lqg_nonlinear_dynamics( ...
    t,z,A,B,C,L,K,...
    M,C_v,K_v,...
    phi_p_0,T_0,V_e_0,...
    alpha_B,alpha_E,...
    gamma_Phi,k_Phi,k_E,...
    C_th,h,tau_e,V_max)

%% ----------------------------------------------------------
% Plant and observer states
% ----------------------------------------------------------

x = z(1:5);
x_hat = z(6:10);

%% ----------------------------------------------------------
% Estimated-state feedback
% ----------------------------------------------------------

u_cmd = -K*x_hat;

%% ----------------------------------------------------------
% Physical actuator saturation
% ----------------------------------------------------------

u = max(min(u_cmd,V_max),-V_max);

%% ----------------------------------------------------------
% Nonlinear plant states
% ----------------------------------------------------------

q     = x(1);
q_dot = x(2);
Phi_p = x(3);
T     = x(4);
V_e   = x(5);

%% ----------------------------------------------------------
% Nonlinear plant dynamics
% ----------------------------------------------------------

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
% Measurements
% ----------------------------------------------------------
%
% The observer is formulated in perturbation coordinates:
%
%   x = [q; qdot; delta_Phi_p; delta_T; delta_V_e]
%
% Therefore the measured output must also be expressed in
% perturbation coordinates.
%

y = C*x_pert;

%% ----------------------------------------------------------
% Estimated measurements
% ----------------------------------------------------------

y_hat = C*x_hat;

%% ----------------------------------------------------------
% Linear Kalman observer
% ----------------------------------------------------------

x_hat_dot = ...
    A*x_hat ...
    + B*u ...
    + L*(y-y_hat);

%% ----------------------------------------------------------
% Combined plant + observer dynamics
% ----------------------------------------------------------

dz = [
    x_dot;
    x_hat_dot
];

end