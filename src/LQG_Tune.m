%% ==========================================================
% LQG Monte Carlo Q/R/W/V Tuning
%% ==========================================================
%
% Stage 1:
%   Randomly sample independent LQR Q/R candidates using
%   equal-probability stratification across specified
%   scale intervals.
%
% Stage 2:
%   Retain the strongest LQR controllers and randomly sample
%   LQE W/V covariance scales using equal-probability
%   stratification across specified scale intervals.
%
% All candidates use the exact same Monte Carlo initial
% conditions for fair comparison.
%
% ==========================================================


%% ==========================================================
% Nominal Model Parameters
%% ==========================================================

M = 1.0;
C_v = 0.3;
K_v = 2.25;

phi_p_0 = 1.0;
T_0     = 0.0;
V_e_0   = 1.0;

alpha_B = 0.20;
alpha_E = 0.12;

gamma_Phi = 0.15;
k_Phi     = 0.04;
k_E       = 0.03;

C_th = 1.0;
h    = 0.25;

tau_e = 0.10;


%% ==========================================================
% Linearized State-Space Model
%% ==========================================================

A = [ ...
    0, 1, 0, 0, 0;
    -K_v/M, -C_v/M, ...
        2*alpha_B*phi_p_0/M, ...
        0, ...
        2*alpha_E*V_e_0/M;
    0, 0, -gamma_Phi, 0, 0;
    0, 0, ...
        2*k_Phi*phi_p_0/C_th, ...
        -h/C_th, ...
        2*k_E*V_e_0/C_th;
    0, 0, 0, 0, -1/tau_e];


B = [ ...
    0;
    0;
    0;
    0;
    1/tau_e];


C = [ ...
    1 0 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];


%% ==========================================================
% Verify Linearized Model Matrices
%% ==========================================================

if ~exist('A','var')
    error( ...
        ['Matrix A was not found. ' ...
         'Load or generate the reduced-order linear model ' ...
         'before running this script.']);
end

if ~exist('B','var')
    error( ...
        ['Matrix B was not found. ' ...
         'Load or generate the reduced-order linear model ' ...
         'before running this script.']);
end

if ~exist('C','var')
    error( ...
        ['Matrix C was not found. ' ...
         'Load or generate the measurement matrix ' ...
         'before running this script.']);
end


%% ==========================================================
% Basic Dimension Checks
%% ==========================================================

if ~isequal(size(A),[5 5])
    error('A must be a 5 x 5 matrix.');
end

if ~isequal(size(B),[5 1])
    error('B must be a 5 x 1 matrix.');
end

if size(C,2) ~= 5
    error('C must have 5 columns.');
end


%% ==========================================================
% Reproducibility
%% ==========================================================
%
% A fixed random seed makes the Monte Carlo tuning study
% reproducible while retaining random sampling.
%
% Change this value to obtain a different realization.
%
% ==========================================================

rng(12345,'twister');


%% ==========================================================
% Simulation Settings
%% ==========================================================

Tsim = 10;

tspan = 0:0.01:Tsim;


%% ----------------------------------------------------------
% Convergence Criterion
%% ----------------------------------------------------------

eps_q     = 1e-2;
eps_q_dot = 1e-2;


%% ----------------------------------------------------------
% Settling Criterion
%% ----------------------------------------------------------

settle_tol = 2e-2;


%% ----------------------------------------------------------
% Physical Actuator-Voltage Saturation Limit
%% ----------------------------------------------------------

V_max = 10;


%% ==========================================================
% Nominal LQE Covariance Matrices
%% ==========================================================
%
% W and V are subsequently scaled by randomly sampled
% dimensionless factors.
%
% W_scale = 1 corresponds to I_5.
%
% V_scale = 1 corresponds to I_4.
%
% ==========================================================

W_nominal = eye(5);

V_nominal = eye(size(C,1));


%% ==========================================================
% Fixed Monte Carlo Initial-Condition Ensemble
%% ==========================================================
%
% Load the exact 1000 initial conditions generated previously
% by MC_initial_conditions_1000.mat.
%
% Every LQR and LQG candidate uses this identical ensemble.
%
% ==========================================================

IC_file = 'MC_initial_conditions_1000.mat';

if ~isfile(IC_file)

    error( ...
        'Required Monte Carlo IC file not found: %s', ...
        IC_file);

end

S = load(IC_file);


%% ----------------------------------------------------------
% Identify Variables in IC File
%% ----------------------------------------------------------

vars = fieldnames(S);

fprintf('\n');
fprintf('Loaded Monte Carlo IC file: %s\n',IC_file);
fprintf('Variables found in file:\n');

for k = 1:numel(vars)

    value = S.(vars{k});

    if isnumeric(value)

        fprintf( ...
            '  %-30s [%s]\n', ...
            vars{k}, ...
            strjoin(string(size(value)),' x '));

    else

        fprintf( ...
            '  %-30s (%s)\n', ...
            vars{k}, ...
            class(value));

    end

end


%% ==========================================================
% Locate Initial-Condition Matrix
%% ==========================================================
%
% Expected format:
%
%   N x 5
%
% with columns:
%
%   [q, qdot, deltaPhi_p, deltaT, deltaV_e]
%
% ==========================================================

ICs = [];

for k = 1:numel(vars)

    value = S.(vars{k});

    if isnumeric(value) && ...
            ismatrix(value) && ...
            size(value,2) == 5

        ICs = value;
        IC_variable_name = vars{k};
        break;

    end

end

if isempty(ICs)

    error([ ...
        'Could not find a numeric N x 5 initial-condition ' ...
        'matrix in %s. Expected columns: ' ...
        '[q, qdot, deltaPhi_p, deltaT, deltaV_e].'], ...
        IC_file);

end


%% ----------------------------------------------------------
% Verify Number of Initial Conditions
%% ----------------------------------------------------------

N_IC = size(ICs,1);

if N_IC ~= 1000

    warning( ...
        'Expected 1000 initial conditions, but found %d.', ...
        N_IC);

end


%% ----------------------------------------------------------
% Force IC Matrix to Double Precision
%% ----------------------------------------------------------

ICs = double(ICs);

fprintf('\n');
fprintf('Monte Carlo IC variable: %s\n', ...
    IC_variable_name);

fprintf('Number of ICs loaded: %d\n', ...
    N_IC);

fprintf('IC matrix dimensions: %d x %d\n', ...
    size(ICs,1),size(ICs,2));

fprintf('\nInitial-condition ranges:\n');

fprintf('q          : [% .6f, % .6f]\n', ...
    min(ICs(:,1)),max(ICs(:,1)));

fprintf('qdot       : [% .6f, % .6f]\n', ...
    min(ICs(:,2)),max(ICs(:,2)));

fprintf('deltaPhi_p  : [% .6f, % .6f]\n', ...
    min(ICs(:,3)),max(ICs(:,3)));

fprintf('deltaT      : [% .6f, % .6f]\n', ...
    min(ICs(:,4)),max(ICs(:,4)));

fprintf('deltaV_e    : [% .6f, % .6f]\n', ...
    min(ICs(:,5)),max(ICs(:,5)));

fprintf('\n');
fprintf( ...
    'All LQR and LQG candidates will use these exact ICs.\n');


%% ==========================================================
% STAGE 1: Stratified-by-Scale Random LQR Q/R Sweep
%% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('STAGE 1: STRATIFIED-BY-SCALE RANDOM LQR Q/R SWEEP\n');
fprintf('=========================================================\n');


%% ----------------------------------------------------------
% Number of Random Q/R Candidates
%% ----------------------------------------------------------

N_QR_candidates = 1000;


%% ----------------------------------------------------------
% Number of Top LQR Controllers Retained
%% ----------------------------------------------------------

N_top_K_requested = 20;


%% ==========================================================
% Stratified-by-Scale Sampling Intervals
%% ==========================================================
%
% For each Q/R parameter:
%
%   1. One of the three intervals is selected with equal
%      probability.
%
%   2. A value is sampled uniformly within that interval.
%
% Each parameter is sampled independently.
%
% No Cartesian product is formed.
%
% ==========================================================

Qq_intervals = [ ...
    1    10;
    10   100;
    100  500];


Qqd_intervals = [ ...
    0.1    1;
    1      10;
    10     100];


Qphi_intervals = [ ...
    0.01   0.1;
    0.1    1;
    1      10];


QT_intervals = [ ...
    0.01   0.1;
    0.1    1;
    1      10];


QVe_intervals = [ ...
    0.01   0.1;
    0.1    1;
    1      10];


R_intervals = [ ...
    0.01   0.1;
    0.1    1;
    1      10];


%% ==========================================================
% Generate Stratified Q/R Candidates
%% ==========================================================

[Qq_random,Qq_interval_index] = ...
    sample_equal_probability_intervals( ...
    Qq_intervals, ...
    N_QR_candidates);


[Qqd_random,Qqd_interval_index] = ...
    sample_equal_probability_intervals( ...
    Qqd_intervals, ...
    N_QR_candidates);


[Qphi_random,Qphi_interval_index] = ...
    sample_equal_probability_intervals( ...
    Qphi_intervals, ...
    N_QR_candidates);


[QT_random,QT_interval_index] = ...
    sample_equal_probability_intervals( ...
    QT_intervals, ...
    N_QR_candidates);


[QVe_random,QVe_interval_index] = ...
    sample_equal_probability_intervals( ...
    QVe_intervals, ...
    N_QR_candidates);


[R_random,R_interval_index] = ...
    sample_equal_probability_intervals( ...
    R_intervals, ...
    N_QR_candidates);


%% ==========================================================
% Q/R Sampling-Balance Check
%% ==========================================================

expected_count = ...
    N_QR_candidates / 3;

fprintf('\n');
fprintf('=========================================================\n');
fprintf('Q/R SAMPLING-BALANCE CHECK\n');
fprintf('=========================================================\n');

fprintf('\nExpected samples per interval: %.2f\n', ...
    expected_count);


%% ----------------------------------------------------------
% Q_q Balance
%% ----------------------------------------------------------

Qq_counts = ...
    accumarray( ...
    Qq_interval_index, ...
    1, ...
    [size(Qq_intervals,1),1]);

fprintf('\nQ_q interval counts:\n');

for interval_index = 1:size(Qq_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        Qq_intervals(interval_index,1), ...
        Qq_intervals(interval_index,2), ...
        Qq_counts(interval_index), ...
        100*Qq_counts(interval_index)/N_QR_candidates);

end


%% ----------------------------------------------------------
% Q_qdot Balance
%% ----------------------------------------------------------

Qqd_counts = ...
    accumarray( ...
    Qqd_interval_index, ...
    1, ...
    [size(Qqd_intervals,1),1]);

fprintf('\nQ_qdot interval counts:\n');

for interval_index = 1:size(Qqd_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        Qqd_intervals(interval_index,1), ...
        Qqd_intervals(interval_index,2), ...
        Qqd_counts(interval_index), ...
        100*Qqd_counts(interval_index)/N_QR_candidates);

end


%% ----------------------------------------------------------
% Q_Phi Balance
%% ----------------------------------------------------------

Qphi_counts = ...
    accumarray( ...
    Qphi_interval_index, ...
    1, ...
    [size(Qphi_intervals,1),1]);

fprintf('\nQ_Phi interval counts:\n');

for interval_index = 1:size(Qphi_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        Qphi_intervals(interval_index,1), ...
        Qphi_intervals(interval_index,2), ...
        Qphi_counts(interval_index), ...
        100*Qphi_counts(interval_index)/N_QR_candidates);

end


%% ----------------------------------------------------------
% Q_T Balance
%% ----------------------------------------------------------

QT_counts = ...
    accumarray( ...
    QT_interval_index, ...
    1, ...
    [size(QT_intervals,1),1]);

fprintf('\nQ_T interval counts:\n');

for interval_index = 1:size(QT_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        QT_intervals(interval_index,1), ...
        QT_intervals(interval_index,2), ...
        QT_counts(interval_index), ...
        100*QT_counts(interval_index)/N_QR_candidates);

end


%% ----------------------------------------------------------
% Q_Ve Balance
%% ----------------------------------------------------------

QVe_counts = ...
    accumarray( ...
    QVe_interval_index, ...
    1, ...
    [size(QVe_intervals,1),1]);

fprintf('\nQ_Ve interval counts:\n');

for interval_index = 1:size(QVe_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        QVe_intervals(interval_index,1), ...
        QVe_intervals(interval_index,2), ...
        QVe_counts(interval_index), ...
        100*QVe_counts(interval_index)/N_QR_candidates);

end


%% ----------------------------------------------------------
% R Balance
%% ----------------------------------------------------------

R_counts = ...
    accumarray( ...
    R_interval_index, ...
    1, ...
    [size(R_intervals,1),1]);

fprintf('\nR interval counts:\n');

for interval_index = 1:size(R_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        R_intervals(interval_index,1), ...
        R_intervals(interval_index,2), ...
        R_counts(interval_index), ...
        100*R_counts(interval_index)/N_QR_candidates);

end


%% ==========================================================
% Q/R Candidate Table
%% ==========================================================

QR_candidate_table = table( ...
    (1:N_QR_candidates).', ...
    Qq_interval_index, ...
    Qq_random, ...
    Qqd_interval_index, ...
    Qqd_random, ...
    Qphi_interval_index, ...
    Qphi_random, ...
    QT_interval_index, ...
    QT_random, ...
    QVe_interval_index, ...
    QVe_random, ...
    R_interval_index, ...
    R_random, ...
    'VariableNames',{ ...
    'Candidate', ...
    'Q_q_Interval', ...
    'Q_q', ...
    'Q_qdot_Interval', ...
    'Q_qdot', ...
    'Q_Phi_Interval', ...
    'Q_Phi', ...
    'Q_T_Interval', ...
    'Q_T', ...
    'Q_Ve_Interval', ...
    'Q_Ve', ...
    'R_Interval', ...
    'R'});


fprintf('\n');
fprintf('=========================================================\n');
fprintf('STRATIFIED Q/R CANDIDATES\n');
fprintf('=========================================================\n');

disp(QR_candidate_table);

fprintf( ...
    'LQR Q/R candidates: %d\n', ...
    N_QR_candidates);

fprintf( ...
    'Initial conditions per candidate: %d\n', ...
    N_IC);

fprintf( ...
    'Stage 1 nonlinear simulations: %d\n', ...
    N_QR_candidates*N_IC);


%% ==========================================================
% Stage 1 Result Storage
%% ==========================================================

lqr_success     = zeros(N_QR_candidates,1);
lqr_RMS         = zeros(N_QR_candidates,1);
lqr_energy      = zeros(N_QR_candidates,1);
lqr_settle      = zeros(N_QR_candidates,1);
lqr_maxVoltage  = zeros(N_QR_candidates,1);
lqr_saturation  = zeros(N_QR_candidates,1);
lqr_satDuration = zeros(N_QR_candidates,1);

lqr_K = zeros(N_QR_candidates,5);


%% ==========================================================
% STAGE 1 LQR Q/R SWEEP
%% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('BEGINNING STAGE 1 LQR SWEEP\n');
fprintf('=========================================================\n');

tic;


for candidate = 1:N_QR_candidates

    %% ------------------------------------------------------
    % Construct Candidate Q and R
    %% ------------------------------------------------------

    Q_test = diag([ ...
        Qq_random(candidate), ...
        Qqd_random(candidate), ...
        Qphi_random(candidate), ...
        QT_random(candidate), ...
        QVe_random(candidate)]);

    R_test = R_random(candidate);


    %% ------------------------------------------------------
    % Compute LQR Gain
    %% ------------------------------------------------------

    try

        K_test = lqr( ...
            A, ...
            B, ...
            Q_test, ...
            R_test);

    catch ME

        warning( ...
            'LQR failed for candidate %d: %s', ...
            candidate, ...
            ME.message);

        lqr_success(candidate)     = NaN;
        lqr_RMS(candidate)         = NaN;
        lqr_energy(candidate)      = NaN;
        lqr_settle(candidate)      = NaN;
        lqr_maxVoltage(candidate)  = NaN;
        lqr_saturation(candidate)  = NaN;
        lqr_satDuration(candidate) = NaN;

        continue;

    end


    %% ------------------------------------------------------
    % Store LQR Gain
    %% ------------------------------------------------------

    lqr_K(candidate,:) = K_test;


    %% ======================================================
    % Monte Carlo Metric Storage
    %% ======================================================

    RMS_i = zeros(N_IC,1);

    energy_i = zeros(N_IC,1);

    settle_i = zeros(N_IC,1);

    maxVoltage_i = zeros(N_IC,1);

    saturated_i = false(N_IC,1);

    saturationDuration_i = ...
        zeros(N_IC,1);

    converged_i = ...
        false(N_IC,1);


    %% ======================================================
    % Monte Carlo Simulations
    %% ======================================================

    for i = 1:N_IC

        %% --------------------------------------------------
        % Fixed Initial Condition
        %% --------------------------------------------------

        dx_0 = ICs(i,:).';

        x_0 = [ ...
            dx_0(1);
            dx_0(2);
            phi_p_0 + dx_0(3);
            T_0     + dx_0(4);
            V_e_0   + dx_0(5)];


        %% --------------------------------------------------
        % Nonlinear LQR Dynamics
        %% --------------------------------------------------

        f_lqr = @(t,x) ...
            lqr_nonlinear_dynamics( ...
            t,x,K_test, ...
            M,C_v,K_v, ...
            phi_p_0,T_0,V_e_0, ...
            alpha_B,alpha_E, ...
            gamma_Phi,k_Phi,k_E, ...
            C_th,h,tau_e,V_max);


        [t,x] = ode45( ...
            f_lqr, ...
            tspan, ...
            x_0);


        %% --------------------------------------------------
        % Extract States
        %% --------------------------------------------------

        q = x(:,1);

        q_dot = x(:,2);


        %% --------------------------------------------------
        % Reconstruct Physical Control
        %% --------------------------------------------------

        x_pert = x;

        x_pert(:,3) = ...
            x(:,3) - phi_p_0;

        x_pert(:,4) = ...
            x(:,4) - T_0;

        x_pert(:,5) = ...
            x(:,5) - V_e_0;

        u_cmd = ...
            -x_pert*K_test.';

        u = max( ...
            min(u_cmd,V_max), ...
            -V_max);


        %% --------------------------------------------------
        % Performance Metrics
        %% --------------------------------------------------

        RMS_i(i) = ...
            sqrt(mean(q.^2));

        energy_i(i) = ...
            trapz(t,u.^2);

        maxVoltage_i(i) = ...
            max(abs(u));


        %% --------------------------------------------------
        % Saturation
        %% --------------------------------------------------

        saturated = ...
            abs(u) >= V_max - 1e-10;

        saturated_i(i) = ...
            any(saturated);

        saturationDuration_i(i) = ...
            100* ...
            trapz( ...
                t, ...
                double(saturated)) / Tsim;


        %% --------------------------------------------------
        % Convergence
        %% --------------------------------------------------

        converged_i(i) = ...
            abs(q(end)) < eps_q && ...
            abs(q_dot(end)) < eps_q_dot;


        %% --------------------------------------------------
        % Settling Time
        %% --------------------------------------------------

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


    %% ======================================================
    % Candidate Summary
    %% ======================================================

    lqr_success(candidate) = ...
        100*mean(converged_i);

    lqr_RMS(candidate) = ...
        mean(RMS_i);

    lqr_energy(candidate) = ...
        mean(energy_i);

    lqr_settle(candidate) = ...
        mean(settle_i);

    lqr_maxVoltage(candidate) = ...
        mean(maxVoltage_i);

    lqr_saturation(candidate) = ...
        100*mean(saturated_i);

    lqr_satDuration(candidate) = ...
        mean(saturationDuration_i);


    %% ------------------------------------------------------
    % Progress
    %% ------------------------------------------------------

    fprintf( ...
        'LQR candidate %4d/%4d | ', ...
        candidate, ...
        N_QR_candidates);

    fprintf( ...
        'Success = %.1f %% | ', ...
        lqr_success(candidate));

    fprintf( ...
        'RMS = %.4g | ', ...
        lqr_RMS(candidate));

    fprintf( ...
        'Energy = %.4g\n', ...
        lqr_energy(candidate));

end


elapsed_stage1 = toc;


fprintf('\n');
fprintf( ...
    'Stage 1 completed in %.2f minutes.\n', ...
    elapsed_stage1/60);


%% ==========================================================
% Stage 1 LQR Results Table
%% ==========================================================

K_results_table = table( ...
    (1:N_QR_candidates).', ...
    Qq_random, ...
    Qqd_random, ...
    Qphi_random, ...
    QT_random, ...
    QVe_random, ...
    R_random, ...
    Qq_interval_index, ...
    Qqd_interval_index, ...
    Qphi_interval_index, ...
    QT_interval_index, ...
    QVe_interval_index, ...
    R_interval_index, ...
    lqr_K(:,1), ...
    lqr_K(:,2), ...
    lqr_K(:,3), ...
    lqr_K(:,4), ...
    lqr_K(:,5), ...
    lqr_success, ...
    lqr_RMS, ...
    lqr_settle, ...
    lqr_energy, ...
    lqr_maxVoltage, ...
    lqr_saturation, ...
    lqr_satDuration, ...
    'VariableNames',{ ...
    'Candidate', ...
    'Q_q', ...
    'Q_qdot', ...
    'Q_Phi', ...
    'Q_T', ...
    'Q_Ve', ...
    'R', ...
    'Q_q_Interval', ...
    'Q_qdot_Interval', ...
    'Q_Phi_Interval', ...
    'Q_T_Interval', ...
    'Q_Ve_Interval', ...
    'R_Interval', ...
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


%% ==========================================================
% Remove Failed LQR Candidates
%% ==========================================================

valid_LQR = ...
    isfinite(K_results_table.SuccessRate);

K_results_table = ...
    K_results_table(valid_LQR,:);


%% ==========================================================
% Sort LQR Results
%% ==========================================================

K_results_table = sortrows( ...
    K_results_table, ...
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
fprintf('TOP LQR CONTROLLERS\n');
fprintf('=========================================================\n');

disp(K_results_table( ...
    1:min(20,height(K_results_table)),:));


%% ==========================================================
% Retain Top LQR Controllers for Stage 2
%% ==========================================================

N_top_K = ...
    min(N_top_K_requested,height(K_results_table));

top_K_table = ...
    K_results_table(1:N_top_K,:);


fprintf('\n');
fprintf('=========================================================\n');
fprintf('LQR CONTROLLERS RETAINED FOR STAGE 2\n');
fprintf('=========================================================\n');

fprintf( ...
    'Top LQR controllers retained: %d\n', ...
    N_top_K);

disp(top_K_table);


%% ==========================================================
% STAGE 2: Stratified-by-Scale Random LQE W/V Sweep
%% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('STAGE 2: STRATIFIED-BY-SCALE RANDOM LQE W/V SWEEP\n');
fprintf('=========================================================\n');


%% ----------------------------------------------------------
% Number of Random W/V Candidates
%% ----------------------------------------------------------

N_WV_candidates = 1000;


%% ==========================================================
% Equal-Probability W/V Scale Intervals
%% ==========================================================
%
% Every W candidate independently selects one interval with
% probability 1/3 and samples uniformly within that interval.
%
% The same is done independently for V.
%
% ==========================================================

W_scale_intervals = [ ...
    0.01  0.10;
    0.10  1.00;
    1.00 10.00];


V_scale_intervals = [ ...
    0.01  0.10;
    0.10  1.00;
    1.00 10.00];


%% ==========================================================
% Generate Common Random W/V Candidates
%% ==========================================================

[W_scale_random,W_interval_index] = ...
    sample_equal_probability_intervals( ...
    W_scale_intervals, ...
    N_WV_candidates);


[V_scale_random,V_interval_index] = ...
    sample_equal_probability_intervals( ...
    V_scale_intervals, ...
    N_WV_candidates);


%% ==========================================================
% W/V Sampling-Balance Check
%% ==========================================================

expected_WV_count = ...
    N_WV_candidates/3;


fprintf('\n');
fprintf('=========================================================\n');
fprintf('W/V SAMPLING-BALANCE CHECK\n');
fprintf('=========================================================\n');

fprintf( ...
    '\nExpected samples per interval: %.2f\n', ...
    expected_WV_count);


%% ----------------------------------------------------------
% W Balance
%% ----------------------------------------------------------

W_counts = ...
    accumarray( ...
    W_interval_index, ...
    1, ...
    [size(W_scale_intervals,1),1]);


fprintf('\nW interval counts:\n');

for interval_index = 1:size(W_scale_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        W_scale_intervals(interval_index,1), ...
        W_scale_intervals(interval_index,2), ...
        W_counts(interval_index), ...
        100*W_counts(interval_index)/N_WV_candidates);

end


%% ----------------------------------------------------------
% V Balance
%% ----------------------------------------------------------

V_counts = ...
    accumarray( ...
    V_interval_index, ...
    1, ...
    [size(V_scale_intervals,1),1]);


fprintf('\nV interval counts:\n');

for interval_index = 1:size(V_scale_intervals,1)

    fprintf( ...
        '  Interval %d [%g, %g]: %4d samples (%.1f %%)\n', ...
        interval_index, ...
        V_scale_intervals(interval_index,1), ...
        V_scale_intervals(interval_index,2), ...
        V_counts(interval_index), ...
        100*V_counts(interval_index)/N_WV_candidates);

end


%% ==========================================================
% W/V Candidate Table
%% ==========================================================

WV_candidate_table = table( ...
    (1:N_WV_candidates).', ...
    W_interval_index, ...
    W_scale_random, ...
    V_interval_index, ...
    V_scale_random, ...
    'VariableNames',{ ...
    'Candidate', ...
    'W_Interval', ...
    'W_Scale', ...
    'V_Interval', ...
    'V_Scale'});


fprintf('\n');
fprintf('=========================================================\n');
fprintf('STRATIFIED W/V OBSERVER CANDIDATES\n');
fprintf('=========================================================\n');

disp(WV_candidate_table);

fprintf( ...
    'LQE W/V candidates per controller: %d\n', ...
    N_WV_candidates);

fprintf( ...
    'Top LQR controllers evaluated: %d\n', ...
    N_top_K);

fprintf( ...
    'Initial conditions per LQG candidate: %d\n', ...
    N_IC);

fprintf( ...
    'Stage 2 nonlinear simulations: %d\n', ...
    N_top_K*N_WV_candidates*N_IC);


%% ==========================================================
% Stage 2 Result Storage
%% ==========================================================

N_LQG_candidates = ...
    N_top_K*N_WV_candidates;


lqg_K_index = ...
    zeros(N_LQG_candidates,1);


lqg_W_scale = ...
    zeros(N_LQG_candidates,1);


lqg_V_scale = ...
    zeros(N_LQG_candidates,1);


lqg_W_interval = ...
    zeros(N_LQG_candidates,1);


lqg_V_interval = ...
    zeros(N_LQG_candidates,1);


lqg_success = ...
    zeros(N_LQG_candidates,1);


lqg_RMS = ...
    zeros(N_LQG_candidates,1);


lqg_energy = ...
    zeros(N_LQG_candidates,1);


lqg_settle = ...
    zeros(N_LQG_candidates,1);


lqg_maxVoltage = ...
    zeros(N_LQG_candidates,1);


lqg_saturation = ...
    zeros(N_LQG_candidates,1);


lqg_satDuration = ...
    zeros(N_LQG_candidates,1);


lqg_L = ...
    zeros(5,4,N_LQG_candidates);


%% ==========================================================
% STAGE 2 W/V SWEEP
%% ==========================================================
%
% Each retained LQR controller is evaluated using the exact
% same set of randomly generated W/V observer candidates.
%
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('BEGINNING STAGE 2 LQG SWEEP\n');
fprintf('=========================================================\n');

tic;


candidate = 0;


for k_index = 1:N_top_K


    %% ------------------------------------------------------
    % Recover LQR Controller
    %% ------------------------------------------------------

    Q_test = diag([ ...
        top_K_table.Q_q(k_index), ...
        top_K_table.Q_qdot(k_index), ...
        top_K_table.Q_Phi(k_index), ...
        top_K_table.Q_T(k_index), ...
        top_K_table.Q_Ve(k_index)]);


    R_test = ...
        top_K_table.R(k_index);


    K_test = lqr( ...
        A, ...
        B, ...
        Q_test, ...
        R_test);


    %% ------------------------------------------------------
    % Evaluate Every Common W/V Candidate
    %% ------------------------------------------------------

    for wv_index = 1:N_WV_candidates

        candidate = candidate + 1;


        W_scale = ...
            W_scale_random(wv_index);


        V_scale = ...
            V_scale_random(wv_index);


        %% --------------------------------------------------
        % Construct W and V
        %% --------------------------------------------------

        W_test = ...
            W_scale*W_nominal;


        V_test = ...
            V_scale*V_nominal;


        %% --------------------------------------------------
        % Compute Kalman Gain
        %% --------------------------------------------------

        try

            L_test = lqe( ...
                A, ...
                eye(5), ...
                C, ...
                W_test, ...
                V_test);

        catch ME

            warning( ...
                'LQE failed for candidate %d: %s', ...
                candidate, ...
                ME.message);


            lqg_K_index(candidate) = ...
                k_index;


            lqg_W_scale(candidate) = ...
                W_scale;


            lqg_V_scale(candidate) = ...
                V_scale;


            lqg_W_interval(candidate) = ...
                W_interval_index(wv_index);


            lqg_V_interval(candidate) = ...
                V_interval_index(wv_index);


            lqg_success(candidate) = NaN;

            lqg_RMS(candidate) = NaN;

            lqg_energy(candidate) = NaN;

            lqg_settle(candidate) = NaN;

            lqg_maxVoltage(candidate) = NaN;

            lqg_saturation(candidate) = NaN;

            lqg_satDuration(candidate) = NaN;


            continue;

        end


        %% ==================================================
        % Monte Carlo Metric Storage
        %% ==================================================

        RMS_i = zeros(N_IC,1);

        energy_i = zeros(N_IC,1);

        settle_i = zeros(N_IC,1);

        maxVoltage_i = zeros(N_IC,1);

        saturated_i = false(N_IC,1);

        saturationDuration_i = ...
            zeros(N_IC,1);

        converged_i = ...
            false(N_IC,1);


        %% ==================================================
        % Monte Carlo Simulations
        %% ==================================================

        for i = 1:N_IC


            %% --------------------------------------------------
            % Fixed Initial Condition
            %% --------------------------------------------------

            dx_0 = ICs(i,:).';


            x_0 = [ ...
                dx_0(1);
                dx_0(2);
                phi_p_0 + dx_0(3);
                T_0     + dx_0(4);
                V_e_0   + dx_0(5)];


            %% --------------------------------------------------
            % Observer Initial Condition
            %% --------------------------------------------------

            xhat_0 = zeros(5,1);


            z_0 = [ ...
                x_0;
                xhat_0];


            %% --------------------------------------------------
            % Nonlinear LQG Dynamics
            %% --------------------------------------------------

            f_lqg = @(t,z) ...
                lqg_nonlinear_dynamics( ...
                t,z,A,B,C,L_test,K_test, ...
                M,C_v,K_v, ...
                phi_p_0,T_0,V_e_0, ...
                alpha_B,alpha_E, ...
                gamma_Phi,k_Phi,k_E, ...
                C_th,h,tau_e,V_max);


            [t,z] = ode45( ...
                f_lqg, ...
                tspan, ...
                z_0);


            %% --------------------------------------------------
            % Extract Plant States
            %% --------------------------------------------------

            x = z(:,1:5);


            q = x(:,1);

            q_dot = x(:,2);


            %% --------------------------------------------------
            % Reconstruct Physical Control
            %% --------------------------------------------------

            x_hat = ...
                z(:,6:10);


            u_cmd = ...
                -x_hat*K_test.';


            u = max( ...
                min(u_cmd,V_max), ...
                -V_max);


            %% --------------------------------------------------
            % Performance Metrics
            %% --------------------------------------------------

            RMS_i(i) = ...
                sqrt(mean(q.^2));


            energy_i(i) = ...
                trapz(t,u.^2);


            maxVoltage_i(i) = ...
                max(abs(u));


            %% --------------------------------------------------
            % Saturation
            %% --------------------------------------------------

            saturated = ...
                abs(u) >= V_max - 1e-10;


            saturated_i(i) = ...
                any(saturated);


            saturationDuration_i(i) = ...
                100* ...
                trapz( ...
                t, ...
                double(saturated)) / Tsim;


            %% --------------------------------------------------
            % Convergence
            %% --------------------------------------------------

            converged_i(i) = ...
                abs(q(end)) < eps_q && ...
                abs(q_dot(end)) < eps_q_dot;


            %% --------------------------------------------------
            % Settling Time
            %% --------------------------------------------------

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
        % Candidate Summary
        %% ==================================================

        lqg_K_index(candidate) = ...
            k_index;


        lqg_W_scale(candidate) = ...
            W_scale;


        lqg_V_scale(candidate) = ...
            V_scale;


        lqg_W_interval(candidate) = ...
            W_interval_index(wv_index);


        lqg_V_interval(candidate) = ...
            V_interval_index(wv_index);


        lqg_L(:,:,candidate) = ...
            L_test;


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
        %% --------------------------------------------------

        fprintf( ...
            'LQG candidate %4d/%4d | K #%2d | ', ...
            candidate, ...
            N_LQG_candidates, ...
            k_index);


        fprintf( ...
            'W = %.4g | V = %.4g | ', ...
            W_scale, ...
            V_scale);


        fprintf( ...
            'Success = %.1f %% | RMS = %.4g\n', ...
            lqg_success(candidate), ...
            lqg_RMS(candidate));

    end

end


elapsed_stage2 = toc;


fprintf('\n');
fprintf( ...
    'Stage 2 completed in %.2f minutes.\n', ...
    elapsed_stage2/60);


%% ==========================================================
% Stage 2 LQG Results Table
%% ==========================================================

Qq_result = ...
    zeros(N_LQG_candidates,1);

Qqd_result = ...
    zeros(N_LQG_candidates,1);

Qphi_result = ...
    zeros(N_LQG_candidates,1);

QT_result = ...
    zeros(N_LQG_candidates,1);

QVe_result = ...
    zeros(N_LQG_candidates,1);

R_result = ...
    zeros(N_LQG_candidates,1);


for i = 1:N_LQG_candidates

    k_index_i = ...
        lqg_K_index(i);


    Qq_result(i) = ...
        top_K_table.Q_q(k_index_i);


    Qqd_result(i) = ...
        top_K_table.Q_qdot(k_index_i);


    Qphi_result(i) = ...
        top_K_table.Q_Phi(k_index_i);


    QT_result(i) = ...
        top_K_table.Q_T(k_index_i);


    QVe_result(i) = ...
        top_K_table.Q_Ve(k_index_i);


    R_result(i) = ...
        top_K_table.R(k_index_i);

end


%% ----------------------------------------------------------
% Create LQG Results Table
%% ----------------------------------------------------------

LQG_results_table = table( ...
    lqg_K_index, ...
    Qq_result, ...
    Qqd_result, ...
    Qphi_result, ...
    QT_result, ...
    QVe_result, ...
    R_result, ...
    lqg_W_interval, ...
    lqg_W_scale, ...
    lqg_V_interval, ...
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
    'W_Interval', ...
    'W_Scale', ...
    'V_Interval', ...
    'V_Scale', ...
    'SuccessRate', ...
    'MeanRMS', ...
    'MeanSettlingTime', ...
    'MeanEnergy', ...
    'MeanMaxVoltage', ...
    'SaturationIncidence', ...
    'MeanSaturationDuration'});


%% ==========================================================
% Remove Failed LQG Candidates
%% ==========================================================

valid_LQG = ...
    isfinite(LQG_results_table.SuccessRate);


LQG_results_table = ...
    LQG_results_table(valid_LQG,:);


%% ==========================================================
% Sort LQG Results
%% ==========================================================

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
%% ==========================================================

maxSuccess = ...
    max(LQG_results_table.SuccessRate);


success_margin = 1.0;


eligible = ...
    LQG_results_table.SuccessRate >= ...
    maxSuccess - success_margin;


eligible_table = ...
    LQG_results_table(eligible,:);


eligible_table = sortrows( ...
    eligible_table, ...
    {'MeanRMS', ...
     'MeanEnergy', ...
     'SaturationIncidence'}, ...
    {'ascend', ...
     'ascend', ...
     'ascend'});


best_LQG = ...
    eligible_table(1,:);


fprintf('\n');
fprintf('=========================================================\n');
fprintf('SELECTED BEST LQG CONFIGURATION\n');
fprintf('=========================================================\n');


fprintf( ...
    'Success Rate      = %.2f %%\n', ...
    best_LQG.SuccessRate);


fprintf( ...
    'Mean RMS          = %.6f\n', ...
    best_LQG.MeanRMS);


fprintf( ...
    'Mean Energy       = %.6f\n', ...
    best_LQG.MeanEnergy);


fprintf( ...
    'Mean Max Voltage  = %.6f\n', ...
    best_LQG.MeanMaxVoltage);


fprintf( ...
    'Saturation Inc.   = %.6f %%\n', ...
    best_LQG.SaturationIncidence);


fprintf( ...
    'Q = diag([%.4g %.4g %.4g %.4g %.4g])\n', ...
    best_LQG.Q_q, ...
    best_LQG.Q_qdot, ...
    best_LQG.Q_Phi, ...
    best_LQG.Q_T, ...
    best_LQG.Q_Ve);


fprintf( ...
    'R                 = %.4g\n', ...
    best_LQG.R);


fprintf( ...
    'W Scale           = %.4g\n', ...
    best_LQG.W_Scale);


fprintf( ...
    'W Interval        = %d\n', ...
    best_LQG.W_Interval);


fprintf( ...
    'V Scale           = %.4g\n', ...
    best_LQG.V_Scale);


fprintf( ...
    'V Interval        = %d\n', ...
    best_LQG.V_Interval);


fprintf('=========================================================\n');


%% ==========================================================
% Recover Best Q, R, W, V, K, and L
%% ==========================================================

Q_best = diag([ ...
    best_LQG.Q_q, ...
    best_LQG.Q_qdot, ...
    best_LQG.Q_Phi, ...
    best_LQG.Q_T, ...
    best_LQG.Q_Ve]);


R_best = ...
    best_LQG.R;


W_best = ...
    best_LQG.W_Scale*W_nominal;


V_best = ...
    best_LQG.V_Scale*V_nominal;


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
% Final LQG Results
%% ==========================================================

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
disp(eig(A - B*K_best));


fprintf('Estimator eigenvalues:\n');
disp(eig(A - L_best*C));


%% ==========================================================
% Final Summary
%% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('FINAL LQG MONTE CARLO SUMMARY\n');
fprintf('=========================================================\n');


fprintf( ...
    'Fixed initial conditions: %d\n', ...
    N_IC);


fprintf( ...
    'Simulation duration: %.2f s\n', ...
    Tsim);


fprintf( ...
    'LQR candidates evaluated: %d\n', ...
    N_QR_candidates);


fprintf( ...
    'Top LQR controllers retained: %d\n', ...
    N_top_K);


fprintf( ...
    'Random W/V candidates per controller: %d\n', ...
    N_WV_candidates);


fprintf( ...
    'Total Stage 1 simulations: %d\n', ...
    N_QR_candidates*N_IC);


fprintf( ...
    'Total Stage 2 simulations: %d\n', ...
    N_top_K*N_WV_candidates*N_IC);


fprintf( ...
    'Total nonlinear simulations: %d\n', ...
    N_QR_candidates*N_IC + ...
    N_top_K*N_WV_candidates*N_IC);


fprintf( ...
    'Best convergence rate: %.2f %%\n', ...
    best_LQG.SuccessRate);


fprintf( ...
    'Best mean RMS error: %.6f\n', ...
    best_LQG.MeanRMS);


fprintf( ...
    'Best mean settling time: %.6f s\n', ...
    best_LQG.MeanSettlingTime);


fprintf( ...
    'Best mean control energy: %.6f\n', ...
    best_LQG.MeanEnergy);


fprintf( ...
    'Best mean maximum voltage: %.6f\n', ...
    best_LQG.MeanMaxVoltage);


fprintf( ...
    'Best saturation incidence: %.2f %%\n', ...
    best_LQG.SaturationIncidence);


fprintf( ...
    'Best mean saturation duration: %.6f %%\n', ...
    best_LQG.MeanSaturationDuration);


%% ==========================================================
% BEST LQR CONFIGURATION
%% ==========================================================

best_K_row = ...
    K_results_table(1,:);


fprintf('\n');
fprintf('=========================================================\n');
fprintf('BEST LQR CONFIGURATION\n');
fprintf('=========================================================\n');


fprintf( ...
    'Success Rate          = %.2f %%\n', ...
    best_K_row.SuccessRate);


fprintf( ...
    'Mean RMS Error        = %.6f\n', ...
    best_K_row.MeanRMS);


fprintf( ...
    'Mean Settling Time    = %.6f s\n', ...
    best_K_row.MeanSettlingTime);


fprintf( ...
    'Mean Control Energy   = %.6f\n', ...
    best_K_row.MeanEnergy);


fprintf( ...
    'Mean Maximum Voltage  = %.6f\n', ...
    best_K_row.MeanMaxVoltage);


fprintf( ...
    'Saturation Incidence  = %.2f %%\n', ...
    best_K_row.SaturationIncidence);


fprintf( ...
    'Mean Saturation Time  = %.6f %%\n', ...
    best_K_row.MeanSaturationDuration);


fprintf('\n');


fprintf( ...
    'Q = diag([%.6g %.6g %.6g %.6g %.6g])\n', ...
    best_K_row.Q_q, ...
    best_K_row.Q_qdot, ...
    best_K_row.Q_Phi, ...
    best_K_row.Q_T, ...
    best_K_row.Q_Ve);


fprintf( ...
    'R = %.6g\n', ...
    best_K_row.R);


fprintf('\n');
fprintf('LQR Gain K:\n');


K_best_stage1 = [ ...
    best_K_row.K_q, ...
    best_K_row.K_qdot, ...
    best_K_row.K_Phi, ...
    best_K_row.K_T, ...
    best_K_row.K_Ve];


disp(K_best_stage1);


%% ==========================================================
% BEST CONFIGURATION BY INDIVIDUAL METRIC
%% ==========================================================

[~,idx_success] = ...
    max(K_results_table.SuccessRate);


[~,idx_rms] = ...
    min(K_results_table.MeanRMS);


[~,idx_settle] = ...
    min(K_results_table.MeanSettlingTime);


[~,idx_energy] = ...
    min(K_results_table.MeanEnergy);


[~,idx_voltage] = ...
    min(K_results_table.MeanMaxVoltage);


[~,idx_sat] = ...
    min(K_results_table.SaturationIncidence);


metric_idx = [ ...
    idx_success;
    idx_rms;
    idx_settle;
    idx_energy;
    idx_voltage;
    idx_sat];


metric_name = { ...
    'Highest Success Rate';
    'Lowest RMS Error';
    'Fastest Settling';
    'Lowest Control Energy';
    'Lowest Maximum Voltage';
    'Lowest Saturation Incidence'};


best_by_metric = table( ...
    metric_name, ...
    K_results_table.SuccessRate(metric_idx), ...
    K_results_table.MeanRMS(metric_idx), ...
    K_results_table.MeanSettlingTime(metric_idx), ...
    K_results_table.MeanEnergy(metric_idx), ...
    K_results_table.MeanMaxVoltage(metric_idx), ...
    K_results_table.SaturationIncidence(metric_idx), ...
    'VariableNames',{ ...
    'Metric', ...
    'SuccessRate', ...
    'MeanRMS', ...
    'MeanSettlingTime', ...
    'MeanEnergy', ...
    'MeanMaxVoltage', ...
    'SaturationIncidence'});


fprintf('\n');
fprintf('=========================================================\n');
fprintf('BEST LQR CONFIGURATION BY INDIVIDUAL METRIC\n');
fprintf('=========================================================\n');


disp(best_by_metric);


%% ==========================================================
% LOCAL FUNCTION:
% Equal-Probability Interval Sampling
%% ==========================================================

function [values,interval_index] = ...
    sample_equal_probability_intervals(intervals,N)

% Number of intervals

N_intervals = ...
    size(intervals,1);


% Randomly select one interval for every sample.
%
% Every interval has exactly 1/N_intervals probability.

interval_index = ...
    randi(N_intervals,N,1);


% Allocate output

values = ...
    zeros(N,1);


% Uniformly sample within selected interval

for i = 1:N

    j = ...
        interval_index(i);


    lower_bound = ...
        intervals(j,1);


    upper_bound = ...
        intervals(j,2);


    values(i) = ...
        lower_bound + ...
        (upper_bound-lower_bound)*rand;

end

end


%% ==========================================================
% LOCAL FUNCTION:
% Nonlinear LQR Dynamics
%% ==========================================================

function dx = ...
    lqr_nonlinear_dynamics( ...
    t,x,K, ...
    M,C_v,K_v, ...
    phi_p_0,T_0,V_e_0, ...
    alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E, ...
    C_th,h,tau_e,V_max)

%% ----------------------------------------------------------
% Force State Vector to Column Orientation
%% ----------------------------------------------------------

x = ...
    x(:);


%% ==========================================================
% Physical States
%% ==========================================================

q = ...
    x(1);


q_dot = ...
    x(2);


Phi_p = ...
    x(3);


T = ...
    x(4);


V_e = ...
    x(5);


%% ==========================================================
% Perturbation Coordinates
%% ==========================================================

x_pert = zeros(5,1);


x_pert(1) = ...
    q;


x_pert(2) = ...
    q_dot;


x_pert(3) = ...
    Phi_p - phi_p_0;


x_pert(4) = ...
    T - T_0;


x_pert(5) = ...
    V_e - V_e_0;


%% ==========================================================
% LQR State Feedback
%% ==========================================================

u_cmd = ...
    -K*x_pert;


%% ==========================================================
% Physical Actuator Saturation
%% ==========================================================

u = max( ...
    min(u_cmd,V_max), ...
    -V_max);


u = ...
    u(1);


%% ==========================================================
% Nonlinear Plant Dynamics
%% ==========================================================

dx = ...
    zeros(5,1);


%% ----------------------------------------------------------
% Mechanical Position
%% ----------------------------------------------------------

dx(1) = ...
    q_dot;


%% ----------------------------------------------------------
% Mechanical Velocity
%% ----------------------------------------------------------

dx(2) = ...
    ( ...
    -C_v*q_dot ...
    -K_v*q ...
    +alpha_B*Phi_p^2 ...
    +alpha_E*V_e^2 ...
    -alpha_B*phi_p_0^2 ...
    -alpha_E*V_e_0^2 ...
    ) / M;


%% ----------------------------------------------------------
% Electromagnetic Flux Perturbation Dynamics
%% ----------------------------------------------------------

dx(3) = ...
    -gamma_Phi*(Phi_p-phi_p_0);


%% ----------------------------------------------------------
% Thermal Dynamics
%% ----------------------------------------------------------

dx(4) = ...
    ( ...
    k_E*V_e^2 ...
    +k_Phi*Phi_p^2 ...
    -h*(T-T_0) ...
    -k_E*V_e_0^2 ...
    -k_Phi*phi_p_0^2 ...
    ) / C_th;


%% ----------------------------------------------------------
% Electrical Actuator Dynamics
%% ----------------------------------------------------------

dx(5) = ...
    -(V_e-V_e_0)/tau_e ...
    +u/tau_e;

end


%% ==========================================================
% LOCAL FUNCTION:
% Nonlinear LQG Dynamics
%% ==========================================================

function dz = ...
    lqg_nonlinear_dynamics( ...
    t,z,A,B,C,L,K, ...
    M,C_v,K_v, ...
    phi_p_0,T_0,V_e_0, ...
    alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E, ...
    C_th,h,tau_e,V_max)

%% ----------------------------------------------------------
% Force Combined State Vector to Column Orientation
%% ----------------------------------------------------------

z = ...
    z(:);


%% ==========================================================
% Plant and Observer States
%% ==========================================================

x = ...
    z(1:5);


x_hat = ...
    z(6:10);


%% ==========================================================
% Estimated-State Feedback
%% ==========================================================

u_cmd = ...
    -K*x_hat;


%% ==========================================================
% Physical Actuator Saturation
%% ==========================================================

u = max( ...
    min(u_cmd,V_max), ...
    -V_max);


u = ...
    u(1);


%% ==========================================================
% Nonlinear Plant States
%% ==========================================================

q = ...
    x(1);


q_dot = ...
    x(2);


Phi_p = ...
    x(3);


T = ...
    x(4);


V_e = ...
    x(5);


%% ==========================================================
% Nonlinear Plant Dynamics
%% ==========================================================

x_dot = ...
    zeros(5,1);


%% ----------------------------------------------------------
% Mechanical Position
%% ----------------------------------------------------------

x_dot(1) = ...
    q_dot;


%% ----------------------------------------------------------
% Mechanical Velocity
%% ----------------------------------------------------------

x_dot(2) = ...
    ( ...
    -C_v*q_dot ...
    -K_v*q ...
    +alpha_B*Phi_p^2 ...
    +alpha_E*V_e^2 ...
    -alpha_B*phi_p_0^2 ...
    -alpha_E*V_e_0^2 ...
    ) / M;


%% ----------------------------------------------------------
% Electromagnetic Flux Perturbation Dynamics
%% ----------------------------------------------------------

x_dot(3) = ...
    -gamma_Phi*(Phi_p-phi_p_0);


%% ----------------------------------------------------------
% Thermal Dynamics
%% ----------------------------------------------------------

x_dot(4) = ...
    ( ...
    k_E*V_e^2 ...
    +k_Phi*Phi_p^2 ...
    -h*(T-T_0) ...
    -k_E*V_e_0^2 ...
    -k_Phi*phi_p_0^2 ...
    ) / C_th;


%% ----------------------------------------------------------
% Electrical Actuator Dynamics
%% ----------------------------------------------------------

x_dot(5) = ...
    -(V_e-V_e_0)/tau_e ...
    +u/tau_e;


%% ==========================================================
% Convert Physical Plant State to Perturbation Coordinates
%% ==========================================================

x_pert = ...
    zeros(5,1);


x_pert(1) = ...
    q;


x_pert(2) = ...
    q_dot;


x_pert(3) = ...
    Phi_p-phi_p_0;


x_pert(4) = ...
    T-T_0;


x_pert(5) = ...
    V_e-V_e_0;


%% ==========================================================
% Measurements
%% ==========================================================

y = ...
    C*x_pert;


%% ----------------------------------------------------------
% Estimated Measurements
%% ----------------------------------------------------------

y_hat = ...
    C*x_hat;


%% ==========================================================
% Linear Kalman Observer
%% ==========================================================
%
% The observer operates in perturbation coordinates:
%
%   x_hat_dot = A*x_hat + B*u + L*(y-y_hat)
%
% The same saturated physical input applied to the nonlinear
% plant is supplied to the observer.
%
% ==========================================================

x_hat_dot = ...
    A*x_hat ...
    +B*u ...
    +L*(y-y_hat);


%% ==========================================================
% Combined Plant + Observer Dynamics
%% ==========================================================

dz = ...
    zeros(10,1);


dz(1:5) = ...
    x_dot;


dz(6:10) = ...
    x_hat_dot;

end
