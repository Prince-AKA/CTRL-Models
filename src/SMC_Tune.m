%% =========================================================
% SMC Random Parameter Robustness Study
%
% Allows eta, phi, and Ks to be:
%   - fixed at user-specified values
%   - varied uniformly over user-specified intervals
%
% Exactly 0, 1, or 2 parameters may be fixed.
%
% N is the NUMBER OF PARAMETER REALIZATIONS.
%
% Every parameter realization is simulated over ALL loaded
% initial conditions in ICs.
%
% Therefore:
%
%   Total closed-loop simulations = N * N_IC
%
% Example:
%
%   N    = 100 parameter realizations
%   N_IC = 1000 initial conditions
%
%   Total simulations = 100,000
%
% Each parameter realization receives its own:
%   - eta
%   - phi
%   - Ks
%
% Each individual IC simulation receives its own:
%   - uncertain physical parameters
%   - disturbance realization
%
% The final table displays the TEN BEST parameter
% realizations ranked by robust convergence percentage.
%
% IMPORTANT:
% All nominal model parameters and simulation data required by
% parfor are explicitly packaged into P before entering the
% parallel loop. This prevents worker-side workspace errors such
% as "Unrecognized function or variable 'Cv'."
% =========================================================

%% =========================================================
% LOAD INITIAL CONDITIONS
% ==========================================================

load('MC_initial_conditions_1000.mat');

%% =========================================================
% VERIFY REQUIRED VARIABLES
% ==========================================================

required_vars = { ...
    'M','Cv','Kv','phi_p_0','T_0','Ve_0', ...
    'alpha_B','alpha_E','gamma_Phi','k_Phi','k_E', ...
    'C_th','h','tau_e', ...
    'A','B','C','L', ...
    'q_ref','tol_q','tol_qdot', ...
    't','dt','Nt','Tsim','N_convergence', ...
    'lambda_e','Ve_min','Ve_max', ...
    'd_A1','d_f1','d_A2','d_f2','d_sigma','d_tau', ...
    'ICs'};

missing_vars = {};

for ii = 1:numel(required_vars)

    if ~exist(required_vars{ii},'var')
        missing_vars{end+1} = required_vars{ii}; 
    end

end

if ~isempty(missing_vars)

    error(['The following required variables are missing from the ', ...
           'workspace:\n\n%s\n\n', ...
           'Run your nominal-model initialization script first, ', ...
           'then run this script.'], ...
           strjoin(missing_vars, ', '));

end

%% =========================================================
% VERIFY INITIAL CONDITIONS
% ==========================================================

N_IC = size(ICs,1);

if N_IC < 1
    error('ICs must contain at least one initial condition.');
end

if size(ICs,2) ~= 5
    error('ICs must be an N_IC x 5 matrix.');
end

%% =========================================================
% VERIFY TIME VECTOR
% ==========================================================

if numel(t) ~= Nt
    error('Length of t (%d) does not match Nt (%d).', ...
        numel(t),Nt);
end

t = t(:);

%% =========================================================
% GRAPHICS SETTINGS
% ==========================================================

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% =========================================================
% PACKAGE ALL MODEL/SIMULATION DATA FOR PARFOR
%
% This is the important worker-compatibility fix.
%
% Instead of relying on individual variables such as Cv being
% visible to the parallel workers, everything required by the
% worker is explicitly placed inside P.
% ==========================================================

P = struct();

%% Nominal physical parameters

P.M         = M;
P.Cv        = Cv;
P.Kv        = Kv;

P.phi_p_0   = phi_p_0;
P.T_0       = T_0;
P.Ve_0      = Ve_0;

P.alpha_B   = alpha_B;
P.alpha_E   = alpha_E;
P.gamma_Phi = gamma_Phi;
P.k_Phi     = k_Phi;
P.k_E       = k_E;

P.C_th      = C_th;
P.h         = h;
P.tau_e     = tau_e;

%% Nominal observer/controller model

P.A         = A;
P.B         = B;
P.C         = C;
P.L         = L;

P.q_ref     = q_ref;
P.tol_q     = tol_q;
P.tol_qdot  = tol_qdot;

%% Time settings

P.t             = t;
P.dt            = dt;
P.Nt            = Nt;
P.Tsim          = Tsim;
P.N_convergence = N_convergence;

%% Actuator settings

P.lambda_e  = lambda_e;
P.Ve_min    = Ve_min;
P.Ve_max    = Ve_max;

%% Disturbance parameters

P.d_A1    = d_A1;
P.d_f1    = d_f1;
P.d_A2    = d_A2;
P.d_f2    = d_f2;
P.d_sigma = d_sigma;
P.d_tau   = d_tau;

%% Initial conditions

P.ICs = ICs;

%% =========================================================
% USER PARAMETER SELECTION
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('SMC Random Parameter Robustness Study\n');
fprintf('=========================================================\n');

fprintf('\nHow many SMC parameters should be VARIED?\n');
fprintf('You may vary 1, 2, or 3 parameters.\n\n');

n_var = input( ...
    'Number of varied parameters [1, 2, or 3]: ');

while isempty(n_var) || ...
        ~isscalar(n_var) || ...
        ~isfinite(n_var) || ...
        ~ismember(n_var,[1 2 3])

    fprintf('\nInvalid selection. Enter 1, 2, or 3.\n');

    n_var = input('Number of varied parameters [1, 2, or 3]: ');

end

fprintf('\n');
fprintf('  1 = eta\n');
fprintf('  2 = phi\n');
fprintf('  3 = Ks\n');

%% =========================================================
% SELECT WHICH PARAMETERS ARE FIXED
% ==========================================================

fixed_idx = [];

if n_var < 3

    fprintf('\n');
    fprintf('Enter the index number of the parameter to fix.\n');
    fprintf('Use:\n');
    fprintf('  1 = eta\n');
    fprintf('  2 = phi\n');
    fprintf('  3 = Ks\n\n');

    if n_var == 1
        fprintf('Example: 1\n\n');
    else
        fprintf('Example: [1 3]\n');
        fprintf('This would fix eta and Ks.\n\n');
    end

    fixed_idx = input( ...
        'Fixed parameter(s): ');

    fixed_idx = fixed_idx(:).';

    valid_selection = ...
        isnumeric(fixed_idx) && ...
        numel(fixed_idx) == n_fixed && ...
        all(isfinite(fixed_idx)) && ...
        all(fixed_idx == round(fixed_idx)) && ...
        all(ismember(fixed_idx,[1 2 3])) && ...
        numel(unique(fixed_idx)) == n_fixed;

    while ~valid_selection

        fprintf('\nInvalid selection.\n');
        fprintf('Use 1 = eta, 2 = phi, 3 = Ks.\n');

        if n_fixed == 1
            fprintf('Enter a single number, e.g. 1.\n\n');
        else
            fprintf('Enter two different numbers, e.g. [1 3].\n\n');
        end

        fixed_idx = input( ...
            'Fixed parameter(s): ');

        fixed_idx = fixed_idx(:).';

        valid_selection = ...
            isnumeric(fixed_idx) && ...
            numel(fixed_idx) == n_fixed && ...
            all(isfinite(fixed_idx)) && ...
            all(fixed_idx == round(fixed_idx)) && ...
            all(ismember(fixed_idx,[1 2 3])) && ...
            numel(unique(fixed_idx)) == n_fixed;

    end

end

%% =========================================================
% PARAMETER NAMES
% ==========================================================

param_names = {'eta','phi','Ks'};

is_fixed = false(1,3);

is_fixed(fixed_idx) = true;

%% =========================================================
% FIXED PARAMETER VALUES
% ==========================================================

eta_fixed = NaN;
phi_fixed = NaN;
Ks_fixed  = NaN;

for i = fixed_idx

    switch i

        case 1

            eta_fixed = input( ...
                '\nEnter fixed eta value: ');

            while isempty(eta_fixed) || ...
                    ~isscalar(eta_fixed) || ...
                    ~isfinite(eta_fixed)

                fprintf('Invalid eta value.\n');

                eta_fixed = input( ...
                    'Enter fixed eta value: ');

            end

        case 2

            phi_fixed = input( ...
                '\nEnter fixed phi value: ');

            while isempty(phi_fixed) || ...
                    ~isscalar(phi_fixed) || ...
                    ~isfinite(phi_fixed) || ...
                    phi_fixed <= 0

                fprintf('Invalid phi value. phi must be > 0.\n');

                phi_fixed = input( ...
                    'Enter fixed phi value: ');

            end

        case 3

            Ks_fixed = input( ...
                '\nEnter fixed Ks value: ');

            while isempty(Ks_fixed) || ...
                    ~isscalar(Ks_fixed) || ...
                    ~isfinite(Ks_fixed)

                fprintf('Invalid Ks value.\n');

                Ks_fixed = input( ...
                    'Enter fixed Ks value: ');

            end

    end

end

%% =========================================================
% BOUNDS FOR VARYING PARAMETERS
% ==========================================================

eta_low  = NaN;
eta_high = NaN;

phi_low  = NaN;
phi_high = NaN;

Ks_low   = NaN;
Ks_high  = NaN;

for i = 1:3

    if ~is_fixed(i)

        switch i

            case 1

                fprintf('\n-----------------------------------------\n');
                fprintf('eta will be VARIED.\n');
                fprintf('-----------------------------------------\n');

                eta_low = input('Lower eta bound: ');
                eta_high = input('Upper eta bound: ');

                while isempty(eta_low) || ...
                        isempty(eta_high) || ...
                        ~isscalar(eta_low) || ...
                        ~isscalar(eta_high) || ...
                        ~isfinite(eta_low) || ...
                        ~isfinite(eta_high) || ...
                        eta_high <= eta_low

                    fprintf('\nInvalid eta interval.\n');
                    fprintf('Upper bound must be greater than lower bound.\n\n');

                    eta_low = input('Lower eta bound: ');
                    eta_high = input('Upper eta bound: ');

                end

            case 2

                fprintf('\n-----------------------------------------\n');
                fprintf('phi will be VARIED.\n');
                fprintf('-----------------------------------------\n');

                phi_low = input('Lower phi bound: ');
                phi_high = input('Upper phi bound: ');

                while isempty(phi_low) || ...
                        isempty(phi_high) || ...
                        ~isscalar(phi_low) || ...
                        ~isscalar(phi_high) || ...
                        ~isfinite(phi_low) || ...
                        ~isfinite(phi_high) || ...
                        phi_low <= 0 || ...
                        phi_high <= phi_low

                    fprintf('\nInvalid phi interval.\n');
                    fprintf('Require 0 < lower bound < upper bound.\n\n');

                    phi_low = input('Lower phi bound: ');
                    phi_high = input('Upper phi bound: ');

                end

            case 3

                fprintf('\n-----------------------------------------\n');
                fprintf('Ks will be VARIED.\n');
                fprintf('-----------------------------------------\n');

                Ks_low = input('Lower Ks bound: ');
                Ks_high = input('Upper Ks bound: ');

                while isempty(Ks_low) || ...
                        isempty(Ks_high) || ...
                        ~isscalar(Ks_low) || ...
                        ~isscalar(Ks_high) || ...
                        ~isfinite(Ks_low) || ...
                        ~isfinite(Ks_high) || ...
                        Ks_high <= Ks_low

                    fprintf('\nInvalid Ks interval.\n');
                    fprintf('Upper bound must be greater than lower bound.\n\n');

                    Ks_low = input('Lower Ks bound: ');
                    Ks_high = input('Upper Ks bound: ');

                end

        end

    end

end

%% =========================================================
% NUMBER OF PARAMETER REALIZATIONS
% ==========================================================

fprintf('\n');

N = input( ...
    'Number of parameter realizations N: ');

while isempty(N) || ...
        ~isscalar(N) || ...
        ~isfinite(N) || ...
        N < 1 || ...
        N ~= round(N)

    fprintf('\nInvalid N. Enter a positive integer.\n');

    N = input( ...
        'Number of parameter realizations N: ');

end

N = round(N);

%% =========================================================
% TOTAL SIMULATION COUNT
% ==========================================================

total_simulations = N * N_IC;

%% =========================================================
% GENERATE SMC PARAMETER REALIZATIONS
% ==========================================================

eta_samples = zeros(N,1);
phi_samples = zeros(N,1);
Ks_samples  = zeros(N,1);

%% eta

if is_fixed(1)

    eta_samples(:) = eta_fixed;

else

    eta_samples = ...
        eta_low + ...
        (eta_high - eta_low).*rand(N,1);

end

%% phi

if is_fixed(2)

    phi_samples(:) = phi_fixed;

else

    phi_samples = ...
        phi_low + ...
        (phi_high - phi_low).*rand(N,1);

end

%% Ks

if is_fixed(3)

    Ks_samples(:) = Ks_fixed;

else

    Ks_samples = ...
        Ks_low + ...
        (Ks_high - Ks_low).*rand(N,1);

end

%% =========================================================
% DISPLAY EXPERIMENT CONFIGURATION
% ==========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('SMC Parameter Study Configuration\n');
fprintf('=========================================================\n');

fprintf('Parameter realizations: %d\n',N);
fprintf('Loaded initial conditions: %d\n',N_IC);
fprintf('Total closed-loop simulations: %d\n\n', ...
    total_simulations);

fprintf('Each parameter realization will be tested against\n');
fprintf('ALL %d loaded initial conditions.\n',N_IC);

fprintf('\neta:\n');

if is_fixed(1)

    fprintf('  FIXED   = %.6g\n',eta_fixed);

else

    fprintf('  VARYING = [%.6g, %.6g]\n', ...
        eta_low,eta_high);

end

fprintf('\nphi:\n');

if is_fixed(2)

    fprintf('  FIXED   = %.6g\n',phi_fixed);

else

    fprintf('  VARYING = [%.6g, %.6g]\n', ...
        phi_low,phi_high);

end

fprintf('\nKs:\n');

if is_fixed(3)

    fprintf('  FIXED   = %.6g\n',Ks_fixed);

else

    fprintf('  VARYING = [%.6g, %.6g]\n', ...
        Ks_low,Ks_high);

end

fprintf('\n=========================================================\n');
fprintf('Beginning %d parameter realizations (%d total simulations)...\n', ...
    N,total_simulations);
fprintf('=========================================================\n\n');

%% =========================================================
% PREALLOCATE RESULTS
%
% Each row corresponds to ONE parameter realization.
%
% Metrics are aggregated over ALL initial conditions.
% ==========================================================

convergence_percentage = zeros(N,1);
num_converged          = zeros(N,1);

mean_energy            = zeros(N,1);
mean_maxV              = zeros(N,1);
saturation_incidence   = zeros(N,1);
mean_sat_fraction      = zeros(N,1);

%% =========================================================
% START PARALLEL POOL
%
% Explicitly make sure the current script/function is available
% to the workers.
% ==========================================================

pool = gcp('nocreate');

if isempty(pool)

    pool = parpool;

end

%% =========================================================
% MONTE CARLO SIMULATIONS
% ==========================================================

parfor j = 1:N

    %% -----------------------------------------------------
    % SMC parameters for THIS parameter realization
    % -----------------------------------------------------

    eta_j = eta_samples(j);
    phi_j = phi_samples(j);
    Ks_j  = Ks_samples(j);

    %% -----------------------------------------------------
    % Local storage for ALL initial conditions
    % -----------------------------------------------------

    success_ic      = false(N_IC,1);
    energy_ic       = zeros(N_IC,1);
    maxV_ic         = zeros(N_IC,1);
    saturated_ic    = false(N_IC,1);
    sat_fraction_ic = zeros(N_IC,1);

    %% =====================================================
    % Run this parameter realization over ALL ICs
    % =====================================================

    for i = 1:N_IC

        %% -------------------------------------------------
        % Random uncertain physical parameters
        % -------------------------------------------------

        M_j = P.M * ...
            (1 + 0.80*(2*rand-1));

        Cv_j = P.Cv * ...
            (1 + 1.00*(2*rand-1));

        Kv_j = P.Kv * ...
            (1 + 0.80*(2*rand-1));

        alpha_B_j = P.alpha_B * ...
            (1 + 0.70*(2*rand-1));

        alpha_E_j = P.alpha_E * ...
            (1 + 0.70*(2*rand-1));

        gamma_j = P.gamma_Phi * ...
            (1 + 0.70*(2*rand-1));

        kPhi_j = P.k_Phi * ...
            (1 + 0.70*(2*rand-1));

        kE_j = P.k_E * ...
            (1 + 0.70*(2*rand-1));

        h_j = P.h * ...
            (1 + 0.70*(2*rand-1));

        tau_j = P.tau_e * ...
            (1 + 0.70*(2*rand-1));

        %% -------------------------------------------------
        % Uncertain plant
        % -------------------------------------------------

        A_j = [ ...
            0, 1, 0, 0, 0;
            -Kv_j/M_j, -Cv_j/M_j, ...
                2*alpha_B_j*P.phi_p_0/M_j, 0, ...
                2*alpha_E_j*P.Ve_0/M_j;
            0, 0, -gamma_j, 0, 0;
            0, 0, ...
                2*kPhi_j*P.phi_p_0/P.C_th, ...
                -h_j/P.C_th, ...
                2*kE_j*P.Ve_0/P.C_th;
            0, 0, 0, 0, -1/tau_j];

        B_j = [ ...
            0;
            0;
            0;
            0;
            1/tau_j];

        %% -------------------------------------------------
        % Initial condition
        % -------------------------------------------------

        dx0_j = P.ICs(i,:).';

        %% -------------------------------------------------
        % Disturbance
        % -------------------------------------------------

        d_profile_j = generate_disturbance( ...
            P.t, ...
            P.dt, ...
            P.d_A1, ...
            P.d_f1, ...
            P.d_A2, ...
            P.d_f2, ...
            P.d_sigma, ...
            P.d_tau);

        %% -------------------------------------------------
        % Initial states
        % -------------------------------------------------

        x = dx0_j;

        xhat = zeros(5,1);

        X       = zeros(P.Nt,1);
        Xdot    = zeros(P.Nt,1);
        U       = zeros(P.Nt,1);
        Ve_hist = zeros(P.Nt,1);

        sat_count = 0;

        %% =================================================
        % CLOSED-LOOP SIMULATION
        % =================================================

        for k = 1:P.Nt

            %% Measurement

            y = P.C*x;

            %% Estimated states

            qhat    = xhat(1);
            qdothat = xhat(2);
            Phihat  = xhat(3);

            %% Sliding surface

            error_hat = qhat - P.q_ref;

            s = qdothat + eta_j*error_hat;

            %% Boundary layer

            sat_s = max(-1,min(1,s/phi_j));

            %% Switching control

            usw = Ks_j*sat_s;

            %% Desired acceleration

            qddot_des = ...
                -eta_j*qdothat ...
                -usw;

            %% Desired physical voltage

            dVe_des = ...
                (P.M/(2*P.alpha_E*P.Ve_0)) * ...
                ( ...
                    qddot_des ...
                    + (P.Cv/P.M)*qdothat ...
                    + (P.Kv/P.M)*qhat ...
                    - (2*P.alpha_B*P.phi_p_0/P.M)*Phihat ...
                );

            Ve_des = P.Ve_0 + dVe_des;

            %% Current actuator state

            Ve = P.Ve_0 + x(5);

            %% Physical actuator command

            Ve_cmd = ...
                Ve + ...
                tau_j*P.lambda_e*(Ve_des - Ve);

            %% Saturation detection

            if Ve_cmd < P.Ve_min || ...
                    Ve_cmd > P.Ve_max

                sat_count = sat_count + 1;

            end

            %% Apply saturation

            Ve_cmd = ...
                max(P.Ve_min, ...
                min(P.Ve_max,Ve_cmd));

            %% Perturbation-coordinate input

            u = Ve_cmd - P.Ve_0;

            %% Disturbance

            d = d_profile_j(k);

            %% Uncertain plant

            B_d_j = ...
                [0;1/M_j;0;0;0];

            x_dot = ...
                A_j*x + ...
                B_j*u + ...
                B_d_j*d;

            x = ...
                x + P.dt*x_dot;

            %% Common nominal observer

            xhat_dot = ...
                P.A*xhat + ...
                P.B*u + ...
                P.L*(y - P.C*xhat);

            xhat = ...
                xhat + P.dt*xhat_dot;

            %% Store

            Ve = P.Ve_0 + x(5);

            X(k)       = x(1);
            Xdot(k)    = x(2);
            U(k)       = u;
            Ve_hist(k) = Ve;

        end

        %% =================================================
        % METRICS FOR THIS INITIAL CONDITION
        % =================================================

        energy_ic(i) = ...
            trapz(P.t,U.^2);

        maxV_ic(i) = ...
            max(Ve_hist);

        saturated_ic(i) = ...
            sat_count > 0;

        sat_fraction_ic(i) = ...
            sat_count/P.Nt;

        %% =================================================
        % ROBUST CONVERGENCE FOR THIS INITIAL CONDITION
        % =================================================

        idx_start = ...
            max(1,P.Nt-P.N_convergence+1);

        q_tail = ...
            X(idx_start:end);

        qdot_tail = ...
            Xdot(idx_start:end);

        q_rms_tail = ...
            rms(q_tail - P.q_ref);

        qdot_rms_tail = ...
            rms(qdot_tail);

        success_ic(i) = ...
            (q_rms_tail < P.tol_q) && ...
            (qdot_rms_tail < P.tol_qdot);

    end

    %% =====================================================
    % AGGREGATE ALL ICs FOR THIS PARAMETER REALIZATION
    % =====================================================

    num_converged(j) = ...
        sum(success_ic);

    convergence_percentage(j) = ...
        100*mean(success_ic);

    mean_energy(j) = ...
        mean(energy_ic);

    mean_maxV(j) = ...
        mean(maxV_ic);

    saturation_incidence(j) = ...
        100*mean(saturated_ic);

    mean_sat_fraction(j) = ...
        100*mean(sat_fraction_ic);

end

%% =========================================================
% RESULTS TABLE
%
% ONE ROW = ONE PARAMETER REALIZATION
%
% Metrics are aggregated over ALL loaded ICs.
% ==========================================================

results_table = table( ...
    (1:N).', ...
    eta_samples, ...
    phi_samples, ...
    Ks_samples, ...
    num_converged, ...
    convergence_percentage, ...
    mean_energy, ...
    mean_maxV, ...
    saturation_incidence, ...
    mean_sat_fraction, ...
    'VariableNames', ...
    {'Run','eta','phi','Ks', ...
     'ConvergedICs','ConvergencePercentage', ...
     'MeanEnergy','MeanMaxVoltage', ...
     'SaturationIncidence','MeanSaturationFraction'});

%% =========================================================
% SORT BY CONVERGENCE PERCENTAGE
% ==========================================================

results_sorted = ...
    sortrows( ...
        results_table, ...
        'ConvergencePercentage', ...
        'descend');

%% =========================================================
% BEST 10 PARAMETER REALIZATIONS
% ==========================================================

n_best = min(10,N);

best_results = ...
    results_sorted(1:n_best,:);

%% =========================================================
% OVERALL SUMMARY
% ==========================================================

fprintf('\n=========================================================\n');
fprintf('SMC Random Parameter Study Results\n');
fprintf('=========================================================\n');

fprintf('\nParameter realizations: %d\n',N);

fprintf('Initial conditions per realization: %d\n',N_IC);

fprintf('Total closed-loop simulations: %d\n', ...
    total_simulations);

fprintf('\nTotal successful closed-loop simulations: %d / %d\n', ...
    sum(num_converged), ...
    total_simulations);

fprintf('Overall convergence rate across all simulations: %.2f %%\n', ...
    100*sum(num_converged)/total_simulations);

fprintf('\nMean control energy across parameter realizations: %.6g\n', ...
    mean(mean_energy));

fprintf('Mean maximum physical voltage: %.6g\n', ...
    mean(mean_maxV));

fprintf('Mean saturation incidence: %.2f %%\n', ...
    mean(saturation_incidence));

fprintf('Mean saturation duration: %.2f %%\n', ...
    mean(mean_sat_fraction));

%% =========================================================
% BEST 10
% ==========================================================

fprintf('\n=========================================================\n');
fprintf('TOP %d PARAMETER REALIZATIONS\n',n_best);
fprintf('Ranked by Convergence Percentage\n');
fprintf('=========================================================\n\n');

disp(best_results);

%% =========================================================
% OPTIONAL SAVE
% ==========================================================

% save('SMC_random_parameter_results.mat', ...
%     'results_table', ...
%     'results_sorted', ...
%     'best_results');

%% =========================================================
% LOCAL DISTURBANCE FUNCTION
%
% Keep this function in the SAME FILE as SMC_Tune.m.
% This makes the disturbance generator available to the
% parallel workers without relying on a separate file.
% ==========================================================

function d_profile = generate_disturbance( ...
    t,dt,d_A1,d_f1,d_A2,d_f2,d_sigma,d_tau)

    t = t(:);

    d_profile = ...
        d_A1*sin(2*pi*d_f1*t) + ...
        d_A2*sin(2*pi*d_f2*t);

    % Gaussian envelope/modulation
    if ~isempty(d_sigma) && ...
            ~isempty(d_tau) && ...
            isfinite(d_sigma) && ...
            isfinite(d_tau) && ...
            d_sigma > 0

        envelope = ...
            exp(-0.5*((t-d_tau)/d_sigma).^2);

        d_profile = ...
            d_profile .* envelope;

    end

end