%% Script 3: PID Robustness Study

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; clc; close all;

%% ============================================================
% Model Parameters
%% ============================================================

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

%% ============================================================
% PID Gains
%% ============================================================

Kp = 5.8;
Ki = 5.3;
Kd = 5.2;

%% ============================================================
% Actuator Limits
%% ============================================================

u_min = 0.0;
u_max = 10.0;

%% ============================================================
% Disturbance Parameters
%% ============================================================
%
% The disturbance is
%
% d(t) =
%     A1*sin(2*pi*f1*t)
%     + A2*sin(2*pi*f2*t)
%     + d_colored(t)
%
% where d_colored is a first-order Ornstein-Uhlenbeck
% colored stochastic process with stationary RMS sigma_d.
%
% The selected values give an expected disturbance RMS of
%
% sqrt((A1^2 + A2^2)/2 + sigma_d^2)
%
% approximately equal to 0.044.
%
% The same disturbance realization is applied to the
% uncontrolled and PID-controlled plants within each
% Monte Carlo realization.

d_A1    = 0.05;
d_f1    = 0.5;

d_A2    = 0.025;
d_f2    = 1.2;

d_sigma = 0.02;
d_tau   = 0.15;

%% ============================================================
% Monte Carlo Initial-Condition Study
%% ============================================================
num_IC = 1000;
t_final = 10;
tspan = linspace(0,t_final,1000);

% Initial-condition bounds
q_max      = 0.8;
q_dot_max  = 0.4;
dphi_max   = 0.3;
dT_max     = 0.24;
dV_e_max   = 0.4;

%% ============================================================
% Generate or Load Fixed Initial-Condition Ensemble
%% ============================================================
IC_file = 'MC_initial_conditions_1000.mat';

if isfile(IC_file)
    load(IC_file, 'ICs');
    fprintf('Loaded fixed Monte Carlo IC ensemble: %s\n', IC_file);
else
    rng(7);

    ICs = [ ...
        q_max     * (2*rand(num_IC,1)-1), ...
        q_dot_max * (2*rand(num_IC,1)-1), ...
        dphi_max  * (2*rand(num_IC,1)-1), ...
        dT_max    * (2*rand(num_IC,1)-1), ...
        dV_e_max  * (2*rand(num_IC,1)-1) ];

    save(IC_file, 'ICs');

    fprintf('Generated and saved fixed Monte Carlo IC ensemble: %s\n', IC_file);
end

%% ============================================================
% Convergence Criteria
%% ============================================================
eps_q     = 1e-2;
eps_q_dot = 1e-2;
eps_norm  = 2e-2;

%% ============================================================
% Preallocate Results
% ============================================================

converged_pid = false(num_IC,1);
converged_nl  = false(num_IC,1);

settle_times_pid   = nan(num_IC,1);
control_energy_pid = nan(num_IC,1);
max_Ve_pid         = nan(num_IC,1);
disturbance_rms    = nan(num_IC,1);

%% ============================================================
% Monte Carlo Initial-Condition and Disturbance Study
%% ============================================================

figure;
hold on;
grid on;

for i = 1:num_IC

    %% --------------------------------------------------------
    % Retrieve Initial Condition
    % ---------------------------------------------------------

    dx_0 = ICs(i,:).';

    %% --------------------------------------------------------
    % Construct Nonlinear Plant Initial Condition
    %% ---------------------------------------------------------

    x_0_nl = [
        dx_0(1);
        dx_0(2);
        phi_p_0 + dx_0(3);
        T_0     + dx_0(4);
        V_e_0   + dx_0(5)
    ];

    %% --------------------------------------------------------
    % PID Augmented State
    %
    % z = [q; qdot; Phi_p; T; Ve; I_e]
    %
    % I_e = integral of tracking error
    % ---------------------------------------------------------

    z_0_pid = [x_0_nl; 0];

    %% --------------------------------------------------------
    % Generate Disturbance Realization
    %% ---------------------------------------------------------

    dt = tspan(2) - tspan(1);

    d_profile = generate_disturbance( ...
        tspan,dt, ...
        d_A1,d_f1, ...
        d_A2,d_f2, ...
        d_sigma,d_tau);

    disturbance_rms(i) = rms(d_profile);

    %% ========================================================
    % Uncontrolled Nonlinear Plant
    %% ========================================================

    f_nl = @(t,x) [
        x(2);

        ( ...
        -C_v*x(2) ...
        -K_v*x(1) ...
        +alpha_B*x(3)^2 ...
        +alpha_E*x(5)^2 ...
        -alpha_B*phi_p_0^2 ...
        -alpha_E*V_e_0^2 ...
        +interp1(tspan,d_profile,t,'linear','extrap') ...
        ) / M;

        -gamma_Phi*(x(3) - phi_p_0);

        ( ...
        k_E*x(5)^2 ...
        +k_Phi*x(3)^2 ...
        -h*(x(4) - T_0) ...
        -k_E*V_e_0^2 ...
        -k_Phi*phi_p_0^2 ...
        ) / C_th;

        -(x(5) - V_e_0)/tau_e
    ];

    [t_nl,x_nl] = ode45( ...
        f_nl,tspan,x_0_nl);

    %% --------------------------------------------------------
    % Uncontrolled Convergence
    %% ---------------------------------------------------------

    q_nl = x_nl(:,1);
    qdot_nl = x_nl(:,2);

    q_nl_final = q_nl(end);
    qdot_nl_final = qdot_nl(end);

    terminal_norm_nl = ...
        norm([q_nl_final;qdot_nl_final],2);

    converged_nl(i) = ...
        abs(q_nl_final) < eps_q && ...
        abs(qdot_nl_final) < eps_q_dot && ...
        terminal_norm_nl < eps_norm;

    %% ========================================================
    % PID-Controlled Nonlinear Plant
    %% ========================================================

    f_pid = @(t,z) pid_dynamics( ...
        t,z, ...
        M,C_v,K_v,alpha_B,alpha_E, ...
        gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
        phi_p_0,T_0,V_e_0, ...
        Kp,Ki,Kd,u_min,u_max, ...
        tspan,d_profile);

    [t_pid,z_pid] = ode45( ...
        f_pid,tspan,z_0_pid);

    %% --------------------------------------------------------
    % Extract PID States
    % ---------------------------------------------------------

    q_pid = z_pid(:,1);
    qdot_pid = z_pid(:,2);
    Ve_pid = z_pid(:,5);

    %% --------------------------------------------------------
    % PID Convergence
    %% ---------------------------------------------------------

    q_pid_final = q_pid(end);
    qdot_pid_final = qdot_pid(end);

    terminal_norm_pid = ...
        norm([q_pid_final;qdot_pid_final],2);

    converged_pid(i) = ...
        abs(q_pid_final) < eps_q && ...
        abs(qdot_pid_final) < eps_q_dot && ...
        terminal_norm_pid < eps_norm;

    %% ========================================================
    % Settling Time
    %% ========================================================
    %
    % A trajectory is considered settled only if it enters
    % the convergence region and remains there for the rest
    % of the simulation.
    % ========================================================

    err_norm = vecnorm([q_pid,qdot_pid],2,2);

    outside_after = ...
        flipud(cummax(flipud(err_norm >= eps_norm)));

    idx_settle = find( ...
        err_norm < eps_norm & ...
        outside_after == 0, ...
        1,'first');

    if ~isempty(idx_settle)
        settle_times_pid(i) = t_pid(idx_settle);
    end

    %% ========================================================
    % Control Metrics
    %% ========================================================

    control_energy_pid(i) = ...
        trapz(t_pid,(Ve_pid - V_e_0).^2);

    max_Ve_pid(i) = ...
        max(abs(Ve_pid));

    %% ========================================================
    % Plot Representative Trajectories
    %% ========================================================

    if i <= 10

        plot(q_nl,qdot_nl,'r');

        plot(q_pid,qdot_pid,'g');

    end

end
%{
%% PID Gain Monte Carlo Sweep
%
% Each gain set is evaluated using the SAME initial conditions and
% SAME disturbance realization for every Monte Carlo realization.
% This provides a controlled comparison between PID gain sets.


rng(17);

num_gain_sets = 10;

%% Gain ranges

Kp_range = [5.6, 6.4];
Ki_range = [4.7, 5.5];
Kd_range = [4.8, 6];

%% Random gain realizations

Kp_sweep = Kp_range(1) + ...
    (Kp_range(2)-Kp_range(1))*rand(num_gain_sets,1);

Ki_sweep = Ki_range(1) + ...
    (Ki_range(2)-Ki_range(1))*rand(num_gain_sets,1);

Kd_sweep = Kd_range(1) + ...
    (Kd_range(2)-Kd_range(1))*rand(num_gain_sets,1);

%% ============================================================
% Generate disturbance realizations ONCE
% ============================================================
%
% The same disturbance realization is used for every gain set
% for a given initial condition.
%
% This ensures that the gain sweep compares controllers under
% identical uncertainty realizations.
% ============================================================

dt = tspan(2) - tspan(1);

disturbance_profiles = zeros(num_IC,length(tspan));

for i = 1:num_IC

    disturbance_profiles(i,:) = generate_disturbance( ...
        tspan,dt, ...
        d_A1,d_f1, ...
        d_A2,d_f2, ...
        d_sigma,d_tau);

end

%% Storage

gain_success_rate = zeros(num_gain_sets,1);

gain_mean_settling = nan(num_gain_sets,1);

gain_mean_energy = nan(num_gain_sets,1);

gain_mean_max_Ve = nan(num_gain_sets,1);

%% ============================================================
% Progress information
%% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('PID Gain Monte Carlo Sweep\n');
fprintf('============================================\n');

fprintf('Gain realizations: %d\n',num_gain_sets);

fprintf('Initial conditions per gain set: %d\n',num_IC);

fprintf('Total closed-loop simulations: %d\n', ...
    num_gain_sets*num_IC);

fprintf('Kp range: [%.2f, %.2f]\n', ...
    Kp_range(1),Kp_range(2));

fprintf('Ki range: [%.2f, %.2f]\n', ...
    Ki_range(1),Ki_range(2));

fprintf('Kd range: [%.2f, %.2f]\n', ...
    Kd_range(1),Kd_range(2));

fprintf('Disturbance RMS: %.4f\n', ...
    mean(rms(disturbance_profiles,2)));

fprintf('============================================\n\n');

%% ============================================================
% Gain sweep
%% ============================================================

for g = 1:num_gain_sets

    %% Current gain set

    Kp_g = Kp_sweep(g);
    Ki_g = Ki_sweep(g);
    Kd_g = Kd_sweep(g);

    %% Per-realization storage

    converged_g = false(num_IC,1);

    settling_g = nan(num_IC,1);

    energy_g = nan(num_IC,1);

    max_Ve_g = nan(num_IC,1);

    %% ========================================================
    % Monte Carlo realizations
    %% ========================================================

    for i = 1:num_IC

        %% Retrieve initial condition

        dx_0 = IC_samples(i,:).';

        %% Construct nonlinear plant initial condition

        x_0_nl = [
            dx_0(1);
            dx_0(2);
            phi_p_0 + dx_0(3);
            T_0     + dx_0(4);
            V_e_0   + dx_0(5)
        ];

        %% PID augmented state
        %
        % z = [q; qdot; Phi_p; T; Ve; I_e]

        z_0_pid = [x_0_nl; 0];

        %% Retrieve FIXED disturbance realization

        d_profile = disturbance_profiles(i,:);

        %% PID controlled nonlinear plant

        f_pid = @(t,z) pid_dynamics( ...
            t,z, ...
            M,C_v,K_v,alpha_B,alpha_E, ...
            gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
            phi_p_0,T_0,V_e_0, ...
            Kp_g,Ki_g,Kd_g,u_min,u_max, ...
            tspan,d_profile);

        [t_pid,z_pid] = ode45( ...
            f_pid,tspan,z_0_pid);

        %% Extract states

        q_pid = z_pid(:,1);

        qdot_pid = z_pid(:,2);

        Ve_pid = z_pid(:,5);

        %% ====================================================
        % Convergence
        %% ====================================================

        q_final = q_pid(end);

        qdot_final = qdot_pid(end);

        terminal_norm = ...
            norm([q_final;qdot_final],2);

        converged_g(i) = ...
            abs(q_final) < eps_q && ...
            abs(qdot_final) < eps_q_dot && ...
            terminal_norm < eps_norm;

        %% ====================================================
        % Settling time
        %% ====================================================

        err_norm = ...
            vecnorm([q_pid qdot_pid],2,2);

        idx = find( ...
            err_norm < eps_norm, ...
            1,'first');

        if ~isempty(idx)

            settling_g(i) = t_pid(idx);

        end

        %% ====================================================
        % Control energy
        %% ====================================================

        energy_g(i) = ...
            trapz( ...
                t_pid, ...
                (Ve_pid - V_e_0).^2);

        %% ====================================================
        % Maximum actuator state
        %% ====================================================

        max_Ve_g(i) = max(abs(Ve_pid));

    end

    %% ========================================================
    % Aggregate results
    %% ========================================================

    gain_success_rate(g) = ...
        100*mean(converged_g);

    gain_mean_settling(g) = ...
        mean(settling_g,'omitnan');

    gain_mean_energy(g) = ...
        mean(energy_g,'omitnan');

    gain_mean_max_Ve(g) = ...
        mean(max_Ve_g,'omitnan');

    %% Progress output

    fprintf( ...
        ['Gain set %3d/%3d: ' ...
         'Kp = %.3f, Ki = %.3f, Kd = %.3f, ' ...
         'Success = %.1f%%\n'], ...
        g,num_gain_sets, ...
        Kp_g,Ki_g,Kd_g, ...
        gain_success_rate(g));

end

%% ============================================================
% Gain Sweep Results Table
%% ============================================================

gain_results = table( ...
    Kp_sweep, ...
    Ki_sweep, ...
    Kd_sweep, ...
    gain_success_rate, ...
    gain_mean_settling, ...
    gain_mean_energy, ...
    gain_mean_max_Ve, ...
    'VariableNames', { ...
    'Kp', ...
    'Ki', ...
    'Kd', ...
    'SuccessRate', ...
    'MeanSettlingTime', ...
    'MeanEnergy', ...
    'MeanMaxVoltage'});

%% Sort primarily by success rate,
% then by control energy

gain_results = sortrows( ...
    gain_results, ...
    {'SuccessRate','MeanEnergy'}, ...
    {'descend','ascend'});

%% Display best 20 gain sets

fprintf('\n');
fprintf('============================================\n');
fprintf('Best PID Gain Sets\n');
fprintf('============================================\n');

disp(gain_results( ...
    1:min(20,height(gain_results)),:));
%}

%% ============================================================
% Plot Formatting
%% ============================================================

xlabel('$q$');
ylabel('$\dot{q}$');

legend( ...
    'Nonlinear Uncontrolled', ...
    'PID Controlled', ...
    'Location','best');

%% ============================================================
% Report Results
%% ============================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf('PID Robustness Study\n');
fprintf('=========================================================\n');

fprintf('\nMonte Carlo parameters:\n');
fprintf('Number of ICs tested: %d\n',num_IC);
fprintf('Simulation time: %.2f s\n',t_final);

fprintf('\nPID gains:\n');
fprintf('Kp = %.3f\n',Kp);
fprintf('Ki = %.3f\n',Ki);
fprintf('Kd = %.3f\n',Kd);

fprintf('\nActuator limits:\n');
fprintf('u_min = %.3f\n',u_min);
fprintf('u_max = %.3f\n',u_max);

fprintf('\nInitial-condition bounds:\n');
fprintf('q      in [-%.3f, %.3f]\n',q_max,q_max);
fprintf('qdot   in [-%.3f, %.3f]\n',q_dot_max,q_dot_max);
fprintf('dPhi   in [-%.3f, %.3f]\n',dphi_max,dphi_max);
fprintf('dT     in [-%.3f, %.3f]\n',dT_max,dT_max);
fprintf('dVe    in [-%.3f, %.3f]\n',dV_e_max,dV_e_max);

fprintf('\nDisturbance parameters:\n');
fprintf('A1 = %.3f\n',d_A1);
fprintf('f1 = %.3f Hz\n',d_f1);
fprintf('A2 = %.3f\n',d_A2);
fprintf('f2 = %.3f Hz\n',d_f2);
fprintf('sigma_d = %.3f\n',d_sigma);
fprintf('tau_d = %.3f s\n',d_tau);

fprintf('\nExpected disturbance RMS: %.4f\n', ...
    sqrt((d_A1^2 + d_A2^2)/2 + d_sigma^2));

fprintf('Mean disturbance RMS: %.4f\n', ...
    mean(disturbance_rms,'omitnan'));

fprintf('\nConvergence criterion:\n');
fprintf('|q(tf)| < %.3e\n',eps_q);
fprintf('|qdot(tf)| < %.3e\n',eps_q_dot);
fprintf('norm([q(tf); qdot(tf)]) < %.3e\n',eps_norm);

fprintf('\nUncontrolled nonlinear convergence rate: %.1f%%\n', ...
    100*mean(converged_nl));

fprintf('PID controlled convergence rate: %.1f%%\n', ...
    100*mean(converged_pid));

fprintf('\nPID performance metrics:\n');

fprintf('Mean settling time: %.3f s\n', ...
    mean(settle_times_pid,'omitnan'));

fprintf('Mean control energy: %.3f\n', ...
    mean(control_energy_pid,'omitnan'));

fprintf('Mean max |Ve|: %.3f\n', ...
    mean(max_Ve_pid,'omitnan'));

fprintf('Mean disturbance RMS: %.4f\n', ...
    mean(disturbance_rms,'omitnan'));

fprintf('\n=========================================================\n');

%% ============================================================
% Representative Single-Trajectory Time Response
%% ============================================================

dx_0 = [
    0.20;
    0.00;
    0.05;
    0.03;
    0.10
];

x_0_nl = [
    dx_0(1);
    dx_0(2);
    phi_p_0 + dx_0(3);
    T_0     + dx_0(4);
    V_e_0   + dx_0(5)
];

z_0_pid = [x_0_nl;0];

%% ------------------------------------------------------------
% Generate Representative Disturbance
%% ------------------------------------------------------------

dt = tspan(2) - tspan(1);

d_profile_rep = generate_disturbance( ...
    tspan,dt, ...
    d_A1,d_f1, ...
    d_A2,d_f2, ...
    d_sigma,d_tau);

%% ------------------------------------------------------------
% Representative Uncontrolled Plant
%% ------------------------------------------------------------

f_nl_rep = @(t,x) [
    x(2);

    ( ...
    -C_v*x(2) ...
    -K_v*x(1) ...
    +alpha_B*x(3)^2 ...
    +alpha_E*x(5)^2 ...
    -alpha_B*phi_p_0^2 ...
    -alpha_E*V_e_0^2 ...
    +interp1( ...
        tspan,d_profile_rep,t, ...
        'linear','extrap') ...
    ) / M;

    -gamma_Phi*(x(3) - phi_p_0);

    ( ...
    k_E*x(5)^2 ...
    +k_Phi*x(3)^2 ...
    -h*(x(4) - T_0) ...
    -k_E*V_e_0^2 ...
    -k_Phi*phi_p_0^2 ...
    ) / C_th;

    -(x(5) - V_e_0)/tau_e
];

%% ------------------------------------------------------------
% Representative PID Plant
%% ------------------------------------------------------------

f_pid_rep = @(t,z) pid_dynamics( ...
    t,z, ...
    M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,V_e_0, ...
    Kp,Ki,Kd,u_min,u_max, ...
    tspan,d_profile_rep);

%% ------------------------------------------------------------
% Simulate Representative Trajectories
%% ------------------------------------------------------------

[t_nl,x_nl] = ode45( ...
    f_nl_rep,tspan,x_0_nl);

[t_pid,z_pid] = ode45( ...
    f_pid_rep,tspan,z_0_pid);

%% ============================================================
% Representative q Response
%% ============================================================

figure;

plot( ...
    t_nl,x_nl(:,1), ...
    'r','LineWidth',2);

hold on;

plot( ...
    t_pid,z_pid(:,1), ...
    'g','LineWidth',2);

yline(1e-2,'k--','LineWidth',1);
yline(-1e-2,'k--','LineWidth',1);

grid on;

xlabel('$t$ (s)');
ylabel('$q(t)$');

legend( ...
    'Nonlinear Uncontrolled', ...
    'PID Controlled', ...
    '$\pm 10^{-2}$ tolerance', ...
    'Location','best');

set(gca,'FontSize',12);

%% ============================================================
% Representative PID Regulation Error
%% ============================================================

figure;

semilogy( ...
    t_pid,abs(z_pid(:,1)), ...
    'g','LineWidth',2);

hold on;

yline( ...
    1e-2,'r--','LineWidth',0.5);

grid on;

xlabel('$t$ (s)');
ylabel('$|q(t)-q_{\mathrm{ref}}|$');

legend( ...
    'PID Controlled', ...
    '$10^{-2}$ tolerance', ...
    'Location','best');

%% ============================================================
% Representative Disturbance
% ============================================================

figure;

plot( ...
    tspan,d_profile_rep, ...
    'k','LineWidth',1.5);

grid on;

xlabel('$t$(s)');
ylabel('$d(t)$');

title('Representative External Disturbance');

%% ============================================================
% Local PID Dynamics Function
% ============================================================

function dz = pid_dynamics( ...
    t,z, ...
    M,C_v,K_v,alpha_B,alpha_E, ...
    gamma_Phi,k_Phi,k_E,C_th,h,tau_e, ...
    phi_p_0,T_0,V_e_0, ...
    Kp,Ki,Kd,u_min,u_max, ...
    tspan,d_profile)

q    = z(1);
qdot = z(2);
Phi  = z(3);
T    = z(4);
Ve   = z(5);
I_e  = z(6);

q_ref = 0;

e = q - q_ref;

dPhi = Phi - phi_p_0;

%% External disturbance

d = interp1( ...
    tspan,d_profile,t, ...
    'linear','extrap');

%% PID acceleration command

qddot_des = ...
    -Kp*e ...
    -Ki*I_e ...
    -Kd*qdot;

%% Linearized inverse electrostatic command

dVe_cmd = ...
    (M/(2*alpha_E*V_e_0)) * ( ...
        qddot_des ...
        + (C_v/M)*qdot ...
        + (K_v/M)*q ...
        - (2*alpha_B*phi_p_0/M)*dPhi);

u = V_e_0 + dVe_cmd;

%% Actuator saturation

u = min(max(u,u_min),u_max);

%% State derivatives

dz = zeros(6,1);

%% Mechanical dynamics

dz(1) = qdot;

dz(2) = ( ...
    -C_v*qdot ...
    -K_v*q ...
    +alpha_B*Phi^2 ...
    +alpha_E*Ve^2 ...
    -alpha_B*phi_p_0^2 ...
    -alpha_E*V_e_0^2 ...
    +d ...
    ) / M;

%% Pinned flux dynamics

dz(3) = ...
    -gamma_Phi*(Phi - phi_p_0);

%% Thermal dynamics

dz(4) = ( ...
    k_E*Ve^2 ...
    +k_Phi*Phi^2 ...
    -h*(T - T_0) ...
    -k_E*V_e_0^2 ...
    -k_Phi*phi_p_0^2 ...
    ) / C_th;

%% Electrical actuator dynamics

dz(5) = ...
    (u - Ve)/tau_e;

%% Integral error dynamics

dz(6) = e;

end

%% ============================================================
% Disturbance Generation Function
% ============================================================

function d_profile = generate_disturbance( ...
    t,dt,A1,f1,A2,f2,sigma_d,tau_d)

%% Zero-disturbance case

if A1 == 0 && ...
   A2 == 0 && ...
   sigma_d == 0

    d_profile = zeros(size(t));

    return;

end

%% Deterministic disturbance

d_det = ...
    A1*sin(2*pi*f1*t) ...
    + A2*sin(2*pi*f2*t);

%% Ornstein-Uhlenbeck colored stochastic disturbance

d_colored = zeros(size(t));

for k = 2:length(t)

    dW = sqrt(dt)*randn;

    d_colored(k) = ...
        d_colored(k-1) ...
        - (d_colored(k-1)/tau_d)*dt ...
        + sqrt(2*sigma_d^2/tau_d)*dW;

end

%% Total disturbance

d_profile = d_det + d_colored;

end




