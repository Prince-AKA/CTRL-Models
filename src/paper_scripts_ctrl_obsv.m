%% Controllability and Observability Analysis

%% State Space Parameter Definitions
M = 1.0;
C = 0.4;
K = 2.0;

alpha_B = 0.8;
alpha_E = 1.2;

gamma_phi = 0.15;

C_th = 1.0;
k_E  = 0.05;
k_Phi = 0.03;
h_th = 0.20;

tau_e = 0.05;

% Operating point
q_0 = 0;
q_dot_0 = 0;
Phi_op = 1.0;
Phi_q = 1.0;    % normalized flux quantum/reference
T0 = 77;
V_e_0 = 0.5;

%% Linearized A matrix
A = [0, 1, 0, 0, 0;
-K/M, -C/M, 2*alpha_B*Phi_op/M, 0, 2*alpha_E*V_e_0/M;
0, 0, -gamma_phi*(2*Phi_op - Phi_q), 0, 0;
0, 0, 2*k_Phi*Phi_op/C_th, -h_th/C_th, 2*k_E*V_e_0/C_th;
0, 0, 0, 0, -1/tau_e];

%% Input matrix
B = [0; 0; 0; 0; 1/tau_e];

%% Output matrix
C = [1 0 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];

D = zeros(size(C,1),1);

%% Controllability / observability
P_o = ctrb(A,B);
P_c = obsv(A,C);

fprintf('Controllability rank = %d out of %d\n',rank(P_o),size(A,1));
fprintf('Observability rank   = %d out of %d\n',rank(P_c),size(A,1));