clear; clc; close all;

%% Nominal Parameters

M = 1.0;
Cv = 0.3;
Kv = 2.25;

phi_p_0 = 1.0;
T_0 = 0.0;
Ve_0 = 1.0;

alpha_B = 0.20;
alpha_E = 0.12;

gamma_Phi = 0.15;
k_Phi = 0.04;
k_E = 0.03;

C_th = 1.0;
h = 0.25;

tau_e = 0.10;


%% Simulation Parameters

dt = 1e-3;
Tsim = 10;

t = 0:dt:Tsim;

N = 1000;        % Monte Carlo trials


x_ref = 1.0;    % desired stress coordinate


%% Storage

LQG_RMS = zeros(N,1);
SMC_RMS = zeros(N,1);

LQG_energy = zeros(N,1);
SMC_energy = zeros(N,1);

LQG_settle = zeros(N,1);
SMC_settle = zeros(N,1);

LQG_Vmax = zeros(N,1);
SMC_Vmax = zeros(N,1);

LQG_success = zeros(N,1);
SMC_success = zeros(N,1);

A_init = [
    0,        1,        0,                              0,              0;
   -Kv/M, -Cv/M, 2*alpha_B*phi_p_0/M,       0,              2*alpha_E*Ve_0/M;
    0,        0,       -gamma_Phi,                        0,              0;
    0,        0,        2*k_Phi*phi_p_0/C_th,          -h/C_th,        2*k_E*Ve_0/C_th;
    0,        0,        0,                              0,             -1/tau_e
    ];

B = [0;
0;
0;
0;
1/tau_e
];

    Cmat = [
1 0 0 0 0;
0 0 1 0 0;
0 0 0 1 0;
0 0 0 0 1
];

%% LQG Controller

Q = diag([100 10 5 5 1]);
R = 1;
K = lqr(A_init,B,Q,R);
W = 0.01*eye(5);
V = 0.01*eye(4);
L = lqe(A_init,eye(5),Cmat,W,V);


%% Monte Carlo Loop

for i = 1:N


    %% ----------------------------
    % Parameter uncertainty
    % -----------------------------

    M_i = M*(1+0.30*(2*rand-1));

    Cv_i = Cv*(1+0.50*(2*rand-1));

    Kv_i = Kv*(1+0.30*(2*rand-1));


    alpha_B_i = alpha_B*(1+0.20*(2*rand-1));

    alpha_E_i = alpha_E*(1+0.20*(2*rand-1));


    gamma_i = gamma_Phi*(1+0.30*(2*rand-1));

    kPhi_i = k_Phi*(1+0.30*(2*rand-1));

    kE_i = k_E*(1+0.30*(2*rand-1));

    h_i = h*(1+0.30*(2*rand-1));
    tau_i = tau_e*(1+0.25*(2*rand-1));


    %% Build uncertain plant

    A_i = [
    0,        1,        0,                              0,              0;
   -Kv_i/M_i, -Cv_i/M_i, 2*alpha_B_i*phi_p_0/M_i,       0,              2*alpha_E_i*Ve_0/M_i;
    0,        0,       -gamma_i,                        0,              0;
    0,        0,        2*kPhi_i*phi_p_0/C_th,          -h/C_th,        2*kE_i*Ve_0/C_th;
    0,        0,        0,                              0,             -1/tau_i
    ];

    x=zeros(5,1);
xhat=zeros(5,1);

    X_lqg=zeros(length(t),1);
    U_lqg=zeros(length(t),1);



    for k=1:length(t)

        y=Cmat*x;


        u=-K*(xhat-[x_ref;0;0;0;0]);


        xdot=A_i*x+B*u;

        x=x+dt*xdot;


        xhat=xhat+dt*(A_i*xhat+B*u+L*(y-Cmat*xhat));


        X_lqg(k)=x(1);

        U_lqg(k)=u;

    end



    %% ======================================================
    % Sliding Mode Controller
    % =======================================================

    lambda = 5;
    eta = 8;


    x=zeros(5,1);


    X_smc=zeros(length(t),1);
    U_smc=zeros(length(t),1);


    for k=1:length(t)


        pos=x(1);
        vel=x(2);


        error=pos-x_ref;


        s=vel+lambda*error;


        % Equivalent control approximation

        ueq = ...
        Kv_i*pos + Cv_i*vel;


        % Switching term

        usw = eta*sign(s);


        u = -(ueq+usw);


        xdot=A_i*x+B*u;

        x=x+dt*xdot;


        X_smc(k)=x(1);

        U_smc(k)=u;


    end



    %% ======================================================
    % Metrics
    % =======================================================

    LQG_RMS(i)=rms(X_lqg-x_ref);

    SMC_RMS(i)=rms(X_smc-x_ref);


    LQG_energy(i)=trapz(U_lqg.^2);

    SMC_energy(i)=trapz(U_smc.^2);


    LQG_Vmax(i)=max(abs(U_lqg));

    SMC_Vmax(i)=max(abs(U_smc));

    % settling time

    idx=find(abs(X_lqg-x_ref)<0.02,1);

if ~isempty(idx)
    LQG_settle(i)=t(idx);
    LQG_success(i)=1;
else
    LQG_settle(i)=NaN;
    LQG_success(i)=0;
end


    idx=find(abs(X_smc-x_ref)<0.02,1);

if ~isempty(idx)
    SMC_settle(i)=t(idx);
    SMC_success(i)=1;
else
    SMC_settle(i)=NaN;
    SMC_success(i)=0;
end


end

LQG_success = mean(~isnan(LQG_settle));
SMC_success = mean(~isnan(SMC_settle));


LQG_vmax_mean = mean(LQG_Vmax);
LQG_vmax_std = std(LQG_Vmax);

SMC_vmax_mean = mean(SMC_Vmax);
SMC_vmax_std = std(SMC_Vmax);

LQG_success_rate = mean(LQG_success)*100;
SMC_success_rate = mean(SMC_success)*100;


LQG_energy_mean = mean(LQG_energy);
LQG_energy_std = std(LQG_energy);

SMC_energy_mean = mean(SMC_energy);
SMC_energy_std = std(SMC_energy);


LQG_settle_mean = mean(LQG_settle,'omitnan');
SMC_settle_mean = mean(SMC_settle,'omitnan');

[max(SMC_Vmax), find(SMC_Vmax==max(SMC_Vmax))]

%% ==========================================================
% Results
% ===========================================================


controller = categorical([
    repmat("LQG",N,1);
    repmat("SMC",N,1)
]);



%% ---------------------------
% RMS Error
% ----------------------------
figure;

boxchart(controller,[LQG_RMS;SMC_RMS])

ylabel('RMS Error')

title('(a) Regulation Error')

grid off

set(gca,'FontSize',9)



%% ---------------------------
% Control Energy
% ----------------------------
figure;

boxchart(controller,[LQG_energy;SMC_energy])

set(gca,'YScale','log')

ylabel('Energy')

title('(b) Control Effort')

grid off

set(gca,'FontSize',9)



%% ---------------------------
% Settling Time
% ----------------------------

figure;

boxchart(controller,[LQG_settle;SMC_settle])

ylabel('Time (s)')

title('(c) Settling Time')

ylim([0 Tsim])

grid off

set(gca,'FontSize',9)



%% Figure size for IEEE column

set(gcf,'Units','inches')
set(gcf,'Position',[1 1 7.16 4.5])