set(groot,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

num_rad = 500;
num_theta = 500;

r = linspace(0,1,num_rad);
theta = linspace(0,2*pi,num_theta);

[Rad,Theta] = meshgrid(r,theta);
x_polar = Rad.*cos(Theta);
y_polar = Rad.*sin(Theta);

%% Number of Snapshots
n_snaps = 500;
t = linspace(0,10,n_snaps);

%% Maxwell Pressure Field Snapshots
% P = f(r, theta, t) is the Maxwell EM pressure tensor
% These values are normal to the controlled surface

%% Normalized Stress Amplitudes
% P_0 is set to unity so all perturbations are fractional deviations
% asym = 0.20 models the dominant dipole/asymmetry mode
% qd = 0.12 models a weaker quadrupole distortion
% arm = 0.06 models a smaller axial/radial gap-like perturbation
% Note that this hierarchy of asym > qd > arm reflects ASSUMED modal energy ordering

P_0 = 1;            % Baseline Pinned Flux Pressure
asym = 0.20;        % Dipole Asymmetry
qd = 0.12;          % Quadrupole distortion
arm = 0.06;         % Axial Radial Mode
noiseAmp = 0.01;

x = zeros(num_rad*num_theta, n_snaps);

for k = 1:n_snaps

    dipoleMode = asym*sin(0.8*t(k))*Rad.*cos(Theta);
    quadMode = qd*cos(1.1*t(k))*Rad.^2.*cos(2*Theta);
    radialMode = arm*sin(0.4*t(k))*(1 - Rad).^2;
    thermalDrift = 0.04*(1-exp(-0.2*t(k)));

    P = P_0 + dipoleMode + quadMode + radialMode + thermalDrift + noiseAmp*randn(size(Rad));

    x(:,k) = P(:);

end

%% Center the Snapshot Matrix
x_mean = mean(x,2);
x_c = x - x_mean;

%% POD via SVD
[U,S,V] = svd(x_c,'econ');
singleVals = diag(S);
modalEnergy = singleVals.^2 / sum(singleVals.^2);
cumuEnergy = cumsum(modalEnergy);

%% Choose Modes
energyThreshold = 0.971; % Play with this
n_keeps = find(cumuEnergy >= energyThreshold, 1, 'first');

snapID = round(0.35*n_snaps);
P_snap = reshape(x(:,snapID),num_theta,num_rad);

fprintf('Number of Modes required for %.1f%% Cumulative Energy: %d\n',100*energyThreshold,n_keeps);

U_r = U(:,1:n_keeps);

% Nominal reduced modal coordinates
q_nom = U_r'*x_c;

% Estimate nominal modal acceleration
dt = mean(diff(t));
qdot_nom  = gradient(q_nom,dt);
qddot_nom = gradient(qdot_nom,dt);

%% Monte Carlo Simulation Setup
num_sims = 1;

% Uncertainty ranges for MC param sweep
P_0_range = P_0*[0.75, 1.33];
asym_range = asym*[0.75, 1.33];
qd_range = qd*[0.75, 1.33];
arm_range = arm*[0.75, 1.33];
noise_range = noiseAmp*[0.50, 2.00];
thermal_range = 0.04*[0.50, 2];

delta_max = 0;
delta_all = zeros(num_sims,1);

for mc = 1:num_sims

    % Randomly sample uncertain parameters
    P0_mc = randUniform(P_0_range);
    asym_mc = randUniform(asym_range);
    qd_mc = randUniform(qd_range);
    arm_mc = randUniform(arm_range);
    noise_mc = randUniform(noise_range);
    therm_mc = randUniform(thermal_range);

    x_mc = zeros(num_rad*num_theta,n_snaps);

    for k = 1:n_snaps

        dipoleMode = asym_mc*sin(0.8*t(k))*Rad.*cos(Theta);
        quadMode   = qd_mc*cos(1.1*t(k))*Rad.^2.*cos(2*Theta);
        radialMode = arm_mc*sin(0.4*t(k))*(1 - Rad).^2;
        thermalDrift = therm_mc*(1-exp(-0.2*t(k)));

        P = P0_mc + dipoleMode + quadMode + radialMode + thermalDrift + noise_mc*randn(size(Rad));

        x_mc(:,k) = P(:);

    end

    % Project perturbed snapshots onto nominal retained POD basis
    q_mc = U_r'*(x_mc - x_mean);

    % Modal acceleration under uncertainty
    qdot_mc  = gradient(q_mc,dt);
    qddot_mc = gradient(qdot_mc,dt);

    % Lumped uncertainty estimate in retained modal dynamics
    delta_mc = qddot_mc - qddot_nom;

    % Worst-Case Modal Acceleration Mismatch
    delta_all(mc) = max(abs(delta_mc),[],'all');
    delta_max = max(delta_max,delta_all(mc));

end

fprintf('Estimated Lumped Uncertainty Bound Delta_max = %.4e\n',delta_max);

% *Probably* Safe SMC Switching Gain
safetyFactor = 1.25;
k_ss = safetyFactor*delta_max;

fprintf('Safe SMC Switching Gain k_s > %.4e\n',k_ss);


% %% Figure 1: Example Maxwell Pressure Heatmap
% figure;
% imagesc(r,theta,P_snap);
% axis xy;
% xlabel('Normalized Radius');
% ylabel('$\theta$  (radians)');
% title('Normal Pressure, $P(r, \theta, t)$');
% colorbar;
% 
% %% Figure 2: Singular value decay
% figure;
% semilogy(singleVals,'o-');
% xlabel('Mode Index');
% ylabel('Singular Value');
% title('Singular Value Decay');
% grid on;
% 
% %% Figure 3: Cumulative Modal Energy
% figure;
% plot(cumuEnergy,'o-');
% hold on;
% yline(0.95,'r--');
% xline(n_keeps,'k--');
% text(50,0.945,'95% Energy','FontSize',10);
% text(6,0.83,sprintf('N = %d',n_keeps),'FontSize',10);
% xlabel('Number of Modes');
% ylabel('Percentage of Cumulative Energy');
% title('Cumulative Energy VS Number of Modes Kept');
% grid on;
% xlim([-10 210]);
% ylim([0.75 1.02]);
% 
% %% Figure 4: First 3 POD Mode Shapes
singVals = diag(S).^2;
energyFrac = singVals/sum(singVals);

num_PlotModes = min(3,size(U,2));

modeLabels = {...
    sprintf('Mode 1 (Retained for Control): %.1f%% Energy',100*energyFrac(1)),...
    sprintf('Mode 2: %.1f%% Energy',100*energyFrac(2)),...
    sprintf('Mode 3: %.1f%% Energy',100*energyFrac(3))};

%% Global Color Limits

cmin = inf;
cmax = -inf;

for i = 1:num_PlotModes
    modeShape = reshape(U(:,i),num_theta,num_rad);
    cmin = min(cmin,min(modeShape(:)));
    cmax = max(cmax,max(modeShape(:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1 - Mode 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w',...
       'Units','normalized',...
       'Position',[0.30 0.10 0.35 0.70]);

tl = tiledlayout(1,1,...
    'TileSpacing','compact',...
    'Padding','compact');

modeShape = reshape(U(:,1),num_theta,num_rad);

ax = nexttile;

surf(ax,...
     x_polar,...
     y_polar,...
     modeShape,...
     'EdgeColor','none');

view(ax,2)

axis(ax,'equal')
axis(ax,'tight')

clim(ax,[cmin cmax])

xlabel(ax,'$x/R$',...
    'Interpreter','latex',...
    'FontName','Times New Roman',...
    'FontSize',12)

ylabel(ax,'$y/R$',...
    'Interpreter','latex',...
    'FontName','Times New Roman',...
    'FontSize',12)

set(ax,...
    'FontName','Times New Roman',...
    'FontSize',12,...
    'LineWidth',1.0,...
    'Box','on');

title(ax,modeLabels{1},...
    'FontName','Times New Roman',...
    'FontWeight','normal',...
    'FontSize',14);

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontName = 'Times New Roman';
cb.FontSize = 12;
cb.Label.String = 'Mode Amplitude';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 10;
exportgraphics(gcf,'dom_pod_mode.png','Resolution',600);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 - Modes 2 & 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w',...
       'Units','normalized',...
       'Position',[0.30 0.05 0.35 0.70]);

tl = tiledlayout(2,1,...
    'TileSpacing','loose',...
    'Padding','compact');

for i = 2:3

    modeShape = reshape(U(:,i),num_theta,num_rad);

    ax = nexttile;

    surf(ax,...
         x_polar,...
         y_polar,...
         modeShape,...
         'EdgeColor','none');

    view(ax,2)

    axis(ax,'equal')
    axis(ax,'tight')

    clim(ax,[cmin cmax])

xlabel(ax,'$x/R$',...
    'Interpreter','latex',...
    'FontName','Times New Roman',...
    'FontSize',12)

ylabel(ax,'$y/R$',...
    'Interpreter','latex',...
    'FontName','Times New Roman',...
    'FontSize',12)

set(ax,...
    'FontName','Times New Roman',...
    'FontSize',12,...
    'LineWidth',1.0,...
    'Box','on');

title(ax,modeLabels{i},...
    'FontName','Times New Roman',...
    'FontWeight','normal',...
    'FontSize',14);

end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontName = 'Times New Roman';
cb.FontSize = 12;
cb.Label.String = 'Mode Amplitude';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 12;
%% Export High Resolution

exportgraphics(gcf,'pod_modes_horz.png','Resolution',600);

% 
% %% Figure 5: Modal Coefficient Time Histories
% figure;
% hold on;
% 
% for i = 1:num_PlotModes
%     plot(t,modalCoeffs(i,:));
% end
% 
% xlabel('Time [sec]');
% ylabel('Modal Coefficient');
% title('Dominant Modal Coordinates');
% legend(arrayfun(@(i)sprintf('Mode %d',i),1:num_PlotModes,'UniformOutput',false));
% grid on;

%% Figure 6: Uncertainty Distribution

figure;
histogram(delta_all,30);
xlabel('Maximum Modal Acceleration Mismatch');
ylabel('Monte Carlo Count');
title('Monte Carlo Estimate of Lumped Uncertainty Bound');
grid on;

%% Figure 7: Params Sweep
figure;
plot(delta_all,'o');
hold on;
yline(delta_max,'r--');
yline(k_ss,'k--');
xlabel('Trial Number');
ylabel('Uncertainty Bound Estimate');
title('Monte Carlo Parameter Sweep for SMC Gain Selection');
legend('$\Delta$  Estimate','$\Delta_{max}$','Suggested  $k_s$','Location','best');
grid on;

%% Figure 8: Reduced-Order Modeling and Control Workflow

figure('Color','w','Position',[200 200 1100 420]);
axis off;

boxW = 0.13;
boxH = 0.16;
y = 0.55;

xpos = [0.04 0.19 0.34 0.49 0.64 0.79];

labels = {
    'Pinned Flux Field'
    'Maxwell Stress Tensor'
    'Stress Distribution'
    'Snapshot Matrix'
    'POD / Modal Truncation'
    'LQR / SMC Control'
};

subtext = {
    'B_p,  \Phi'
    'T_{EM}(E, B)'
    'P(r, \theta, t)'
    'X = [x_1 ... x_N]'
    'U_r ,  q(t)'
    'u(t) ,  V_e(t)'
};

for i = 1:length(xpos)
    annotation('rectangle',[xpos(i) y boxW boxH]);

    annotation('textbox',[xpos(i) y+0.055 boxW boxH*0.45], ...
        'String',labels{i}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'EdgeColor','none', ...
        'FontSize',10, ...
        'FontWeight','bold');

    annotation('textbox',[xpos(i) y+0.005 boxW boxH*0.35], ...
        'String',subtext{i}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'EdgeColor','none', ...
        'FontSize',10);
end

for i = 1:length(xpos)-1
    x1 = xpos(i)+boxW;
    x2 = xpos(i+1);
    annotation('arrow',[x1 x2],[y+boxH/2 y+boxH/2]);
end

% Feedback Loop from Controller to Actuation
annotation('line',[0.86 0.86],[0.54 0.30]);
annotation('line',[0.86 0.10],[0.30 0.30]);
annotation('arrow',[0.10 0.10],[0.30 0.54]);

annotation('textbox',[0.32 0.18 0.35 0.08], ...
    'String','Closed-Loop Electrostatic Modulation and Stress Regulation', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none', ...
    'FontSize',10, ...
    'FontAngle','italic');
% 
% % Optional export
% exportgraphics(gcf,'workflow_reduced_order_control.png','Resolution',300);
% exportgraphics(gcf,'workflow_reduced_order_control.pdf','ContentType','vector');

%% Generate a Random Value in Range
function vals = randUniform(range)
    vals = range(1) + (range(2)-range(1))*rand;
end