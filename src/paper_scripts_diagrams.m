%% Superconductor for 2.1
figure('Color','w','Position',[200 120 900 650]);
hold on; 
R = 1.0;
th = linspace(0,2*pi,400);

% Superconducting domain
fill(R*cos(th),R*sin(th),[0.92 0.94 1.00], ...
    'EdgeColor','k','LineWidth',1.5);

% Main label moved upward away from vortices
text(0,0.52,'Type-II Superconductor', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',14);

text(0,0.25,'Pinned Abrikosov Vortex Ensemble', ...
    'HorizontalAlignment','center', ...
    'FontAngle','italic', ...
    'FontSize',11);
%% Vortex cores in loose lattice, avoiding center text
[xgrid,ygrid] = meshgrid([-0.55 -0.25 0.05 0.35 0.65],[-0.55 -0.25 0.05]);
xv = xgrid(:);
yv = ygrid(:);

% remove points too close to central text region
keep = ~(abs(xv)<0.45 & yv > 0.10);
xv = xv(keep);
yv = yv(keep);

% keep only inside disk
keep = (xv.^2 + yv.^2) < 0.72^2;
xv = xv(keep);
yv = yv(keep);

for k = 1:length(xv)
    % Vortex Core
    plot(xv(k),yv(k),'ko','MarkerFaceColor','k','MarkerSize',8);

    % Small Circular Supercurrent Loop
    plot(xv(k)+0.060*cos(th), yv(k)+0.060*sin(th),'k-','LineWidth',1);

    % Tiny tangential Arrow for Circulation
    ang = pi/4;
    quiver(xv(k)+0.03*cos(ang), ...
           yv(k)+0.03*sin(ang), ...
          -0.040*sin(ang), ...
           0.040*cos(ang), ...
           0, ...
           'Color','r', ...
           'LineWidth',0.7, ...
           'MaxHeadSize',2);
end

%% Vortex label outside disk with leader arrow
text(-1.6,0.45,'Pinned Vortex Cores', ...
    'FontSize',10, ...
    'HorizontalAlignment','left');

plot([-1.25 -0.25],[0.4 -0.26],'k-','LineWidth',0.8);
plot([-1.25 0.05],[0.4 -0.24],'k-','LineWidth',0.8);

%% Applied magnetic bias field
quiver(0,1.45,0,-0.32,0, ...
    'LineWidth',2, ...
    'MaxHeadSize',0.8);

text(0,1.52,'Applied Magnetic Field  $B_p$', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'FontSize',12);

%% Segmented electrostatic gates
gateY = -1.45;
gateW = 0.32;
gateH = 0.16;
gateX = [-0.75 -0.25 0.25 0.75];

for k = 1:4
    rectangle('Position',[gateX(k)-gateW/2 gateY gateW gateH], ...
        'FaceColor',[0.86 0.86 0.86], ...
        'EdgeColor','k', ...
        'LineWidth',1.1);

    text(gateX(k),gateY+gateH/2,sprintf('$V_%d$',k), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',12);
end

text(0,gateY-0.16,'Segmented Electrostatic Gate Array', ...
    'HorizontalAlignment','center', ...
    'FontSize',11, ...
    'FontWeight','bold');

%% Multiple Electrostatic Modulation Arrows from Gates
for k = 1:4
    quiver(gateX(k),gateY+gateH+0.05,0,0.34,0,'Color','k',...
        'LineWidth',1.5, ...
        'MaxHeadSize',0.65);
end

text(0,-1.08,'Localized Electrostatic Modulation', ...
    'HorizontalAlignment','center', ...
    'FontSize',10);

%% Schematic field-shaping arcs above gates
arcR = 0.55;
arcTheta = linspace(210*pi/180,330*pi/180,80);
plot(arcR*cos(arcTheta), -0.72 + 0.18*sin(arcTheta), ...
    'k--','LineWidth',0.9);

text(0,-0.78,'*Stress Field Shaping Region*', ...
    'HorizontalAlignment','center', ...
    'FontSize',9);

%% Maxwell stress output
quiver(1.08,0,0.36,0,0, ...
    'LineWidth',2, ...
    'MaxHeadSize',0.7);

text(1.7,0.09,'Electromagnetic Stress', ...
    'HorizontalAlignment','center', ...
    'FontSize',10);

text(1.7,0,'$\mathbf{T}_{EM}(\mathbf{E},B_p)$', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'FontSize',12);

%% Hall Sensor
rectangle('Position',[1.34,-0.62,0.48,0.23], ...
    'FaceColor',[0.96 0.96 0.96], ...
    'EdgeColor','k', ...
    'LineWidth',1.1);

text(1.58,-0.505,'Hall Sensor', ...
    'HorizontalAlignment','center', ...
    'FontSize',10);

quiver(1.17,-0.35,0.18,-0.08,0, ...
    'LineWidth',1.2, ...
    'MaxHeadSize',0.8);

text(1.6,-0.74,'Field Measurement', ...
    'HorizontalAlignment','center', ...
    'FontSize',9);

%% Controlled surface callout
text(-1.55,-0.85,'Controlled Surface', ...
    'FontSize',10, ...
    'HorizontalAlignment','left');

plot([-0.995 -0.47],[-0.83 -0.68],'k-','LineWidth',0.8);

xlim([-1.85 2.05]);
ylim([-1.70 1.70]);

axis equal;
axis off;

set(gca,'Clipping','off');


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 2.4: Reduced-Order Pinned-Flux Damped Mode
%clear; close all; clc;

figure('Color','w','Position',[100 100 850 500]);
axis equal; axis off; hold on;

rectangle('Position',[1 1 6 3.2],'Curvature',[0.15 0.15], ...
    'FaceColor',[0.88 0.92 1.0],'EdgeColor','k','LineWidth',1.5);

text(4,3.95,'Type-II Superconductor / Pinned Vortex Ensemble', ...
    'HorizontalAlignment','center','FontSize',12,'Interpreter','latex');

x0 = [2.0 3.0 4.0 5.0 6.0 2.5 3.5 4.5 5.5];
y0 = [2.7 2.7 2.7 2.7 2.7 1.8 1.8 1.8 1.8];

scale = 0.35;
dx = scale * sin((x0-1)/6*pi);
dy = 0.15 * cos((x0-1)/6*pi);

for i = 1:length(x0)
    plot(x0(i),y0(i),'ko','MarkerSize',12,'LineWidth',1.2);
    plot(x0(i),y0(i),'k.','MarkerSize',18);

    plot(x0(i)+dx(i),y0(i)+dy(i),'ro','MarkerSize',8,'LineWidth',1.2);

    quiver(x0(i),y0(i),dx(i),dy(i),0,'Color',[0.6 0 0], ...
        'LineWidth',1.2,'MaxHeadSize',0.8);
end

text(4,1.25,'Dominant modal displacement: $q(t)$', ...
    'HorizontalAlignment','center','FontSize',12,'Interpreter','latex');

quiver(0.9,2.6,0.8,0,0,'LineWidth',1.8,'Color',[0.0 0.3 0.8], ...
    'MaxHeadSize',0.8);
text(0.65,2.85,'$F_{\mathrm{EM}}$', ...
    'FontSize',12,'Color',[0.0 0.3 0.8],'Interpreter','latex');

quiver(7.25,2.6,-0.8,0,0,'LineWidth',1.8,'Color',[0.8 0.2 0.0], ...
    'MaxHeadSize',0.8);
text(7.45,2.85,'$K_v q$', ...
    'FontSize',12,'Color',[0.8 0.2 0.0],'Interpreter','latex');

quiver(7.25,1.8,-0.8,0,0,'LineWidth',1.8,'Color',[0.3 0.3 0.3], ...
    'MaxHeadSize',0.8);
text(7.45,2.05,'$C_v \dot{q}$', ...
    'FontSize',12,'Color',[0.3 0.3 0.3],'Interpreter','latex');

text(4,0.45, ...
    '$M_v\ddot{q} + C_v\dot{q} + K_vq = F_{\mathrm{EM}} + d$', ...
    'HorizontalAlignment','center','FontSize',14,'Interpreter','latex');

% %% Section 3.3: Microscopic Origin of Damped Elastic Mode
% clear; close all; clc;
% 
% set(groot,'defaultTextInterpreter','latex');
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure('Color','w','Position',[100 100 850 520]);
% axis equal; axis off; hold on;
% 
% % Pinning potential well
% x = linspace(-1.8,1.8,400);
% U = 0.5*x.^2 + 0.08*x.^4;
% plot(x+4,U+0.85,'k','LineWidth',0.1);
% 
% text(4,2.6,'Pinning Potential Well $U_p(r)$', ...
%     'HorizontalAlignment','center','FontSize',12.5);
% 
% % Equilibrium and displaced positions
% xeq = 4.0;
% yeq = 0.85;
% 
% xd = 4.65;
% yd = 0.5*(xd-4).^2 + 0.08*(xd-4).^4+0.85;
% 
% % Neighboring vortices, smaller and deliberately labeled
% neighbors = [2.8 1.35; 5.36 1.33; 3.25 2; 4.08 2.31; 4.94 2.15];
% 
% for i = 1:size(neighbors,1)
%     plot(neighbors(i,1),neighbors(i,2),'ko','MarkerSize',8,'LineWidth',1.0);
%     plot(neighbors(i,1),neighbors(i,2),'k.','MarkerSize',12);
%     plot([neighbors(i,1) xd],[neighbors(i,2) yd], ...
%         'k--','LineWidth',0.4);
% end
% 
% text(3.76,2.13,'Neighboring Vortices', ...
%     'FontSize',10,'HorizontalAlignment','center');
% 
% % Equilibrium vortex
% plot(xeq,yeq,'ko','MarkerSize',13,'LineWidth',1.3);
% plot(xeq,yeq,'k.','MarkerSize',20);
% text(xeq,yeq-0.08,'Pinned Equilibrium', ...
%     'HorizontalAlignment','center','FontSize',10);
% 
% % Displaced vortex
% plot(xd,yd,'ro','MarkerSize',13,'LineWidth',1.4);
% plot(xd,yd,'r.','MarkerSize',20);
% text(4.7,1,'Displaced Vortex', ...
%     'FontSize',10,'Color',[1 0 0]);
% 
% % Displacement vector
% quiver(xeq,yeq-0.015,xd-xeq,yd-yeq-0.015,0, 'r--','LineWidth',0.4, ...
%     'MaxHeadSize',1e-10);
% text(4.33,0.904,'$\Delta r \sim q$', ...
%     'FontSize',11,'Color',[1 0 0]);
% 
% % Restoring force vector, tangent toward well minimum
% quiver(xd,yd,-0.65,-0.22,0, ...
%     'Color',[0 0.25 0.85],'LineWidth',1.5, ...
%     'MaxHeadSize',0.27);
% text(xd-0.88,yd+0.01,'Restoring Force', ...
%     'FontSize',10,'Color',[0 0.25 0.85]);
% text(xd-0.85,yd-0.09,'$F_k=-K_vq$', ...
%     'FontSize',10,'Color',[0 0.25 0.85]);
% 
% % Damping force vector, opposing velocity/motion
% quiver(xd,yd,-0.35,0.42,0, ...
%     'Color',[0.3 0.3 0.3],'LineWidth',1.5, ...
%     'MaxHeadSize',0.28);
% text(xd-0.67,yd+0.575,'Damping Force', ...
%     'FontSize',10,'Color',[0.3 0.3 0.3]);
% text(xd-0.64,yd+0.475,'$F_c=-C_v\dot q$', ...
%     'FontSize',10,'Color',[0.3 0.3 0.3]);
% 
% % Compact explanatory line
% text(4,0.48,'$M_v\ddot q+C_v\dot q+K_vq=F_{\mathrm{EM}}+d$','HorizontalAlignment','center','FontSize',13);
% 
% 
% xlim([2.1 6.1]);
% ylim([0.0 3.75]);