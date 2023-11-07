close all; 
clear;
clc;

%% User Defines

l = 0.01;
a = 0.04;
k = 100000;
b = 0;
J = 0.1;
Fz = 50;
mu = 0.6;
mu_slip = 0.4;
n = 1000;
epsilon = 4e-10;
kp = 25;
kd = 4;
ki = 0;
psidesprofile = [0 0;
    0.01, -0.7;
    3, -0.7;
    3.01, 0.6;
    6, 0.6;
    6.01, -0.5
    9, -0.5;
    9.01, 0.4
    12, 0.4
    ];

%% Simulation parameters

SimParms = struct;
SimParms.l = l;
SimParms.k = k;
SimParms.J = J;
SimParms.n = n;
SimParms.epsilon = epsilon;
SimParms.kp = kp;
SimParms.kd = kd;
SimParms.ki = ki;
SimParms.psidesprofile = psidesprofile;

SimParms.dx = 2*a/n;

x = linspace(-a,a,n+1);
SimParms.x = x;

qcrp = 3*Fz*mu/(4*a^3*k)*(a^2-x.^2);
q_slip = 3*Fz*mu_slip/(4*a^3*k)*(a^2-x.^2);
SimParms.qcrp = qcrp;
SimParms.q_slip = q_slip;

%% Time discretization

% kezdeti ertekek / initial conditions
y0=[0;0;((1:n+1)==0)';0];

% lepeskoz / timestep size
dt=0.01;%dx/v/2

% kezdeti idopillanat / initial time value
t0=0;

% vegso idopillanat / final time value
tfinal = psidesprofile(end,1);
t = t0:dt:tfinal;
pont = length(t);

%% Simulation

y = zeros(pont,n+4);
y(1,:) = y0;
tic;
ii = 1;
wb = waitbar(0,'Simulation in process...');

for i=2:pont
    [y(i,:), y(i-1,:)] = solver(t(i),y(i-1,:),dt,SimParms);
    if i > (fix(pont/10)*ii + 1)
        ii = ii + 1;
    end
    waitbar(i/pont,wb)
end
simtime = toc
close(wb);
q = y(:,3:n+3);

%% Animation

if true
for i=1:pont
        figure(3)
        plot(x(:),y(i,3:n+3),'.-')
        hold on;
        plot(x, qcrp, 'k--');
        plot(x, q_slip, 'k:');
        plot(x, -qcrp, 'k--');
        plot(x, -q_slip, 'k:');
        hold off;
        grid on;
        ylim([-0.01,0.01]);
        drawnow;
end
end

%% Plots

% elfordulas az idoben
figure;
plot(t,y(:,1),'r');
hold on;
plot(psidesprofile(:,1), psidesprofile(:,2), 'k');
hold off;
xlabel('t [sec]');
ylabel('\psi [rad]');
legend('Simulated', 'Desired');
grid on;

% figure;
% plot(t,y(:,end),'-');
% grid on;

%% plots for TDK

if true
    
    Psi = y(:,1);
    Psides = zeros(size(Psi));
    Omega = y(:,2);
    Z = y(:,end);
    
    for i = 1:length(Psides)
        psides = [];
        t_i = t(i);
        for ii = 1:length(psidesprofile(:,1))
            if t_i <= psidesprofile(ii,1)
                psides = psidesprofile(ii,2);
                break;
            end
        end

        if isempty(psides)
            psides = psidesprofile(end,2);
        end
        Psides(i) = psides;
    end
    
    M = -kp*(Psi - Psides)-kd*Omega-ki*Z;
    
    figure;
    subplot(2,1,1);
    plot(t, M, 'k', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Steering torque [Nm]');
    grid on;
    ylim([-15, 15]);
    
    subplot(2,1,2);
    plot(t, Psi, 'k', 'LineWidth', 1);
    hold on;
    plot([2, 5, 8, 11], [Psi(201), Psi(501), Psi(801), Psi(1101)], 'ok', 'MarkerFaceColor', 'k');
    hold off;
    xlabel('Time [s]');
    ylabel('Yaw angle [rad]');
    grid on;
    ylim([-0.75, 0.75])
    
    ymax = max(qcrp);
    
    figure;
    subplot(2, 2, 1);
    plot(x(:),y(201,3:n+3),'k', 'LineWidth', 1.2)
    hold on;
    plot(x, qcrp, 'k--');
    plot(x, q_slip, 'k:');
    plot(x, -qcrp, 'k--');
    plot(x, -q_slip, 'k:');
    hold off;
    grid on;
    ylim([-1.1*ymax, 1.1*ymax]);
    xlabel('x [m]');
    ylabel('q(x,t) [m]');
    title('$t_1=2$ s', 'Interpreter', 'Latex');
    
    subplot(2, 2, 2);
    plot(x(:),y(501,3:n+3),'k', 'LineWidth', 1.2)
    hold on;
    plot(x, qcrp, 'k--');
    plot(x, q_slip, 'k:');
    plot(x, -qcrp, 'k--');
    plot(x, -q_slip, 'k:');
    hold off;
    grid on;
    ylim([-1.1*ymax, 1.1*ymax]);
    xlabel('x [m]');
    ylabel('q(x,t) [m]');
    title('$t_2=5$ s', 'Interpreter', 'Latex');
    
    subplot(2, 2, 3);
    plot(x(:),y(801,3:n+3),'k', 'LineWidth', 1.2)
    hold on;
    plot(x, qcrp, 'k--');
    plot(x, q_slip, 'k:');
    plot(x, -qcrp, 'k--');
    plot(x, -q_slip, 'k:');
    hold off;
    grid on;
    ylim([-1.1*ymax, 1.1*ymax]);
    xlabel('x [m]');
    ylabel('q(x,t) [m]');
    title('$t_3=8$ s', 'Interpreter', 'Latex');
   
    subplot(2, 2, 4);
    plot(x(:),y(1101,3:n+3),'k', 'LineWidth', 1.2)
    hold on;
    plot(x, qcrp, 'k--');
    plot(x, q_slip, 'k:');
    plot(x, -qcrp, 'k--');
    plot(x, -q_slip, 'k:');
    hold off;
    grid on;
    ylim([-1.1*ymax, 1.1*ymax]);
    xlabel('x [m]');
    ylabel('q(x,t) [m]');
    title('$t_4=11$ s', 'Interpreter', 'Latex');
end
