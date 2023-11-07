close all; 
clear;
clc;

anim = 0;
fric_calc = 1;

%% User Defines
addpath('Old_meas');
load fem6

l = e;
b = 6.8509;
J = 0.145228 + 2.513*(0.042 + e)^2;
mu = 1.674; 
mu_slip = 1.467;
n = 1000;
epsilon = 4e-10;

%% Simulation parameters

SimParms = struct;
SimParms.l = l;
SimParms.k = k;
SimParms.J = J;
SimParms.n = n;
SimParms.epsilon = epsilon;
SimParms.Psi = Psi;
SimParms.Time = time;
SimParms.Fz = Fz;
SimParms.Omega = [0, diff(Psi)./diff(time)];
SimParms.mu = mu;
SimParms.mu_slip = mu_slip;
SimParms.a = a;

SimParms.dx = 2*a/n;

x = linspace(-a,a,n+1);
SimParms.x = x;

SimParms.Mz = Mz;

SimParms.OL = 1;

%% Time discretization

% kezdeti ertekek / initial conditions
y0=[0;((1:n+1)==0)'];

% lepeskoz / timestep size
dt=0.01;%dx/v/2

% kezdeti idopillanat / initial time value
t0=0;

% vegso idopillanat / final time value
tfinal = time(end);
t = t0:dt:tfinal;
pont = length(t);

%% Simulation

y = zeros(pont,n+2);
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
toc;
close(wb);
q = y(:,2:end);

%% Animation

if anim
for i=1:pont                                                                
        figure(3)
        plot(x(:),y(i,2:end),'.-');
        hold on;
        plot(x, 3*Fz(i)*mu/(4*a^3*k)*(a^2-x.^2), 'k--');
        plot(x, -3*Fz(i)*mu/(4*a^3*k)*(a^2-x.^2), 'k--');
        plot(x, 3*Fz(i)*mu_slip/(4*a^3*k)*(a^2-x.^2), 'k:');
        plot(x, -3*Fz(i)*mu_slip/(4*a^3*k)*(a^2-x.^2), 'k:');
        hold off;
        grid on;
        ylim([-0.01,0.01]);
        drawnow;
end
end

%% Plots

SIM_RMSE = sqrt(mean((Mz' - y(:,1)).^2))

% elfordulas az idoben
figure;
plot(t,y(:,1),'r');
hold on;
plot(time, Mz, 'k');
hold off;
xlabel('Time [s]');
ylabel('Self-aliging torque [Nm]');
grid on;
xlim([0, time(end)]);
legend('Sim.', 'Meas.');

psi = interp1(time, Psi, t);
figure;
plot(psi,y(:,1),'r');
hold on;
plot(Psi, Mz, 'k');
grid on;
hold off;
xlabel('Yaw angle [rad]');
ylabel('Self-aliging torque [Nm]');
grid on;
legend('Sim.', 'Meas.');

%% Friction calculation
if fric_calc
mu_approx = friction_calc(SimParms);

figure;
plot(time, mu.*ones(size(time)), 'k', 'LineWidth', 1.2);
hold on;
plot(time, 1.1*mu.*ones(size(time)), '--k');
plot(time, 0.9*mu.*ones(size(time)), '--k');
plot(time, movmean(mu_approx, [15,0]), 'r', 'LineWidth', 1.2);
xlim([0, time(end)]);
hold off;
legend('Meas.', '10% error', '', 'Estim.', 'Location', 'best');
xlabel('Time [s]');
ylabel('Coeff. of friction [1]');
grid on;

RMSE = sqrt(mean((movmean(mu_approx, [15,0]) - mu.*ones(size(time))).^2))
end


