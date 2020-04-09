clc
clear
close all

%% Dati

M1 = 1; M2 = 0.35;
J1 = 0.005;
R1 = 0.1; R2 = 0.3;
k1 = 18; k2 = 25;
c1 = 0.7; c2 = 1.2;
g = 9.81;

%% 1) Equazione del moto: equazione di Lagrange

% a)

syms F

m_g = J1 + M2 * R2^2;
c1_g = c1 * R2^2 + c1 * R1^2;
k_g = k1 * R2^2 + k2 * R1^2;
C_g = F * R2 + M2 * g * R2;

% b)

csi1 = c1_g / (2 * sqrt(m_g * k_g));

% c)

omega_n = sqrt(k_g / m_g);
alpha = csi1 * omega_n;
omega_d = sqrt(omega_n^2 - alpha^2);

%% 2) Moto libero

% a)

lambda1 = roots([m_g c1_g k_g]);

% Condizioni iniziali
theta0 = 2;
omega0 = 16;

t = 0:0.01:15;

X2 = (omega0 - lambda1(1) * theta0)/(lambda1(2) - lambda1(1));
X1 = theta0 - X2;

theta1 = X1*exp(lambda1(1)*t) + X2*exp(lambda1(2)*t);

% Grafico
figure
plot(t,theta1);
axis([-inf, inf, min(theta1)*1.1, max(theta1)*1.1]);
title('Time response of system free motion');
xlabel('time [s]'); ylabel('displacement [rad]');
grid on

% Metodo alternativo usando funzioni per risolvere equazioni differenziali
%
% syms theta(t) t
% 
% % Parameters
% theta0 = 2;
% omega0 = 16;
% C_g_lib = vpa(subs(C_g,F,0),3); % Moto libero
% 
% Dtheta = diff(theta,t,1);
% Ddtheta = diff(theta,t,2);
% 
% ode = Ddtheta * m_g + Dtheta * c_g + theta * k_g == C_g_lib;
% cond1 = theta(0) == theta0;
% cond2 = Dtheta(0) == omega0;
% conds = [cond1, cond2];
% theta(t) = dsolve(ode,conds);
% theta = simplify(theta);
% 
% fplot(theta,[0,10]);

% b)

c2_g = c1_g / 2;
csi2 = c2_g / (2 * sqrt(m_g * k_g));

lambda2 = roots([m_g c2_g k_g]);

X2 = (omega0 - lambda2(1) * theta0)/(lambda2(2) - lambda2(1));
X1 = theta0 - X2;

theta2 = X1*exp(lambda2(1)*t) + X2*exp(lambda2(2)*t);

% Grafico
figure
plot(t,theta2);
axis([-inf, inf, min(theta2)*1.1, max(theta2)*1.1]);
title('Time response of system free motion (halved damping ratio)');
xlabel('time [s]'); ylabel('displacement [rad]');
grid on

% c)

csi3 = 1.2;
c3_g = csi3 * 2 * m_g * omega_n;

lambda3 = roots([m_g c3_g k_g]);

X2 = (omega0 - lambda2(1) * theta0)/(lambda2(2) - lambda2(1));
X1 = theta0 - X2;

theta3 = X1*exp(lambda3(1)*t) + X2*exp(lambda3(2)*t);

% Grafico
figure
plot(t,theta3);
axis([-inf, inf, min(real(theta3))*1.1, max(theta3)*1.1]);
title('Time response of system free motion (damping ratio = 1.2)');
xlabel('time [s]'); ylabel('displacement [rad]');
grid on

%% 3) Moto forzato

% a)

omega = 0:0.01:5*omega_d;

% Caso 1: 2a)

FRF1 = 1./(-m_g*omega.^2 + 1i*omega*c1_g + k_g);

figure

subplot(2,1,1)
plot(omega, abs(FRF1));
axis([-inf, inf, -inf, inf]);
title('Frequency Response Function modulus');
xlabel('\omega [rad/sample]'); ylabel('$|\frac{\tilde{X_0}(\omega)}{F_0}|$', 'Interpreter', 'LaTeX');
grid on

subplot(2,1,2)
plot(omega, angle(FRF1));
axis([-inf, inf, -pi, pi]);
title('Frequency Response Function phase');
xlabel('\omega [rad/sample]'); ylabel('$\angle(\frac{\tilde{X_0}(\omega)}{F_0}) [rad]$', 'Interpreter', 'LaTeX');
grid on

% Caso 2: 2b)

FRF2 = 1./(-m_g*omega.^2 + 1i*omega*c2_g + k_g);

figure

subplot(2,1,1)
plot(omega, abs(FRF2));
axis([-inf, inf, -inf, inf]);
title('Frequency Response Function modulus (halved damping ratio)');
xlabel('\omega [rad/sample]'); ylabel('$|\frac{\tilde{X_0}(\omega)}{F_0}|$', 'Interpreter', 'LaTeX');
grid on

subplot(2,1,2)
plot(omega, angle(FRF2));
axis([-inf, inf, -pi, pi]);
title('Frequency Response Function phase (halved damping ratio)');
xlabel('\omega [rad/sample]'); ylabel('$\angle(\frac{\tilde{X_0}(\omega)}{F_0}) [rad]$', 'Interpreter', 'LaTeX');
grid on

% Caso 3: 2c)

FRF3 = 1./(-m_g*omega.^2 + 1i*omega*c3_g + k_g);

figure

subplot(2,1,1)
plot(omega, abs(FRF3));
axis([-inf, inf, -inf, inf]);
title('Frequency Response Function modulus (damping ratio = 1.2)');
xlabel('\omega [rad/sample]'); ylabel('$|\frac{\tilde{X_0}(\omega)}{F_0}|$', 'Interpreter', 'LaTeX');
grid on

subplot(2,1,2)
plot(omega, angle(FRF3));
axis([-inf, inf, -pi, pi]);
title('Frequency Response Function phase (damping ratio = 1.2)');
xlabel('\omega [rad/sample]'); ylabel('$\angle(\frac{\tilde{X_0}(\omega)}{F_0}) [rad]$', 'Interpreter', 'LaTeX');
grid on

% b)

% Parametri dell'eccitazione armonica per oscillazione totale (libera +
% forzata)
A = 2.5;
f_case1 = 0.15; f_case2 = 4.5; % 2 casi
phi = pi/3;

% c)

% Parametri dell'eccitazione armonica per oscillazione forzata
B1 = 1.2; B2 = 0.5; B3 = 5;
f1 = 0.1; f2 = 0.6; f3 = 3.3;
phi1 = pi/4; phi2 = pi/5; phi3 = pi/6;
