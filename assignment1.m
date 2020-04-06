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
c_g = c1 * R2^2 + c1 * R1^2;
k_g = k1 * R2^2 + k2 * R1^2;
C_g = F * R2 + M2 * g * R2;

% b)

csi = c_g / (2 * sqrt(m_g * k_g));

% c)

omega_n = sqrt(k_g / m_g);
alpha = csi * omega_n;
omega_d = sqrt(omega_n^2 - alpha^2);

%% 2) Moto libero

% a)

lambda = roots([m_g c_g k_g]);

% Condizioni iniziali
theta0 = 2;
omega0 = 16;

t = 0:0.01:10;

X2 = (omega0 - lambda(1) * theta0)/(lambda(2) - lambda(1));
X1 = theta0 - X2;

theta = X1*exp(lambda(1)*t) + X2*exp(lambda(2)*t);

% Grafico
figure
plot(t,theta);
axis([-inf, inf, min(theta)*1.1, max(theta)*1.1]);
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

% c)


%% 3) Moto forzato

% a)

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
