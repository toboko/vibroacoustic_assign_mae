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

% Passo 1: definizione energie

% E_c = (1/2) * J1 * omega1^2 + (1/2) * M2 * v2;

% V_e = (1/2) * k1 * DELTA_L1^2 + (1/2) * k2 * DELTA_L2^2;
% V_g = M2 * g * h2;
% V = V_e + V_g;

% D = (1/2) * c1 * DELTA_L1_dot^2 + (1/2) * c2 * DELTA_L2_dot^2;

% delta_W = F * delta_y2; % Assumiamo una forza verticale diretta verso
% l'alto applicata su M2

% Passo 2: espressione variabili dipendenti in funzione di quelle
% indipendenti

% omega1 = theta_dot;
% v2 = y2_dot = omega1 * R2 = theta_dot * R2;
% DELTA_L2 = (Rivals) = - theta * R1;
% DELTA_L2_dot = - theta_dot * R1;
% DELTA_L1 = (Rivals) = theta * R2;
% DELTA_L1_dot = theta_dot * R2;
% h2 = (livello potenziale gravitazionale = 0 al punto di equilibrio) = y2
% = theta * R2;
% delta_y2 = delta_theta * R2;

% Passo 3: sostituzione nell'equazione di Lagrange

% delta_W = Q_theta * delta_theta = F * delta_y2 = F * delta_theta * R2 => Q_theta = F * R2;

% d/d_t(d/d_theta_dot(E_c)) - d/d_theta(E_c) + d/d_theta(V) + d/d_theta_dot(D) =
% Q_theta

% b)

% c)


%% 2) Moto libero

% a)

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
