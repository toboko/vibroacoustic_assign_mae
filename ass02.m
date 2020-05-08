%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 02       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, clear all, clc

%% DATA

M_1 = 5; J_1 = 2.5; R_1 = 1;
M_2 = 1.25; J_2 = 0.16; R_2 = 0.5; 
M_3 = 10;  

k1 = 1000; c1 = 0.5; 
k2 = 100; c2 = 0.5;
k3 = 560; c3 = 1;
k4 = 800; c4 = 4;

A_1 = 15; A_2 = 7; f_1 = 1.5; f_2 = 3.5; f_0 = 0.75;

%{
x = [ theta1 theta2 x3]
%}

% Init Conditions
x3_0 = 0.1; theta1_0 = pi/12; theta2_0 = -pi/12; 
x_0 = [x3_0; theta1_0; theta2_0];
x3_0_dot = 1; theta1_0_dot = 0.5; theta2_0_dot = 2;
x_0_dot = [x3_0_dot; theta1_0_dot; theta2_0_dot];

% Matrices

%da aggiornare
jac_M = [0   R_2  1;
         1   0   0;
         0   0   1;
         0   1   0; 
         0   0   1 ]; 

M = diag([M_1 J_1 M_2 J_2 M_3]);
M_gen = jac_M'*M*jac_M;

% da aggiornare
jac_K = [0  -R_2    -1;
         -R_1 2*R_2  0;
         R_1   R_2   0;
         0     0     1];

K = diag([k1 k2 k3 k4]);
K_gen = jac_K'*K*jac_K;


jac_C = jac_K; 

C = diag([c1 c2 c3 c4]);
C_gen = jac_C'*C*jac_C;


%% 1. Equations of motion and system matrices

%% Eigenfrequencies & eigenvectors

%% 1.b
% Evaluate the eigenfrequencies and corresponding 
% eigenvectors in case of undamped and damped system.

%% UNDAMPED CASE

% Method: eigenvalues-eigenvectors problem 
% lambda = i*omega;
[V_und,D_und] = eig(inv(M_gen)*K_gen); % V are the eigenvectors, D are the eigenvalues

w_nat = sqrt(diag(D_und)); %natural frequencies
 
%% Normalization
Vn_und = [V_und]./[V_und(1,:)];

%% Sorting the solutions
w_nat_ord = [w_nat(1); -w_nat(1); w_nat(2); -w_nat(2); w_nat(3); -w_nat(3)];
V_ord = [Vn_und(:,1), Vn_und(:,1), Vn_und(:,2), Vn_und(:,2), Vn_und(:,3), Vn_und(:,3)]; %all vibration modes

%% DAMPED CASE 

% in [A_damp] the submatrix [C] is full
% use state form matrix: new variable z = [vel; pos];
A_damp = - inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*[C_gen K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp);
V_damp = V_damp(3:4,:); %consider only the displacements
lambda = diag(D_damp); %complex and conjugate

%% 1.c
% Assuming Rayleigh damping, evaluate alpha and beta to approximate the generalized damping
% matrix [C*] to be of the form alpha[M]+beta[K]

for i = 1:length(w_nat)
    xi(i,1) = C(i,i)/(2*w_nat(i)*M(i,i));
end


tmp = [1/(2*w_nat(1)) w_nat(1)/2;
       1/(2*w_nat(2)) w_nat(2)/2;
       1/(2*w_nat(3)) w_nat(3)/2];

ab = pinv(tmp)*xi;
alpha = ab(1); beta = ab(2);
C_Rayleigh = alpha.*M_gen + beta.*K_gen; 


%% 2. Free motion of the system (considering the Rayleigh damping as in 1.c)



%% 2.A - Plot and comment the free motion of the system

