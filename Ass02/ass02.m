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

% x = [ theta1 theta2 x3]

% Init Conditions
theta1_0 = pi/12; theta2_0 = -pi/12; x3_0 = 0.1;
x_0 = [theta1_0; theta2_0; x3_0];
theta1_0_dot = 0.5; theta2_0_dot = 2; x3_0_dot = 1;
x_0_dot = [theta1_0_dot; theta2_0_dot; x3_0_dot;];

% Matrices

jac_M = [0   R_2  1;
         1   0   0;
         0   0   1;
         0   1   0; 
         0   0   1 ]; 

M = diag([M_1 J_1 M_2 J_2 M_3]);
M_gen = jac_M'*M*jac_M;

jac_K = [0  -R_2    -1;
         -R_1 2*R_2  0;
         R_1   R_2   0;
         0     0     1];

K = diag([k1 k2 k3 k4]);
K_gen = jac_K'*K*jac_K;


jac_C = jac_K; 

C = diag([c1 c2 c3 c4]);
C_gen = jac_C'*C*jac_C;

%Utils
Ndof = 3;
T_max = 20;
f_s = 100;
t = 0:1/f_s:T_max;

%wMax = pi*f_s;
wMax = 50;
w = 0:0.01:wMax; 

%% 1. Equations of motion and system matrices

%% Eigenfrequencies & eigenvectors

%% 1.b
% Evaluate the eigenfrequencies and corresponding 
% eigenvectors in case of undamped and damped system.

%% UNDAMPED CASE

% Method: eigenvalues-eigenvectors problem 
% lambda = i*omega;

[V_und,D_und] = eig(M_gen\K_gen); % Equivalent to eig(-inv(M_gen)*K_gen), but more efficient. 
% V are the eigenvectors, D are the eigenvalues


w_nat = sqrt(diag(D_und)); %natural frequencies
 
%% Normalization
Vn_und = [V_und]./[V_und(1,:)];

%% Sorting the solutions
w_nat_ord = [w_nat(1); -w_nat(1); w_nat(2); -w_nat(2); w_nat(3); -w_nat(3)];
V_ord = [Vn_und(:,1), Vn_und(:,1), Vn_und(:,2), Vn_und(:,2), Vn_und(:,3), Vn_und(:,3)]; %all vibration modes

%% DAMPED CASE 

% in [A_damp] the submatrix [C_gen] is full
% use state form matrix: new variable z = [vel; pos];
A_damp = -inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*[C_gen K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp); % Modeshapes and eigenfrequencies
V_damp = V_damp(4:6,:); %consider only the displacements, all vibration modes
lambda = diag(D_damp); %complex and conjugate

%% Normalization
Vn_damp = V_damp./V_damp(1,:);

%% 1.c
% Assuming Rayleigh damping, evaluate alpha and beta to approximate the generalized damping
% matrix [C*] to be of the form alpha[M]+beta[K]

for i = 1:Ndof
    xi(i,1) = C_gen(i,i)/(2*w_nat(i)*M_gen(i,i));
end


tmp = [1/(2*w_nat(1)) w_nat(1)/2;
       1/(2*w_nat(2)) w_nat(2)/2;
       1/(2*w_nat(3)) w_nat(3)/2];

ab = pinv(tmp)*xi;
alpha = ab(1); beta = ab(2);

C_Rayleigh = alpha.*M_gen + beta.*K_gen; 


%% 2. Free motion of the system (considering the Rayleigh damping as in 1.c)

lambdas = [lambda(1); lambda(3); lambda(5)]; % Extracts each mode's eigenvalue
alphas = abs(real(lambdas)); % Extracts each mode's damping factor
omegas_d = imag(lambdas); % Extracts each mode's damped frequency
modeshapes = [Vn_damp(:,1), Vn_damp(:,3), Vn_damp(:,5)]; % Extracts each mode's modeshape

%% 2.A - Plot and comment the free motion of the system

init_cond = [x_0; x_0_dot];

S = [abs(modeshapes) .* cos(angle(modeshapes))                                                  abs(modeshapes) .* sin(angle(modeshapes))
     abs(modeshapes) .* (-alphas'.*cos(angle(modeshapes))-omegas_d'.*sin(angle(modeshapes)))    abs(modeshapes) .* (-alphas'.*sin(angle(modeshapes))+omegas_d'.*cos(angle(modeshapes)))];

%S = [modeshapes             zeros(3);
%    -modeshapes.*alphas'    modeshapes.*omegas_d']; % Known part (depending on system parameters) of the free response's sinusoidal and cosinusoidal terms coefficients

AB = S \ init_cond; % Inverse formula of init_cond = S * AB. N.B.: S\init_cond =inv(S)*init_cond. AB is the unknown part (A_i and B_i, depending on initial conditions) of the free response's sinusoidal and cosinusoidal terms coefficients

A = AB(1:3); %A_i coefficients
B = AB(4:end); %B_i coefficients


% Alternative method using automatic equation system solver
% 
% syms A1 A2 A3 B1 B2 B3
% 
% eqn1 = theta1_0 == modeshapes(1,1)*A1 + modeshapes(1,2)*A2 + modeshapes(1,3)*A3;
% eqn2 = theta2_0 == modeshapes(2,1)*A1 + modeshapes(2,2)*A2 + modeshapes(2,3)*A3;
% eqn3 = x3_0 == modeshapes(3,1)*A1 + modeshapes(3,2)*A2 + modeshapes(3,3)*A3;
% 
% eqn4 = theta1_0_dot == modeshapes(1,1)*(-alphas(1)*A1 + omegas_d(1)*B1) + modeshapes(1,2)*(-alphas(2)*A2 + omegas_d(2)*B2) + modeshapes(1,3)*(-alphas(3)*A3 + omegas_d(3)*B3);
% eqn5 = theta2_0_dot == modeshapes(2,1)*(-alphas(1)*A1 + omegas_d(1)*B1) + modeshapes(2,2)*(-alphas(2)*A2 + omegas_d(2)*B2) + modeshapes(2,3)*(-alphas(3)*A3 + omegas_d(3)*B3);
% eqn6 = x3_0_dot == modeshapes(3,1)*(-alphas(1)*A1 + omegas_d(1)*B1) + modeshapes(3,2)*(-alphas(2)*A2 + omegas_d(2)*B2) + modeshapes(3,3)*(-alphas(3)*A3 + omegas_d(3)*B3);
% 
% [A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [A1, A2, A3, B1, B2, B3]);
% 
% AB = vpa(linsolve(A,B));

% Building up the free response manually

% Free response of each independent variable: one matrix for each of them
% with one row for the response of each mode
theta1_free = zeros(3,length(t));
theta2_free = zeros(3,length(t));
x3_free = zeros(3,length(t));

for i=1:3
    theta1_free(i,:) = exp(-alphas(i)*t) * abs(modeshapes(1,i)) .* (A(i)*cos(omegas_d(i)*t + angle(modeshapes(1,i))) + B(i)*sin(omegas_d(i)*t + angle(modeshapes(1,i))));
    theta2_free(i,:) = exp(-alphas(i)*t) * abs(modeshapes(2,i)) .* (A(i)*cos(omegas_d(i)*t + angle(modeshapes(2,i))) + B(i)*sin(omegas_d(i)*t + angle(modeshapes(2,i))));
    x3_free(i,:) = exp(-alphas(i)*t) * abs(modeshapes(3,i)) .* (A(i)*cos(omegas_d(i)*t + angle(modeshapes(3,i))) + B(i)*sin(omegas_d(i)*t + angle(modeshapes(3,i))));
end

% Free response of each independent variable, summing together each mode's
% response (so, columnwise) for each independent variable
x_free = [sum(theta1_free,1);
         sum(theta2_free,1);
         sum(x3_free,1)];
     
% Grafico
figure('Name','2.A')

subplot(3,1,1)
plot(t,x_free(1,:));
axis([-inf, inf, min(x_free(1,:))*1.1, max(x_free(1,:))*1.1]);
title('Time response of \theta_1 free motion','FontSize',6);
xlabel('time t [s]','FontSize',6); ylabel('rotation \theta_{1_{free}}(t)   [rad]','FontSize',6);
grid on

subplot(3,1,2)
plot(t,x_free(2,:));
axis([-inf, inf, min(x_free(2,:))*1.1, max(x_free(2,:))*1.1]);
title('Time response of \theta_2 free motion','FontSize',6);
xlabel('time t [s]','FontSize',6); ylabel('rotation \theta_{2_{free}}(t)   [rad]','FontSize',6);
grid on

subplot(3,1,3)
plot(t,x_free(3,:));
axis([-inf, inf, min(x_free(3,:))*1.1, max(x_free(3,:))*1.1]);
title('Time response of x_3 free motion','FontSize',6);
xlabel('time t [s]','FontSize',6); ylabel('displacement x_{3_{free}}(t)   [m]','FontSize',6);
grid on

%% 2.b - Impose particular initial conditions so that only one mode contributes to the free motion of the system.

% A_i and B_i coefficients exciting 2nd mode only
A_II = [0; 1; 0];
B_II = [0; 1; 0];

AB_II = [A_II; B_II];

init_cond_II = S * AB_II;

% Building up the free response for each variable
theta1_free_II = zeros(3,length(t));
theta2_free_II = zeros(3,length(t));
x3_free_II = zeros(3,length(t));

for i=1:3
    theta1_free_II(i,:) = modeshapes(1,i) * (exp(-alphas(i)*t) .* (A_II(i)*cos(omegas_d(i)*t) + B_II(i)*sin(omegas_d(i)*t)));
    theta2_free_II(i,:) = modeshapes(2,i) * (exp(-alphas(i)*t) .* (A_II(i)*cos(omegas_d(i)*t) + B_II(i)*sin(omegas_d(i)*t)));
    x3_free_II(i,:) = modeshapes(3,i) * (exp(-alphas(i)*t) .* (A_II(i)*cos(omegas_d(i)*t) + B_II(i)*sin(omegas_d(i)*t)));
end

x_free_II = [sum(theta1_free_II,1);
            sum(theta2_free_II,1);
            sum(x3_free_II,1)];
     
% Grafico
figure('Name','2.b')

subplot(3,1,1)
plot(t,real(x_free_II(1,:)));
axis([-inf, inf, min(real(x_free_II(1,:)))*1.1, max(real(x_free_II(1,:)))*1.1]);
title('Time response of \theta_1 free motion, 2^{nd} mode forced only','FontSize',6);
xlabel('time t [s]','FontSize',6); ylabel('rotation \theta_{1_{free}}(t)   [rad]','FontSize',6);
grid on

subplot(3,1,2)
plot(t,real(x_free_II(2,:)));
axis([-inf, inf, min(real(x_free_II(2,:)))*1.1, max(real(x_free_II(2,:)))*1.1]);
title('Time response of \theta_2 free motion, 2^{nd} mode forced only','FontSize',6);
xlabel('time t [s]','FontSize',6); ylabel('rotation \theta_{2_{free}}(t)   [rad]','FontSize',6);
grid on

subplot(3,1,3)
plot(t,real(x_free_II(3,:)));
axis([-inf, inf, min(real(x_free_II(3,:)))*1.1, max(real(x_free_II(3,:)))*1.1]);
title('Time response of x_3 free motion, 2^{nd} mode forced only','FontSize',6);
xlabel('time t [s]','FontSize',6); ylabel('displacement x_{3_{free}}(t)   [m]','FontSize',6);
grid on

%% 3. Forced motion of the system 
%     (considering the Rayleigh damping as in 1.c)

%% a. Plot and comment the elements of the Frequency Response Matrix 

% DAMPED - Frequency Responce Function (FRF)

% [H(w)] = [D(w)]^-1 = [-w^2*[M]+i*w*[C]+[K]]^-1
for ii = 1:length(w)
    FRF = inv(-w(ii)^2*M_gen+1i*w(ii)*C_Rayleigh+K_gen);
    FRF11(ii) = FRF(1,1); FRF12(ii) = FRF(1,2); FRF13(ii) = FRF(1,3);
    FRF21(ii) = FRF(2,1); FRF22(ii) = FRF(2,2); FRF23(ii) = FRF(2,3);
    FRF31(ii) = FRF(3,1); FRF32(ii) = FRF(3,2); FRF33(ii) = FRF(3,3);
end

% Plot of the FRF modulus
figure('Name', 'Case 3.a (1)')
subplot(3,3,1); plot(w/2/pi,abs(FRF11),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_1')
subplot(3,3,2); plot(w/2/pi,abs(FRF12),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_2| [m/N]'); title('H_1_2')
subplot(3,3,3); plot(w/2/pi,abs(FRF13),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_3| [m/N]'); title('H_1_3')
subplot(3,3,4); plot(w/2/pi,abs(FRF21),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_1| [m/N]'); title('H_2_1')
subplot(3,3,5); plot(w/2/pi,abs(FRF22),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_2| [m/N]'); title('H_2_2')
subplot(3,3,6); plot(w/2/pi,abs(FRF23),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_3| [m/N]'); title('H_2_3')
subplot(3,3,7); plot(w/2/pi,abs(FRF31),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_1| [m/N]'); title('H_3_1')
subplot(3,3,8); plot(w/2/pi,abs(FRF32),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_2| [m/N]'); title('H_3_2')
subplot(3,3,9); plot(w/2/pi,abs(FRF33),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [m/N]'); title('H_3_3')
sgtitle('Modulus of FRF')

% Plot of the FRF phase
figure('Name', 'Case 3.a (2)')
subplot(3,3,1); plot(w/2/pi,angle(FRF11)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_1 [deg]'); title('H_1_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,2); plot(w/2/pi,angle(FRF12)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_2 [deg]'); title('H_1_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,3); plot(w/2/pi,angle(FRF13)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_3 [deg]'); title('H_1_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,4); plot(w/2/pi,angle(FRF21)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_1 [deg]'); title('H_2_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,5); plot(w/2/pi,angle(FRF22)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_2 [deg]'); title('H_2_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,6); plot(w/2/pi,angle(FRF23)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_3 [deg]'); title('H_2_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,7); plot(w/2/pi,angle(FRF31)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_1 [deg]'); title('H_3_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,8); plot(w/2/pi,angle(FRF32)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_2 [deg]'); title('H_3_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,9); plot(w/2/pi,angle(FRF33)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_3 [deg]'); title('H_3_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Phase of FRF')

%% b. Plot the co-located FRF of point A at the centre of the disk 
%  (co-located meaning the FRF between the displacement in A and a force applied in A).

% Lagrangian component of the force 
deltaA = [0; R_2; 1]; 

% Co-located FRF of point A

for ii = 1:length(w)
    FRF = inv(-w(ii)^2*M_gen+1i*w(ii)*C_Rayleigh+K_gen);
    FRF_A(ii) = deltaA'*FRF*deltaA;
end

% Plot of the FRF modulus Co-located in point A
figure('Name', 'Case 3.b')
subplot(2,1,1); 
plot(w/2/pi,abs(FRF_A),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_A| [m/N]'); 

% Plot of the FRF phase
subplot(2,1,2); 
plot(w/2/pi,angle(FRF_A)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('\angle H_AA [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Co-located FRF in point A')


%% c. Plot the co-located FRF between the rotation of the disk of radius R_2 and 
% the torque applied onto the disk itself.

% Lagrangian component of the force 
deltaR2 = [0; R_2; 0]; 

% Co-located FRF of the disk 2

for ii = 1:length(w)
    FRF = inv(-w(ii)^2*M_gen+1i*w(ii)*C_Rayleigh+K_gen);
    FRF_R2(ii) = deltaR2'*FRF*deltaR2;
end

% Plot of the FRF modulus Co-located of the disk 2
figure('Name', 'Case 3.c')
subplot(2,1,1); 
plot(w/2/pi,abs(FRF_R2),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_R_2| [m/N]'); 

% Plot of the FRF phase
subplot(2,1,2); 
plot(w/2/pi,angle(FRF_R2)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('\angle H_R_2 [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Co-located FRF of the disk 2')

%% d. Starting from the initial condition defined in (2.a), evaluate the 
% COMPLETE time response of the system for the three degrees of freedom 
% to the horizontal force applied in A, considering that 
% the force is harmonic of the form:
% F_A = A_1*cos(2*pi*f_1*t) + A_2*cos(2*pi*f_2*t)

%t = 0:1/f_s:20;
omega_1 = 2*pi*f_1;
omega_2 = 2*pi*f_2;

omega_1_index = round(omega_1 / (1/f_s)) + 1;
omega_2_index = round(omega_2 / (1/f_s)) + 1;

for ii = 1:length(w)
    FRF = inv(-w(ii)^2*M_gen+1i*w(ii)*C_Rayleigh+K_gen)*deltaA;
    FRF11_D(ii) = FRF(1,1); FRF21_D(ii) = FRF(2,1); FRF31_D(ii) = FRF(3,1);
end

%Complete time response
theta1_forced = A_1*R_1*abs(FRF11_D(omega_1_index))*cos(omega_1*t + angle(FRF11_D(omega_1_index))) + ...
    A_2*R_1*abs(FRF11_D(omega_2_index))*cos(omega_2*t + angle(FRF11_D(omega_2_index)));
theta1_caseD = sum(theta1_free) + theta1_forced;

theta2_forced = A_1*R_2*abs(FRF21_D(omega_1_index))*cos(omega_1*t + angle(FRF21_D(omega_1_index))) + ...
    A_2*R_2*abs(FRF21_D(omega_2_index))*cos(omega_2*t + angle(FRF21_D(omega_2_index)));
theta2_caseD = sum(theta2_free) + theta2_forced;

x3_forced = A_1*abs(FRF31_D(omega_1_index))*cos(omega_1*t + angle(FRF31_D(omega_1_index))) + ...
    A_2*abs(FRF31_D(omega_2_index))*cos(omega_2*t + angle(FRF31_D(omega_2_index)));
x3_caseD = sum(x3_free) + x3_forced;

% PLOT
figure('Name', 'Case 3.D')
sgtitle('Complete Time Response'); 

subplot(3,1,1)
plot(t, theta1_caseD)
xlabel('Time [s]')
ylabel('\theta_1')

subplot(3,1,2)
plot(t, theta2_caseD)
xlabel('Time [s]')
ylabel('\theta_2')

subplot(3,1,3)
plot(t, x3_caseD)
xlabel('Time [s]')
ylabel('x_3 [m]')



%% e. Evaluate the steady-state response of the system of horizontal displacement of point A to
% the horizontal force applied in A considering that the force is a periodic triangular wave, with
% fundamental frequency ?0, 
% f_0 = 0.75;

for kk = 0:4
    omegaE(kk+1) =2*pi*f_0*(2*kk+1);
    omegaE_index(kk+1) = round(omegaE(kk+1) / (1/f_s)) + 1;
    aE(kk+1) = (-1)^kk / (2*kk+1);
end

aE = aE* (8/pi^2);

t_ss = t(round(length(t)/2):end);

tr_wave = zeros(1,length(t_ss));
%triangular wave 
for kk = 1:5
    tr_wave = tr_wave + aE(kk)*sin(omegaE(kk)*t_ss);
end

% Calcolare steady state response di A
HA_caseE = FRF_A(omegaE_index);

xA_free_caseE = zeros(1, length(t_ss));
for kk = 1:5
    xA_free_caseE = xA_free_caseE + aE(kk)*abs(HA_caseE(kk))*sin( ...
        omegaE(kk)*t_ss + angle(HA_caseE(kk)));
end


% plottare xa e plottare triang wave
figure('Name', 'CASE 3.E')
sgtitle('Steady-State Time Response'); 
subplot(2,1,1)
plot(t_ss,tr_wave, 'Color', 'green');
xlabel('Time [s]','FontSize',8); ylabel('F(t) - (triangular wave) [N]','FontSize',8);
subplot(2,1,2)
plot(t_ss,xA_free_caseE, 'Color', 'magenta');
xlabel('Time [s]','FontSize',8); ylabel('x_A [m]','FontSize',8);

%% 4) Modal approach (considering the Rayleigh damping as in 1.c)

%% a. Derive the equations of motion in modal coordinates and plot the elements of the corresponding Frequency Response Matrix 

%matrix of the mode shapes
%phi = [V_1(:,1), V_1(:,2), V_1(:,3)]; 
phi = [Vn_damp(:,1), Vn_damp(:,3), Vn_damp(:,5)]; 
%phi = [Vn_damp(:,2), Vn_damp(:,4), Vn_damp(:,6)]; 

% Diagonalization of the matrix
M_mod = phi'*M_gen*phi;
K_mod = phi'*K_gen*phi;

C_mod = phi'*C_Rayleigh*phi;

% Put equal to 0 the extra diagonal terms (< eps = 2.2204e-16)
diag_m = diag(M_mod);
MM = [diag(diag_m)];
diag_k = diag(K_mod);
KK = [diag(diag_k)];
diag_c = diag(C_mod);
CC = [diag(diag_c)];

%% b. Reconstruct the co-located FRF of point A employing a modal approach and compare with the one obtained using physical coordinates in (3.b).
% w = 0:0.01:30;
%deltaA = [0; R_2; 1];
deltaA_q = inv(phi)*deltaA;
deltaR2_q = inv(phi)*deltaR2;

% Diagonal matrix with eps in extra diagonal positions
for ii = 1:length(w)
    FRF_q = inv(-w(ii)^2*M_mod+1i*w(ii)*C_mod+K_mod);
    FRF_qA(ii) = deltaA_q'*FRF_q*deltaA_q;
    FRF_qR2(ii) = deltaR2_q'*FRF_q*deltaR2_q;


    FRF11_q(ii) = FRF_q(1,1); FRF12_q(ii) = FRF_q(1,2); FRF13_q(ii) = FRF_q(1,3);
    FRF21_q(ii) = FRF_q(2,1); FRF22_q(ii) = FRF_q(2,2); FRF23_q(ii) = FRF_q(2,3);
    FRF31_q(ii) = FRF_q(3,1); FRF32_q(ii) = FRF_q(3,2); FRF33_q(ii) = FRF_q(3,3);
end

% Plot of the FRF_q modulus (4.A)
figure('Name', 'Case 4.A (1)')
subplot(3,3,1); plot(w/2/pi,abs(FRF11_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_1')
subplot(3,3,2); plot(w/2/pi,abs(FRF12_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_2| [m/N]'); title('H_1_2')
subplot(3,3,3); plot(w/2/pi,abs(FRF13_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_3| [m/N]'); title('H_1_3')
subplot(3,3,4); plot(w/2/pi,abs(FRF21_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_1| [m/N]'); title('H_2_1')
subplot(3,3,5); plot(w/2/pi,abs(FRF22_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_2| [m/N]'); title('H_2_2')
subplot(3,3,6); plot(w/2/pi,abs(FRF23_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_3| [m/N]'); title('H_2_3')
subplot(3,3,7); plot(w/2/pi,abs(FRF31_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_1| [m/N]'); title('H_3_1')
subplot(3,3,8); plot(w/2/pi,abs(FRF32_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_2| [m/N]'); title('H_3_2')
subplot(3,3,9); plot(w/2/pi,abs(FRF33_q),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [m/N]'); title('H_3_3')
sgtitle('Modulus of FRF - Modal Approach')

% Plot of the FRF_q phase (4.A)
figure('Name', 'Case 4.A (2)')
subplot(3,3,1); plot(w/2/pi,angle(FRF11_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_1 [deg]'); title('H_1_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,2); plot(w/2/pi,angle(FRF12_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_2 [deg]'); title('H_1_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,3); plot(w/2/pi,angle(FRF13_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_3 [deg]'); title('H_1_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,4); plot(w/2/pi,angle(FRF21_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_1 [deg]'); title('H_2_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,5); plot(w/2/pi,angle(FRF22_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_2 [deg]'); title('H_2_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,6); plot(w/2/pi,angle(FRF23_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_3 [deg]'); title('H_2_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,7); plot(w/2/pi,angle(FRF31_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_1 [deg]'); title('H_3_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,8); plot(w/2/pi,angle(FRF32_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_2 [deg]'); title('H_3_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,9); plot(w/2/pi,angle(FRF33_q)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_3 [deg]'); title('H_3_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Phase of FRF - Modal Approach')

% Plot 4.B
figure('Name', 'Case 4.B')
subplot(2,1,1); 
plot(w/2/pi,abs(FRF_qA),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_A| [m/N]'); 
subplot(2,1,2);
plot(w/2/pi,angle(FRF_qA)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('\angle H_A [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Co-located FRF in point A - Modal Approach)')

%% 4.C
% Plot of the FRF modulus Co-located of the disk 2 (4.c)
figure('Name', 'Case 4.c')
subplot(2,1,1); 
plot(w/2/pi,abs(FRF_qR2),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_R_2| [m/N]'); 

% Plot of the FRF phase Co-located of the disk 2 (4.c)
subplot(2,1,2); 
plot(w/2/pi,angle(FRF_qR2)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('\angle H_R_2 [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Co-located FRF of the disk 2 - Modal Approach')
