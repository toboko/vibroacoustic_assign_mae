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
x3_0 = 0.1; theta1_0 = pi/12; theta2_0 = -pi/12; 
x_0 = [x3_0; theta1_0; theta2_0];
x3_0_dot = 1; theta1_0_dot = 0.5; theta2_0_dot = 2;
x_0_dot = [x3_0_dot; theta1_0_dot; theta2_0_dot];

% Matrices

% da aggiornare
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

stepsize = 100;


%% 1. Equations of motion and system matrices

%% Eigenfrequencies & eigenvectors

%% 1.b
% Evaluate the eigenfrequencies and corresponding 
% eigenvectors in case of undamped and damped system.

%% UNDAMPED CASE

% Method: eigenvalues-eigenvectors problem 
% lambda = i*omega;
%[V_und,D_und] = eig(inv(M_gen)*K_gen); % V are the eigenvectors, D are the eigenvalues
[V_und,D_und] = eig(-M_gen\K_gen);

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

%% 3. Forced motion of the system 
%     (considering the Rayleigh damping as in 1.c)

%% a. Plot and comment the elements of the Frequency Response Matrix 

% DAMPED - Frequency Responce Function (FRF)
wMax = 50;
w = 0:0.01:wMax; 

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


%% c. Plot the co-located FRF between the rotation of the disk of radius ?2 and 
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

t = 0:1/stepsize:20;
omega_1 = 2*pi*f_1;
omega_2 = 2*pi*f_2;

omega_1_index = round(omega_1 / (1/stepsize)) + 1;
omega_2_index = round(omega_2 / (1/stepsize)) + 1;

for ii = 1:length(w)
    FRF = inv(-w(ii)^2*M_gen+1i*w(ii)*C_Rayleigh+K_gen)*deltaA;
    FRF11_D(ii) = FRF(1,1); FRF21_D(ii) = FRF(2,1); FRF31_D(ii) = FRF(3,1);
end

%FREE MOTION FAKE
theta1_free = 0;
theta2_free = 0;
x3_free = 0 ;


%Complete time response
theta1_forced = A_1*R_1*abs(FRF11_D(omega_1_index))*cos(omega_1*t + angle(FRF11_D(omega_1_index))) + ...
    A_2*R_1*abs(FRF11_D(omega_2_index))*cos(omega_2*t + angle(FRF11_D(omega_2_index)));
theta1_caseD = theta1_free + theta1_forced;

theta2_forced = A_1*R_2*abs(FRF21_D(omega_1_index))*cos(omega_1*t + angle(FRF21_D(omega_1_index))) + ...
    A_2*R_2*abs(FRF21_D(omega_2_index))*cos(omega_2*t + angle(FRF21_D(omega_2_index)));
theta2_caseD = theta2_free + theta2_forced;

x3_forced = A_1*abs(FRF31_D(omega_1_index))*cos(omega_1*t + angle(FRF31_D(omega_1_index))) + ...
    A_2*abs(FRF31_D(omega_2_index))*cos(omega_2*t + angle(FRF31_D(omega_2_index)));
x3_caseD = x3_free + x3_forced;

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
    omegaE_index(kk+1) = round(omegaE(kk+1) / (1/stepsize)) + 1;
    aE(kk+1) = (-1)^kk / (2*kk+1);
end

aE = aE* (8/pi^2);

t = 10:1/stepsize:30;

tr_wave = zeros(1,length(t));
%triangular wave 
for kk = 1:5
    tr_wave = tr_wave + aE(kk)*sin(omegaE(kk)*t);
end

% Calcolare steady state response di A
HA_caseE = FRF_A(omegaE_index);

xA_free_caseE = zeros(1, length(t));
for kk = 1:5
    xA_free_caseE = xA_free_caseE + aE(kk)*abs(HA_caseE(kk))*sin( ...
        omegaE(kk)*t + angle(HA_caseE(kk)));
end


% plottare xa e plottare triang wave
figure('Name', 'CASE 3.E')
sgtitle('Steady-State Time Response'); 
subplot(2,1,1)
plot(t,tr_wave, 'Color', 'green');
xlabel('Time [s]','FontSize',8); ylabel('F(t) - (triangular wave) [N]','FontSize',8);
subplot(2,1,2)
plot(t,xA_free_caseE, 'Color', 'magenta');
xlabel('Time [s]','FontSize',8); ylabel('x_A [m]','FontSize',8);

