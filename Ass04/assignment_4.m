%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNDAMENTALS OF VIBRATION ANALYSIS AND VIBROACOUSTICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULE 2 - VIBROACOUSTICS OF MUSICAL INSTRUMENTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ASSIGNMENT 1 - AXIAL VIBRATION OF UNDAMPED AND DAMPED BARS  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

% Data
rho = 2700;
E = 70e9;
L = 2;
h = 0.05;
b = h;
f_max = 10000;

% Utils
%T = 10;
%f_s = 2 * f_max;
lambda_s = 0.005;
% %omega_max = 2*pi*f_max;
% omega_s = 2*pi*f_s;
% c = sqrt(E/rho);
% %k_max = omega_max/c;
% k_s = omega_s/c;
% %k_s = 2 * k_max;
% lambda_s = 2*pi/k_s;

%t = 0:1/f_s:T;
x = 0:lambda_s:L;


%% 1) Natural frequencies and mode shapes in free-fixed configuration

% Number of vibration modes up to f_max
i_free_fixed_max = floor(2*L * f_max * sqrt(rho/E) + 1/2);

% Mechanical waves velocity
c = sqrt(E/rho);

% Natural frequencies
f_nat_free_fixed = zeros(1,i_free_fixed_max);

for i=1:i_free_fixed_max
    f_nat_free_fixed(i) = (2*i-1)/(4*L) * c;
end

% Mode shapes
phi_free_fixed = zeros(i_free_fixed_max, length(x));

for i=1:i_free_fixed_max
    k_free_fixed_i = (2*i-1)*pi / (2*L);
    phi_free_fixed(i,:) = cos(k_free_fixed_i * x);
end

% Plots
figure('Name', '1) Mode shapes for free-fixed configuration')

subplot(i_free_fixed_max/2,2,1)
plot(x, real(phi_free_fixed(1,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('1^{st} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(1)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,2)
plot(x, real(phi_free_fixed(2,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('2^{nd} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(2)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,3)
plot(x, real(phi_free_fixed(3,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('3^{rd} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(3)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,4)
plot(x, real(phi_free_fixed(4,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('4^{th} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(4)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,5)
plot(x, real(phi_free_fixed(5,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('5^{th} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(5)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,6)
plot(x, real(phi_free_fixed(6,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('6^{th} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(6)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,7)
plot(x, real(phi_free_fixed(7,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('7^{th} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(7)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,8)
plot(x, real(phi_free_fixed(8,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('8^{th} mode shape - free-fixed bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(8)}}(x)','FontSize',12);
grid on


%% 2) Natural frequencies and mode shapes in free-free configuration

% Number of vibration modes up to f_max
i_free_free_max = floor(2*L * f_max * sqrt(rho/E));

% Natural frequencies
f_nat_free_free = zeros(1,i_free_free_max + 1);

for i=1:i_free_free_max + 1
    f_nat_free_free(i) = (i-1)/(2*L) * c;
end

% Mode shapes
phi_free_free = zeros(i_free_free_max, length(x));

for i=1:i_free_free_max + 1
    k_free_free_i = (i-1)*pi / L;
    phi_free_free(i,:) = cos(k_free_free_i * x);
end

% Plots
figure('Name', '2) Mode shapes for free-free configuration')

subplot(ceil(i_free_free_max/2),2,1)
plot(x, real(phi_free_free(1,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('1^{st} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(1)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,2)
plot(x, real(phi_free_free(2,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('2^{nd} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(2)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,3)
plot(x, real(phi_free_free(3,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('3^{rd} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(3)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,4)
plot(x, real(phi_free_free(4,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('4^{th} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(4)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,5)
plot(x, real(phi_free_free(5,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('5^{th} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(5)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,6)
plot(x, real(phi_free_free(6,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('6^{th} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(6)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,7)
plot(x, real(phi_free_free(7,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('7^{th} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(7)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,8)
plot(x, real(phi_free_free(8,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('8^{th} mode shape - free-free bar','FontSize',8);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(8)}}(x)','FontSize',12);
grid on


%% 3) FRFs in free-fixed configuration

% Custom function to get row and col length
rows = @(x) size(x,1); 
cols = @(x) size(x,2);

% Data
xL1 = L/2;
xL2 = L/5;
xL = [xL1, xL2];
S   = b * h;

% Case 1: undamped bar (standing wave Solution)

T_s   = 0.05;
f_s   = 1/T_s;
ll    = 10000;
f     = linspace(0, f_max, ll);
FRF1   = zeros(size(f));
FRF2   = zeros(size(f));

for x = 1:ll
    k      = 2 * pi * f(x) / c;
    FRF1(x) = sin(k * (L - xL(1))) / (k * E * S * cos(k * L));
    FRF2(x) = sin(k * (L - xL(2))) / (k * E * S * cos(k * L));
end

% Plot Undamped Bar (SWS)
figure('Name', 'Undamped Bar (SWS)');
subplot(2,2,1)
semilogy(f, abs(FRF1));
title('x = L/2');
xlabel('frequency [Hz]'); ylabel('FRF_x1 (f)');
grid on
subplot(2,2,2)
semilogy(f, abs(FRF2));
% axis([0, f_max, -0.00005, 0.00005]);
title('x = L/5');
xlabel('frequency [Hz]'); ylabel('FRF_x2 (f)');
grid on

subplot(2,2,3)
plot(f, angle(FRF1)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x1 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on
subplot(2,2,4)
plot(f, angle(FRF2)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x2 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on

% Case 2: damped bar (standing wave solution)

eta      = 0.01;
E_damp   = E * (1 + 1i * eta);

FRF1_damp = zeros(size(f));
FRF2_damp = zeros(size(f));

for x = 1:ll
    k_damp = (2 * pi* f(x) / c) * (1 - 1i * (eta / 2));
    FRF1_damp(x) = sin(k_damp * (L - xL(1))) / (k_damp * E_damp * S * cos(k_damp * L));
    FRF2_damp(x) = sin(k_damp * (L - xL(2))) / (k_damp * E_damp * S * cos(k_damp * L));
end

% Plot Damped Bar (SWS)
figure('Name', 'Damped Bar (SWS)');
subplot(2,2,1)
semilogy(f, abs(FRF1_damp));
title('x = L/2');
xlabel('frequency [Hz]'); ylabel('FRF_x1 (f)');
grid on
subplot(2,2,2)
semilogy(f, abs(FRF2_damp));
% axis([0, f_max, -0.00005, 0.00005]);
title('x = L/5');
xlabel('frequency [Hz]'); ylabel('FRF_x2 (f)');
grid on

subplot(2,2,3)
plot(f, angle(FRF1_damp)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x1 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on
subplot(2,2,4)
plot(f, angle(FRF2_damp)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x2 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on

% Damped bar (Modal Superposition Approach)
%f     = linspace(0.1, f_max, ll);
m = rho * S; % mass per unit length
M = zeros(i_free_fixed_max, i_free_fixed_max);
K = zeros(i_free_fixed_max, i_free_fixed_max);

k_nat = 2 * pi * f_nat_free_fixed / c;

for x = 1:i_free_fixed_max
    for y = 1:i_free_fixed_max
        kr = k_nat(x);
        ks = k_nat(y);
        if (x == y)
            M(x,y) = m * (L/2 + (((sin(2 * ks * L)) / (4 * ks))));
        else
            M(x,y) = m * ((sin(L*(kr - ks))/(2*(kr - ks))) + (sin(L*(kr + ks))/(2*(kr + ks))));
        end
    end
end

for x = 1:i_free_fixed_max
    for y = 1:i_free_fixed_max
        kr = k_nat(x);
        ks = k_nat(y);
        if (x == y)
            K(x,y) = (E*S*kr*ks) * (L/2 - sin(2*ks*L)/(4*ks));
        else
            K(x,y) = (E*S*kr*ks) * (sin(L*(kr - ks))/(2*(kr - ks)) - sin(L*(kr + ks))/(2*(kr + ks)));
        end
    end
end

FRF1_MSA   = zeros(size(f));
FRF2_MSA   = zeros(size(f));

for x = 1:ll
    for s = 1:i_free_fixed_max
        FRF1_MSA(x) = FRF1_MSA(x) + ((cos(k_nat(s)*xL(1)))/(((-1)*((f(x)*2*pi)^2)*M(s,s))+((1+1i*eta)*K(s,s))));
        FRF2_MSA(x) = FRF2_MSA(x) + ((cos(k_nat(s)*xL(2)))/(((-1)*((f(x)*2*pi)^2)*M(s,s))+((1+1i*eta)*K(s,s))));        
    end
end


% Plot Modal Superposition Approach
figure('Name', 'FRF Modal Superposition Approach');
subplot(2,2,1)
semilogy(f, abs(FRF1_MSA));
title('x = L/2');
xlabel('frequency [Hz]'); ylabel('FRF_x1 (f)');
grid on
subplot(2,2,2)
semilogy(f, abs(FRF2_MSA));
% axis([0, f_max, -0.00005, 0.00005]);
title('x = L/5');
xlabel('frequency [Hz]'); ylabel('FRF_x2 (f)');
grid on

subplot(2,2,3)
plot(f, angle(FRF1_MSA)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x1 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on
subplot(2,2,4)
plot(f, angle(FRF2_MSA)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x2 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on

% Comparison between SWS and MSA

figure('Name', 'Comparison between SWS and MSA ');
subplot(2,2,1)
semilogy(f, abs(FRF1_MSA), 'g', f, abs(FRF1_damp), 'm');
title('x = L/2');
xlabel('frequency [Hz]'); ylabel('FRF_x1 (f)');
grid on
subplot(2,2,2)
semilogy(f, abs(FRF2_MSA), 'g', f, abs(FRF2_damp), 'm');
% axis([0, f_max, -0.00005, 0.00005]);
title('x = L/5');
xlabel('frequency [Hz]'); ylabel('FRF_x2 (f)');
grid on

legend('MSA','SWS');

subplot(2,2,3)
plot(f, angle(FRF1_MSA)*180/pi,'b', f, angle(FRF1_damp)*180/pi,'r'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x1 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on
subplot(2,2,4)
plot(f, angle(FRF2_MSA)*180/pi,'b', f, angle(FRF2_damp)*180/pi,'r'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle FRF_x2 [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])

legend('MSA','SWS');

%% Case 4) 
%{ 
For the free-fixed bar excited at the free-end, plot and comment the driving-point
impedance in the following cases:
- damped bar (wave propagation solution) – loss factor ? = 0.01;
- damped bar (modal superposition approach) – loss factor ? = 0.01
%}

omega = 2*pi.*f;
x_resp = 0;

%% WPS
E_loss = E*(1+1i*eta);
k = sqrt(rho/E).*omega;
k_loss = sqrt(rho/E)*(1-1i*eta/2).*omega;
% k_loss = sqrt(rho/E_loss)*(1-1i*eta/2).*omega;
S = b*h;

%Z_wps_noloss = (-1i*sqrt(rho/E) * E * b * h).*(cos(sqrt(rho/E)*L.*omega)./sin((sqrt(rho/E).*omega * (L - x_resp))));
Z_wps = (-1i*sqrt(rho/E)*(1-1i*eta/2) * E_loss * S).*(cos(k_loss.*L)./sin(k_loss .* (L - x_resp)));

figure('Name', 'Case 4 - WPS')
subplot(2,1,1); 
semilogy(omega/2/pi, abs(Z_wps),'r'); 
grid off; xlabel('Frequency [Hz]'); ylabel('|Z| [Nm/s]'); ylim([0,10e7]); xlim([0,10e3]);

subplot(2,1,2); 
plot(omega/2/pi, angle(Z_wps)*180/pi,'r'); grid off; xlabel('Frequency [Hz]'); ylabel('\angle Z [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])
 
%% MSA
Nmodes = [5,8];
Z5_msa = 0;
Z8_msa = 0;
denZ5 =0;
denZ8 =0;

for ii = 1:Nmodes(1)
    denZ5 = denZ5 + 1./(-(omega.^2)*M(ii,ii)+ (1 + 1i*eta).*K(ii,ii));    
end
denZ5 = omega.*denZ5;
Z5_msa = (-1i)./denZ5;

for ii = 1:Nmodes(2)
    denZ8 = denZ8 + 1./(-(omega.^2)*M(ii,ii)+ (1 + 1i*eta).*K(ii,ii));    
end
denZ8 = omega.*denZ8;
Z8_msa = (-1i)./denZ8;

figure('Name', 'Case 4 - MSA (8 modes)')
subplot(2,1,1); 
semilogy(omega/2/pi, abs(Z8_msa),'r'); 
grid off; xlabel('Frequency [Hz]'); ylabel('|Z| [Nm/s]'); ylim([0,10e7]); xlim([0,10e3]);

subplot(2,1,2); 
plot(omega/2/pi, angle(Z8_msa)*180/pi,'r'); grid off; xlabel('Frequency [Hz]'); ylabel('\angle Z [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])

%% Comparison between 8 modes and 5 modes

figure('Name', 'Comparison between SWS, MSA 8 modes and MSA 5 modes ');
subplot(2,1,1); 
semilogy(omega/2/pi, abs(Z_wps),'r', omega/2/pi, abs(Z5_msa),'b-', omega/2/pi, abs(Z8_msa),'g-'); 
grid off; xlabel('Frequency [Hz]'); ylabel('|Z| [Nm/s]'); ylim([0,10e7]); xlim([0,10e3]);

subplot(2,1,2); 
plot(omega/2/pi, angle(Z_wps)*180/pi,'r', omega/2/pi, angle(Z5_msa)*180/pi,'b-', omega/2/pi, angle(Z8_msa)*180/pi,'g-'); 
grid off; xlabel('Frequency [Hz]'); ylabel('\angle Z [deg]'); yticks([-180 -90 0 90 180]); ylim([-200,200])

legend('SWS', 'MSA 5 modes', 'MSA 8 modes');


