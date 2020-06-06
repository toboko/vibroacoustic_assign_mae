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
E = 70e09;
L = 2;
h = 0.05;
b = h;
f_max = 10000;

% Utils
T = 10;
f_s = 2 * f_max;
lambda_s = 0.005;
% %omega_max = 2*pi*f_max;
% omega_s = 2*pi*f_s;
% c = sqrt(E/rho);
% %k_max = omega_max/c;
% k_s = omega_s/c;
% %k_s = 2 * k_max;
% lambda_s = 2*pi/k_s;

t = 0:1/f_s:T;
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
title('1^{st} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(1)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,2)
plot(x, real(phi_free_fixed(2,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('2^{nd} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(2)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,3)
plot(x, real(phi_free_fixed(3,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('3^{rd} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(3)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,4)
plot(x, real(phi_free_fixed(4,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('4^{th} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(4)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,5)
plot(x, real(phi_free_fixed(5,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('5^{th} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(5)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,6)
plot(x, real(phi_free_fixed(6,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('6^{th} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(6)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,7)
plot(x, real(phi_free_fixed(7,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('7^{th} mode shape - free-fixed bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fx}^{(7)}}(x)','FontSize',12);
grid on

subplot(i_free_fixed_max/2,2,8)
plot(x, real(phi_free_fixed(8,:)));
axis([-inf, inf, min([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1, max([real(phi_free_fixed(1,:)), real(phi_free_fixed(2,:)), real(phi_free_fixed(3,:)), real(phi_free_fixed(4,:)), real(phi_free_fixed(5,:)), real(phi_free_fixed(6,:)), real(phi_free_fixed(7,:)), real(phi_free_fixed(8,:))])*1.1]);
title('8^{th} mode shape - free-fixed bar','FontSize',16);
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
title('1^{st} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(1)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,2)
plot(x, real(phi_free_free(2,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('2^{nd} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(2)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,3)
plot(x, real(phi_free_free(3,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('3^{rd} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(3)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,4)
plot(x, real(phi_free_free(4,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('4^{th} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(4)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,5)
plot(x, real(phi_free_free(5,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('5^{th} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(5)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,6)
plot(x, real(phi_free_free(6,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('6^{th} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(6)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,7)
plot(x, real(phi_free_free(7,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('7^{th} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(7)}}(x)','FontSize',12);
grid on

subplot(ceil(i_free_free_max/2),2,8)
plot(x, real(phi_free_free(8,:)));
axis([-inf, inf, min([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1, max([real(phi_free_free(1,:)), real(phi_free_free(2,:)), real(phi_free_free(3,:)), real(phi_free_free(4,:)), real(phi_free_free(5,:)), real(phi_free_free(6,:)), real(phi_free_free(7,:)), real(phi_free_free(8,:))])*1.1]);
title('8^{th} mode shape - free-free bar','FontSize',16);
xlabel('position x [m]','FontSize',12); ylabel('\Phi^{{fr-fr}^{(8)}}(x)','FontSize',12);
grid on
