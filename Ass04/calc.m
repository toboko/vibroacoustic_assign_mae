%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 04       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% Custom function to get row and col length
rows = @(x) size(x,1); 
cols = @(x) size(x,2);

% Importing experimental data
p   = 2700;
E   = 70e9;
L   = 2;
h   = 0.05;
xL1 = L/2;
xL2 = L/5;
S   = h * h;

% Undamped bar (Standing Wave Solution)
T_s   = 0.05;
f_s   = 1/T_s;
f_max = 10000;

ll    = 10000;
f     = linspace(0, f_max, ll);
FRF   = zeros(size(f));

for x = 1:ll
    k      = 2 * pi * f(x) * sqrt(p ./ E);
    FRF(x) = sin(k * (L - xL1)) ./ (k * E * S * cos(k * L));
end

% Plots
figure('Name', 'FRF Output');
subplot(4,1,1)
semilogy(f, abs(FRF));
% axis([0, f_max, -0.00005, 0.00005]);
title('Undamped Bar');
xlabel('frequency [Hz]'); ylabel('FRF(f)');
grid on

subplot(4,1,2)
plot(f, angle(FRF)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle H [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on

% Damped bar (Standing Wave Solution)
f     = linspace(0.1, f_max, ll);
FRF   = zeros(size(f));
n     = 0.01;
E     = E * (1 + 1i * n);

for x = 1:ll
    k      = (2 * pi* f(x) * sqrt(p ./ E)) * (1 - 1i * (n ./ 2));
    FRF(x) = sin(k * (L - xL1)) ./ (k * E * S * cos(k * L));
end

% Plots
subplot(4,1,3)
semilogy(f, abs(FRF));
% axis([0, f_max, -0.00005, 0.00005]);
title('Damped Bar');
xlabel('frequency [Hz]'); ylabel('FRF(f)');

subplot(4,1,4)
plot(f, -angle(FRF)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle H [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on

% Damped bar (Modal Superposition Approach)
f     = linspace(0.1, f_max, ll);
FRF   = zeros(size(f));
E     = 70e9;

i_free_fixed_max = floor(2*L * f_max * sqrt(p/E) + 1/2);
% Natural frequencies
f_nat_free_fixed = zeros(1,i_free_fixed_max);

for i=1:i_free_fixed_max
    f_nat_free_fixed(i) = (2*i-1)/(4*L) * sqrt(E/p);
end

m = p * S * L;
M = zeros(i_free_fixed_max, i_free_fixed_max);
K = zeros(i_free_fixed_max, i_free_fixed_max);

K_nat = sqrt(p / E) * 2 * pi .* f_nat_free_fixed;


for x = 1:i_free_fixed_max
    for y = 1:i_free_fixed_max
        kr = K_nat(x);
        ks = K_nat(y);
        if (x == y)
            M(x,y) = m * (L/2 + (((sin(2 * ks * L)) / (4 * ks))));
        else
            M(x,y) = m * ((sin(L*(kr - ks))/(2*(kr - ks))) + (sin(L*(kr + ks))/(2*(kr + ks))));
        end
    end
end

for x = 1:i_free_fixed_max
    for y = 1:i_free_fixed_max
        kr = K_nat(x);
        ks = K_nat(y);
        if (x == y)
            K(x,y) = (E * S * kr * ks) * (L/2 - ((( sin(2 * ks * L)) / (4 * ks))));
        else
            K(x,y) = (E * S * kr * ks) * ((sin(L*(kr - ks))/(2*(kr - ks))) - (sin(L*(kr + ks))/(2*(kr + ks))));
        end
    end
end

FRF_MSA = zeros(size(f));

for x = 1:ll
    for s = 1:i_free_fixed_max
        FRF_MSA(x) = FRF_MSA(x) + (cos(K_nat(s)*xL1)/((-(f(x)*2*pi)*M(s,s))+((1+1i*n)*K(s,s))));
    end
end

% Plots
figure('Name', 'FRF Modal');
subplot(2,1,1)
semilogy(f, abs(FRF_MSA));
% axis([0, f_max, -0.00005, 0.00005]);
title('Undamped Bar');
xlabel('frequency [Hz]'); ylabel('FRF(f)');
grid on

subplot(2,1,2)
plot(f, angle(FRF_MSA)*180/pi,'b'); 
grid minor; 
xlabel('Frequency [Hz]'); 
ylabel('\angle H [deg]'); 
yticks([-180 -90 0 90 180]); 
ylim([-200,200])
grid on

