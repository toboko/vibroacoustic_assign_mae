%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNDAMENTALS OF VIBRATION ANALYSIS AND VIBROACOUSTICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     MODULE 1 - FUNDAMENTALS OF VIBRATION ANALYSIS     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ASSIGNMENT 3 - MODAL PARAMETER IDENTIFICATION     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all


%% 1) Experimental FRFs

% Importing experimental data
data = load("Data.mat");

t  = data.Data(:,1)';
F  = data.Data(:,2)';

N_dof = size(data.Data,2) - 2;

x = zeros(N_dof, size(data.Data,1));

for i = 1:N_dof
    x(i,:) = data.Data(:,i+2)';
end

% FFT

N_up = 1; % Oversampling factor
N_fft = size(x,2) * N_up; % FFT points

F_fft = fft(F, N_fft) / length(F);

x_fft = zeros(size(x,1), N_fft);

for i = 1:size(x,1)
    x_fft(i,:) = fft(x(i,:), N_fft) / size(x,2);
end

% FRF

H_exp = zeros(size(x_fft,1), round(size(x_fft,2)/2));

for i = 1:size(H_exp,1)
    tmp = x_fft(i,:) ./ F_fft;
    H_exp(i,:) = tmp(1:length(H_exp(i,:))); % Considering half spectrum (till Nyquist frequency)
end

T_s = t(5) - t(4);
f_s = 1/T_s;
f_max = 5;
f = linspace(0, f_s/2, size(H_exp,2));
%i_max = find(abs(f-f_max) < 1e-4);
f(1) = 1e-4; % Preventing division by zero

% Plots

% Finding maximum value among all the plots (for y axis limits in magnitude
% plot)
y_max_H_exp = max(abs(H_exp(1,:)));

for i = 1:size(H_exp,1)
    y_max_H_exp = max([abs(H_exp(i,:)), y_max_H_exp]);
end

figure('Name', 'Experimental FRFs'' modulus');
sgtitle('Experimental FRFs'' modulus', 'FontSize', 20)

for i = 1:size(H_exp,1)
    subplot(ceil(size(H_exp,1)/2),2,i)
    plot(f, abs(H_exp(i,:)));
    axis([0, f_max, 0, y_max_H_exp*1.1]);
    title(sprintf('H^{exp}_%d(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]', 'FontSize', 12);
    ylabel(sprintf('|H^{exp}_%d(f)|   [m/N]', i), 'FontSize', 12);
    grid minor
end

figure('Name', 'Experimental FRFs'' phase');
sgtitle('Experimental FRFs'' phase', 'FontSize', 20)

for i = 1:size(H_exp,1)
    subplot(ceil(size(H_exp,1)/2),2,i)
    plot(f, angle(H_exp(i,:)));
    axis([0, f_max, -pi*1.1, pi*1.1]);
    title(sprintf('H^{exp}_%d(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]', 'FontSize', 12);
    ylabel(sprintf('$ \\angle\\bigl(H^\\textup{exp}_%d(f)\\bigr) $ \\, [deg]', i), 'Interpreter', 'LaTeX', 'FontSize', 12);
    set(gca, 'YTick', f);
    set(gca, 'YTickLabel', ["-\pi","-\pi/2","0","\pi/2","\pi"]);
    yticks([-pi, -pi/2, 0, pi/2 pi]);
    grid minor
end


%% 2) Finding natural frequencies, damping ratios and mode shapes with simplified methods

% Finding peaks and their indices

[peaks_H1, indices_H1] = findpeaks(abs(H_exp(1,:)));
[peaks_H2, indices_H2] = findpeaks(abs(H_exp(2,:)));
[peaks_H3, indices_H3] = findpeaks(abs(H_exp(3,:)));
[peaks_H4, indices_H4] = findpeaks(abs(H_exp(4,:)));

% Retaining only peaks at least 1 order of magnitude lower than maximum
% peak

indices_H1(peaks_H1 < max(peaks_H1)*1e-1) = [];
peaks_H1(peaks_H1 < max(peaks_H1)*1e-1) = [];
indices_H2(peaks_H2 < max(peaks_H2)*1e-1) = [];
peaks_H2(peaks_H2 < max(peaks_H2)*1e-1) = [];
indices_H3(peaks_H3 < max(peaks_H3)*1e-1) = [];
peaks_H3(peaks_H3 < max(peaks_H3)*1e-1) = [];
indices_H4(peaks_H4 < max(peaks_H4)*1e-1) = [];
peaks_H4(peaks_H4 < max(peaks_H4)*1e-1) = [];

% Each row corresponds to a single measurement
peaks_H = [peaks_H1; peaks_H2; peaks_H3; peaks_H4];
indices_H = [indices_H1; indices_H2; indices_H3; indices_H4];

% Natural frequencies
f_nat = [f(indices_H1); f(indices_H2); f(indices_H3); f(indices_H4);];
omega_nat = 2*pi*f_nat;

f_nat_mean = mean(f_nat,1);

% Applying the Half Power Bandwidth method, we estimated the factor having the 
% two frequencies f1, f2 in which the FRF halves the power

% The half power approach for estimating damping is based
% on finding the bandwidth for each mode

% Custom function to get row and column length
rows = @(x) size(x,1); 
cols = @(x) size(x,2);

% Init csi matrix
csi = zeros(rows(peaks_H), cols(peaks_H));

for x = 1:rows(csi)
    for y = 1:cols(csi)
        thresh_mag = peaks_H(x,y)*0.5*sqrt(2);
        i         = indices_H(x,y);
        f0        = f_nat(x,y);
        mag = abs(H_exp(x,i));
        % Search for f1
        while (mag > thresh_mag)
            %   Looking for index before f0
            i = i - 1;
            mag = abs(H_exp(x,i));
        end
        f1 = f(i + 1);
        
        % Reset index and magnitude under evaluation
        i  = indices_H(x,y);
        mag = abs(H_exp(x,i));
        
        % Search for f2
        while (mag > thresh_mag)
            %   Looking for index after f0
            i = i + 1;
            mag = abs(H_exp(x,i));
        end
        f2 = f(i - 1);
        
%       Calc csi
        csi(x,y) = (f2^2 - f1^2)/(4*f0^2);
    end
end

% Excluding values out of the mean
csi_red = csi;
csi_red(2,2) = NaN;
csi_red(3,3) = NaN;

csi_1 = [];
csi_2 = [];
csi_3 = [];
csi_4 = [];

for i = 1:N_dof
    if ~isnan(csi_red(i,1))
        csi_1 = [csi_1 csi_red(i,1)];
    end
end

for i = 1:N_dof
    if ~isnan(csi_red(i,2))
        csi_2 = [csi_2 csi_red(i,2)];
    end
end

for i = 1:N_dof
    if ~isnan(csi_red(i,3))
        csi_3 = [csi_3 csi_red(i,3)];
    end
end

for i = 1:N_dof
    if ~isnan(csi_red(i,4))
        csi_4 = [csi_4 csi_red(i,4)];
    end
end

csi_mean = [mean(csi_1), mean(csi_2), mean(csi_3), mean(csi_4)];

% Modal mass and modal damping matrices
m_q = ones(N_dof);
c_q = zeros(N_dof);

for k=1:N_dof
    for i=1:N_dof
        c_q(k,i) = 2 * csi(k,i) * m_q(k,i) * omega_nat(k,i);
    end
end

% Modal matrix
phi = zeros(N_dof);

for k=1:N_dof
    for i=1:N_dof
        phi(k,i) = -imag(H_exp(k,indices_H(k,i))) * omega_nat(k,i) * c_q(k,i);
    end
end

phi_norm = phi ./ phi(1,:);


%% OPTIONAL
%% i)

% find our frequency ranges (local minima)
rf = zeros(rows(indices_H), cols(indices_H)+1);
i_max = find(f == f_max);
wp_dx = 1;

for yy = 1:cols(rf) % over the m(+1) peaks
    wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro

    if yy>cols(indices_H) %(limite destro)
        wp_dx = i_max; %length(H1);
    else
        wp_dx = indices_H(1,yy); % w del picco successivo (limite destro)
    end
    [~, min_index] = min(H_exp(:,(wp_sx:wp_dx)),[],2);
    rf(:,yy) = min_index + wp_sx-1;
end

% Plots showing minima and maxima of each experimental FRF
figure('Name', 'Experimental FRFs'' modulus - Minima and maxima');
sgtitle(['Experimental FRFs'' modulus' newline 'Minima and maxima'], 'FontSize', 20)

for i = 1:size(H_exp,1)
    subplot(ceil(size(H_exp,1)/2),2,i)
    plot(f, abs(H_exp(i,:))); hold on
    plot(f(rf(i,:)), abs(H_exp(i,rf(i,:))), 'g*'); hold on
    plot(f(indices_H(i,:)), abs(H_exp(i,indices_H(i,:))), 'r*');
    axis([0, f_max, 0, y_max_H_exp*1.1]);
    title(sprintf('H^{exp}_%d(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]', 'FontSize', 12);
    ylabel(sprintf('|H^{exp}_%d(f)|   [m/N]', i), 'FontSize', 12);
    grid minor
end

clear wp_sx wp_dx min_index

% Preparing the figures for the identified FRFs
fig_1 = figure('Name', 'i) FRFs'' modulus - Analytical vs experimental');
fig_2 = figure('Name', 'i) FRFs'' phase - Analytical vs experimental');

R = rows(csi);
C = cols(csi);

vpar = zeros(R,C,9);
csi_min = zeros(R,C);
omega_nat_min = zeros(R,C);

disp('Minimization started...')

for mm = 1:R % Over the n measurements
    for pp = 1:C % Over the m peaks
        
        i_peak = indices_H(mm,pp); 
        w_peak = omega_nat(mm,pp); % w del picco
        csi0i = csi(mm,pp);
        r0i=2*w_peak*csi0i;
        
        % Function reduced to the range of interest
        iini= rf(mm,pp); ifin = rf(mm,pp+1);
        rfHjki = f(iini:ifin);
        Hjkiexp = H_exp(mm, iini:ifin); 
        
        % Constant parameter guessing (initial guess of Aj)
        Aj0i=-imag(H_exp(mm,i_peak))*w_peak*r0i; % all the other constants are set equal to 0

        % Filling the vector xpar (initial guess)
        xpar0 = [csi0i ; w_peak; Aj0i; zeros(5,1)]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
        
        % Identification: single channel
        option=optimset('fminsearch');
        options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
        options.Display = 'none';
        xpar=fminsearch(@(xpar) err_i(xpar , rfHjki , Hjkiexp), xpar0, options);

        % Plot results of identification
        vpar(mm,pp,:)=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
        %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
        csi_min(mm,pp) = xpar(1);
        omega_nat_min(mm,pp) = xpar(2);

        % Reconstruction
        Hjkiid=funHjki(vpar(mm,pp,:), rfHjki);
        
        %FRFreco(mm, iini:ifin) = Hjkiid; %Giustaapposteh
        
        % Plots
        
        % Finding maximum value among all the plots (for y axis limits in magnitude plot)
        y_max = max([abs(Hjkiexp), abs(Hjkiid)]);
        
        ngraph = (mm-1)*C + pp;
        
        figure(fig_1)
        subplot(R,C,ngraph);
        plot(rfHjki, abs(Hjkiid), 'r-', 'LineWidth', 1.2); hold on
        plot(rfHjki, abs(Hjkiexp), 'b');
        axis([-inf, inf, 0, y_max*1.1]);
        title(sprintf('H^{(%d)}_%d(f)', pp, mm), 'FontSize', 12);
        xlabel('frequency f [Hz]', 'FontSize', 10);
%         if pp == 1
% %             newrange = rfHjki;
% %             newrange(1) = 0;
% %             set(gca, 'XTick', rfHjki);
% %             set(gca, 'XTickLabel', newrange);
% %             set(gca, 'XTick', 1e-4);
% %             set(gca, 'XTickLabel', 0);
% %             xticks([1e-4]);
% %             xticklabels([0]);
%             set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
%         end
        ylabel(sprintf('|H^{(%d)}_%d(f)|   [m/N]', pp, mm), 'FontSize', 10);
        grid minor
        
        figure(fig_2)
        subplot(R,C,ngraph);
        plot(rfHjki, angle(Hjkiid)*180/pi, 'r-', 'LineWidth', 1.2); hold on
        plot(rfHjki, angle(Hjkiexp)*180/pi, 'b');
        axis([-inf, inf, -200, 200]);
        title(sprintf('H^{(%d)}_%d(f)', pp, mm), 'FontSize', 12);
        xlabel('frequency f [Hz]', 'FontSize', 10);
%         if pp == 1
% %             newrange = rfHjki;
% %             newrange(1) = 0;
% %             set(gca, 'XTick', rfHjki);
% %             set(gca, 'XTickLabel', newrange);
% %             set(gca, 'XTick', 1e-4);
% %             set(gca, 'XTickLabel', 0);
% %             xticks([1e-4]);
% %             xticklabels([0]);
%             set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
%         end
        ylabel(sprintf('$ \\angle\\bigl(H^{(%d)}_%d(f)\\bigr) $ \\, [deg]', pp, mm), 'Interpreter', 'LaTeX', 'FontSize', 10);
        yticks([-180, -90, 0, 90, 180]);
        grid minor
    end
    disp("Measurement " + num2str(mm) + " done." )
end

figure(fig_1)
sgtitle(['FRFs'' modulus' newline 'Analytical vs experimental'], 'FontSize', 16);
legend('Analytical', 'Experimental');

figure(fig_2)
sgtitle(['FRFs'' phase' newline 'Analytical vs experimental'], 'FontSize', 16);
legend('Analytical', 'Experimental');


%% ii) Modal parameter comparison

% Getting one value for each parameter (mean along measurements)
f_nat_min = omega_nat_min/(2*pi);

f_nat_min_mean = mean(f_nat_min,1);
csi_min_mean = mean(csi_min,1);

phi_min = vpar(:,:,4);
phi_min_norm = phi_min ./ phi_min(1,:);


%% iii) FRF reconstruction

% FRF reconstruction
H_rec = reco(vpar, f);

% Plots

% Finding maximum value among all the plots (for y axis limits in magnitude
% plot)
y_max_H_rec = max(abs(H_rec(1,:)));

for i = 1:size(H_exp,1)
    y_max_H_rec = max([abs(H_rec(i,:)), y_max_H_rec]);
end

figure('Name', 'iii) FRFs'' modulus - Reconstructed vs experimental');
sgtitle(['FRFs'' modulus' newline 'Reconstructed vs experimental'], 'FontSize', 20)

for i = 1:size(H_rec,1)
    subplot(ceil(size(H_rec,1)/2),2,i)
    plot(f, abs(H_rec(i,:)), f, abs(H_exp(i,:)));
    axis([0, f_max, 0, max([y_max_H_rec, y_max_H_exp])*1.1]);
    title(sprintf('H_%d(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]', 'FontSize', 12);
    ylabel(sprintf('|H_%d(f)|   [m/N]', i), 'FontSize', 12);
    grid minor
end

legend('Reconstructed', 'Experimental')

figure('Name', 'iii) FRFs'' phase - Reconstructed vs experimental');
sgtitle(['FRFs'' phase' newline 'Reconstructed vs experimental'], 'FontSize', 20)

for i = 1:size(H_rec,1)
    subplot(ceil(size(H_rec,1)/2),2,i)
    plot(f, angle(H_rec(i,:)), f, angle(H_exp(i,:)));
    axis([0, f_max, -pi*1.1, pi*1.1]);
    title(sprintf('H_%d(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]', 'FontSize', 12);
    ylabel(sprintf('$ \\angle\\bigl(H_%d(f)\\bigr) $ \\, [deg]', i), 'Interpreter', 'LaTeX', 'FontSize', 12);
    set(gca, 'YTick', f);
    set(gca, 'YTickLabel', ["-\pi","-\pi/2","0","\pi/2","\pi"]);
    yticks([-pi, -pi/2, 0, pi/2 pi]);
    grid minor
end

legend('Reconstructed', 'Experimental')


%% EXTRA
%% iv) Co-located FRF for x2

% Finding the modal FRM first

% Getting one value for each parameter (mean along measurements)
m_mean = mean(vpar(:,:,1),1);
c_mean = mean(vpar(:,:,2),1);
k_mean = mean(vpar(:,:,3),1);

M_q = diag(m_mean); % Modal mass matrix
C_q = diag(c_mean); % Modal damping matrix
K_q = diag(k_mean); % Modal stiffness matrix

% Modal FRM construction
omega = 2*pi*f;

FRM_q = zeros(N_dof, N_dof, length(omega)); % Modal Frequency Response Matrix has one modal FRF (vector in the 3rd dimension) in each of the Ndof^2 spots in the first 2 dimensions
D_q = zeros(N_dof, N_dof, length(omega)); % Modal mechanical impedance matrix

for i = 1:N_dof
    for j = 1:N_dof
        D_q(i,j,:) = -omega.^2 * M_q(i,j) + 1i*omega * C_q(i,j) + K_q(i,j);
    end
end

for i = 1:length(omega)
    FRM_q(:,:,i) = inv(D_q(:,:,i));
end

% Plots

% Finding maximum value among all the plots (for y axis limits in magnitude plot)
y_max = max(abs(FRM_q(1,1,:)));

for i = 1:size(FRM_q, 1)
    for j = 1:size(FRM_q, 2)
        y_max = max([max(abs(FRM_q(i,j,:))), y_max]);
    end
end

figure('Name', 'iv) Modal FRFs'' modulus')
sgtitle('Modal FRFs'' modulus', 'FontSize', 20);

for i = 1:size(FRM_q, 1) % Rows: loop through modal variables
    for j = 1:size(FRM_q, 2) % Columns: loop through modal forces
        k = j + size(FRM_q, 1) * (i-1); % Retrieve loop counter from row and column indices
        subplot(N_dof,N_dof,k);
        plot(f, abs(squeeze(FRM_q(i,j,:))), 'b');
        axis([0, f_max, 0, y_max*1.1]);
        title(sprintf('H_{q_{%d%d}}(f)', i, j), 'FontSize', 16);
        xlabel('frequency f [Hz]', 'FontSize', 12);
        ylabel(sprintf('|H_{q_{%d%d}}(f)|   [m/N]', i, j), 'FontSize', 12);
        grid minor;
    end
end

figure('Name', 'iv) Modal FRFs'' phase')
sgtitle('Modal FRFs'' phase', 'FontSize', 20);

for i = 1:size(FRM_q, 1) % Rows: loop through modal variables
    for j = 1:size(FRM_q, 2) % Columns: loop through modal forces
        k = j + size(FRM_q, 1) * (i-1); % Retrieve loop counter from row and column indices
        subplot(N_dof,N_dof,k);
        plot(f, angle(squeeze(FRM_q(i,j,:)))*180/pi, 'b');
        axis([0, f_max, -200, 200]);
        title(sprintf('H_{q_{%d%d}}(f)', i, j), 'FontSize', 16);
        xlabel('frequency f [Hz]', 'FontSize', 12);
        ylabel(sprintf('$ \\angle\\bigl(H_{q_{%d%d}}(f)\\bigr) $ \\, [deg]', i, j), 'Interpreter', 'LaTeX', 'FontSize', 12);
        yticks([-180, -90, 0, 90, 180]);
        grid minor;
    end
end

% Finding the Frequency Response Matrix of numerical model from modal FRM
FRM = zeros(N_dof, N_dof, length(omega));

for i = 1:length(omega)
    FRM(:,:,i) = phi_min_norm * FRM_q(:,:,i) * phi_min_norm';
end

% Plots

% Finding maximum value among all the plots (for y axis limits in magnitude plot)
y_max = max(abs(FRM(1,1,:)));

for i = 1:size(FRM, 1)
    for j = 1:size(FRM, 2)
        y_max = max([max(abs(FRM(i,j,:))), y_max]);
    end
end

% Plots

figure('Name', 'iv) Modulus of system''s FRFs')
sgtitle('Modulus of system''s FRFs', 'FontSize', 20);

for i = 1:size(FRM, 1) % Rows: loop through independent variables
    for j = 1:size(FRM, 2) % Columns: loop through generalized forces
        k = j + size(FRM, 1) * (i-1); % Retrieve loop counter from row and column indices
        subplot(N_dof,N_dof,k);
        plot(f, abs(squeeze(FRM(i,j,:))), 'b');
        axis([0, f_max, 0, y_max*1.1]);
        title(sprintf('H_{%d%d}(f)', i, j), 'FontSize', 16);
        xlabel('frequency f [Hz]', 'FontSize', 12);
        ylabel(sprintf('|H_{%d%d}(f)|   [m/N]', i, j), 'FontSize', 12);
        grid minor;
    end
end

figure('Name', 'iv) Phase of system''s FRFs')
sgtitle('Phase of system''s FRFs', 'FontSize', 20);

for i = 1:size(FRM, 1) % Rows: loop through independent variables
    for j = 1:size(FRM, 2) % Columns: loop through generalized forces
        k = j + size(FRM, 1) * (i-1); % Retrieve loop counter from row and column indices
        subplot(N_dof,N_dof,k);
        plot(f, angle(squeeze(FRM(i,j,:)))*180/pi, 'b');
        axis([0, f_max, -200, 200]);
        title(sprintf('H_{%d%d}(f)', i, j), 'FontSize', 16);
        xlabel('frequency f [Hz]', 'FontSize', 12);
        ylabel(sprintf('$ \\angle\\bigl(H_{%d%d}(f)\\bigr) $ \\, [deg]', i, j), 'Interpreter', 'LaTeX', 'FontSize', 12);
        yticks([-180, -90, 0, 90, 180]);
        grid minor;
    end
end

% Comparison between first column of system's FRM and experimental FRFs

% Finding maximum value among all the plots (for y axis limits in magnitude plot)
y_max_H = max(abs(FRM(1,1,:)));

for i = 1:size(FRM, 1)
    y_max_H = max([max(abs(FRM(i,1,:))), y_max_H]);
end

figure('Name', 'iv) FRFs'' modulus - Numerical model vs experimental')
sgtitle(['FRFs'' modulus' newline 'Numerical model vs experimental'], 'FontSize', 20);
for i = 1:size(FRM, 1) % Rows: loop through independent variables
    subplot(N_dof,1,i);
    plot(f, abs(squeeze(FRM(i,1,:))), 'b', f, abs(H_exp(i,:)), 'r:', 'LineWidth', 1.2);
    axis([0, f_max, 0, max([y_max_H, y_max_H_exp])*1.1]);
    title(sprintf('H_{%d1}(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]','FontSize',12);
    ylabel(sprintf('|H_{%d1}(f)|   [m/N]', i), 'FontSize', 12);
    grid minor;
end

legend('Numerical model', 'Experimental');

figure('Name', 'iv) FRFs'' phase - Numerical model vs experimental')
sgtitle(['FRFs'' phase' newline 'Numerical model vs experimental'], 'FontSize', 20);

for i = 1:size(FRM, 1) % Rows: loop through independent variables
    subplot(N_dof,1,i);
    plot(f, angle(squeeze(FRM(i,1,:)))*180/pi, 'b', f, angle(H_exp(i,:))*180/pi, 'r:', 'LineWidth', 1.2);
    axis([0, f_max, -200, 200]);
    title(sprintf('H_{%d1}(f)', i), 'FontSize', 16);
    xlabel('frequency f [Hz]', 'FontSize', 12);
    ylabel(sprintf('$ \\angle\\bigl(H_{%d1}(f)\\bigr) $ \\, [deg]', i), 'Interpreter', 'LaTeX', 'FontSize', 12);
    yticks([-180, -90, 0, 90, 180]);
    grid minor;
end

legend('Numerical model', 'Experimental');

% Co-located FRF for point 2 is the second element on the FRM's diagonal

figure('Name', 'iv) Co-located FRF of point 2')
sgtitle('Co-located FRF of point 2', 'FontSize', 20);

subplot(2,1,1);
plot(f, abs(squeeze(FRM(2,2,:))));
axis([0, f_max, 0, max(abs(FRM(2,2,:)))*1.1]);
title('Modulus', 'FontSize', 16);
xlabel('frequency f [Hz]', 'FontSize', 12);
ylabel('|H^{coloc}_2(f)|   [m/N]', 'FontSize', 12);
grid minor;

subplot(2,1,2);
plot(f, angle(squeeze(FRM(2,2,:)))*180/pi);
axis([0, f_max, -200, 200]);
title('Phase', 'FontSize', 16);
xlabel('frequency f [Hz]', 'FontSize', 12);
ylabel('$ \angle\bigl(H^{coloc}_2(f)\bigr) $ \, [deg]', 'Interpreter', 'LaTeX', 'FontSize', 12);
yticks([-180, -90, 0, 90, 180]);
grid minor;


%% v) Generalized matrices of numerical model

% Modes are broadly orthogonal: M* = (phi^T)^(-1) M_q phi^(-1)

M_gen = phi_min_norm' \ M_q / phi_min_norm;
C_gen = phi_min_norm' \ C_q / phi_min_norm;
K_gen = phi_min_norm' \ K_q / phi_min_norm;

% Proof that the FRM built from generalized matrices is the same as the one
% obtained from modal FRM

% System's FRM construction
FRM_prime = zeros(N_dof, N_dof, length(omega)); % Frequency Response Matrix has one FRF (vector in the 3rd dimension) in each of the N_dof^2 spots in the first 2 dimensions
D = zeros(N_dof, N_dof, length(omega)); % Mechanical impedance matrix

for i = 1:N_dof
    for j = 1:N_dof
        D(i,j,:) = -omega.^2 * M_gen(i,j) + 1i*omega * C_gen(i,j) + K_gen(i,j);
    end
end

for i = 1:length(omega)
    FRM_prime(:,:,i) = inv(D(:,:,i));
end

% Plots

% Finding maximum value among all the plots (for y axis limits in magnitude plot)
y_max = max(abs(FRM_prime(1,1,:)));

for i = 1:size(FRM_prime, 1)
    for j = 1:size(FRM_prime, 2)
        y_max = max([max(abs(FRM_prime(i,j,:))), y_max]);
    end
end

figure('Name', 'Modulus of system''s FRFs')
sgtitle('Modulus of system''s FRFs', 'FontSize', 20);

for i = 1:size(FRM_prime, 1) % Rows: loop through independent variables
    for j = 1:size(FRM_prime, 2) % Columns: loop through generalized forces
        k = j + size(FRM_prime, 1) * (i-1); % Retrieve loop counter from row and column indices
        subplot(N_dof,N_dof,k);
        plot(f, abs(squeeze(FRM_prime(i,j,:))), 'b');
        axis([0, f_max, 0, y_max*1.1]);
        title(sprintf('H_{%d%d}(f)', i, j), 'FontSize', 16);
        xlabel('frequency f [Hz]','FontSize',12);
        ylabel(sprintf('|H_{%d%d}(f)|   [m/N]', i, j), 'FontSize', 12);
        grid minor;
    end
end

figure('Name', 'Phase of system''s FRFs')
sgtitle('Phase of system''s FRFs', 'FontSize', 20);

for i = 1:size(FRM_prime, 1) % Rows: loop through independent variables
    for j = 1:size(FRM_prime, 2) % Columns: loop through generalized forces
        k = j + size(FRM_prime, 1) * (i-1); % Retrieve loop counter from row and column indices
        subplot(N_dof,N_dof,k);
        plot(f, angle(squeeze(FRM_prime(i,j,:)))*180/pi, 'b');
        axis([0, f_max, -200, 200]);
        title(sprintf('H_{%d%d}(f)', i, j), 'FontSize', 16);
        xlabel('frequency f [Hz]', 'FontSize', 12);
        ylabel(sprintf('$ \\angle\\bigl(H_{%d%d}(f)\\bigr) $ \\, [deg]', i, j), 'Interpreter', 'LaTeX', 'FontSize', 12);
        yticks([-180, -90, 0, 90, 180]);
        grid minor;
    end
end