%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 03       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% 1) Experimental FRFs

% Importing experimental data
st = load("Data.mat");

t  = st.Data(:,1)';
F  = st.Data(:,2)';
x1 = st.Data(:,3)';
x2 = st.Data(:,4)';
x3 = st.Data(:,5)';
x4 = st.Data(:,6)';

% FFT
fft_F  = fft(F);
fft_x1 = fft(x1);
fft_x2 = fft(x2);
fft_x3 = fft(x3);
fft_x4 = fft(x4);

% FRF
H1 =  fft_x1./fft_F;
H2 =  fft_x2./fft_F;
H3 =  fft_x3./fft_F;
H4 =  fft_x4./fft_F;

% Considering half spectrum (till Nyquist frequency)
H1 = H1(1:round(length(H1)/2));
H2 = H2(1:round(length(H2)/2));
H3 = H3(1:round(length(H3)/2));
H4 = H4(1:round(length(H4)/2));

H  = [H1; H2; H3; H4];

n = 4;
T_s = t(5) - t(4);
f_s = 1/T_s;
f_max = 5;
f = linspace(0,f_s/2,length(H1));

% Plots

figure('Name', 'Experimental FRFs'' modulus');

subplot(2,2,1)
plot(f,abs(H1));
axis([0, f_max, 0, max(abs(H1))*1.1]);
title('FRF between x_1 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_1(f)|','FontSize',12);
grid on

subplot(2,2,2)
plot(f,abs(H2));
axis([0, f_max, 0, max(abs(H2))*1.1]);
title('FRF between x_2 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_2(f)|','FontSize',12);
grid on

subplot(2,2,3)
plot(f,abs(H3));
axis([0, f_max, 0, max(abs(H3))*1.1]);
title('FRF between x_3 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_3(f)|','FontSize',12);
grid on

subplot(2,2,4)
plot(f,abs(H4));
axis([0, f_max, 0, max(abs(H4))*1.1]);
title('FRF between x_4 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_4(f)|','FontSize',12);
grid on


figure('Name', 'Experimental FRFs'' phase');

subplot(2,2,1)
plot(f,angle(H1));
axis([0, f_max, -pi, pi]);
title('FRF between x_1 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_1(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

subplot(2,2,2)
plot(f,angle(H2));
axis([0, f_max, -pi, pi]);
title('FRF between x_2 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_2(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

subplot(2,2,3)
plot(f,angle(H3));
axis([0, f_max, -pi, pi]);
title('FRF between x_3 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_3(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

subplot(2,2,4)
plot(f,angle(H4));
axis([0, f_max, -pi, pi]);
title('FRF between x_4 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_4(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on


%% 2) Finding natural frequencies, damping ratios and mode shapes with simplified methods

% Finding peaks and their indices
[peaks_H1, indices_H1] = findpeaks(abs(H1));
[peaks_H2, indices_H2] = findpeaks(abs(H2));
[peaks_H3, indices_H3] = findpeaks(abs(H3));
[peaks_H4, indices_H4] = findpeaks(abs(H4));

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
peaks_H   = [peaks_H1; peaks_H2; peaks_H3; peaks_H4];
indices_H = [indices_H1; indices_H2; indices_H3; indices_H4];

% Natural frequencies
f_nat = [f(indices_H1); f(indices_H2); f(indices_H3); f(indices_H4);];
omega_nat = 2*pi.*f_nat;

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
        mag = abs(H(x,i));
        % Search for f1
        while (mag > thresh_mag)
            %   Looking for index before f0
            i = i - 1;
            mag = abs(H(x,i));
        end
        f1 = f(i + 1);
        
        % Reset index and magnitude under evaluation
        i  = indices_H(x,y);
        mag = abs(H(x,i));
        
        % Search for f2
        while (mag > thresh_mag)
            %   Looking for index after f0
            i = i + 1;
            mag = abs(H(x,i));
        end
        f2 = f(i - 1);
        
%       Calc csi
        csi(x,y) = (f2^2 - f1^2)/(4*f0^2);
    end
end

% Modal mass and modal damping matrices
m_q = ones(n);
c_q = zeros(n);

for k=1:n
    for i=1:n
        c_q(k,i) = 2 * csi(k,i) * m_q(k,i) * omega_nat(k,i);
    end
end

% Modal matrix
phi = zeros(n);

for k=1:n
    for i=1:n
        H_exp = H(k,:);
        phi(k,i) = -imag(H_exp(indices_H(k,i))) * omega_nat(k,i) * c_q(k,i);
    end
end

phi_norm = phi./phi(1,:);

%% Opzionale

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
        [~, min_index] = min(H(:,(wp_sx:wp_dx)),[],2);
        rf(:,yy) = min_index + wp_sx-1;
    end

clear wp_sx wp_dx min_index

f1 = figure('Name','MAGNITUDE - Just dont explode plis');
f2 = figure('Name','PHASE - Just dont explode plis');
    
R = rows(csi);
C = cols(csi);

vpar = zeros(R,C,9);
FRFreco = zeros(size(H));

disp('Minimization started...')
for mm = 1:R %over the n measurement
    for pp = 1:C % over the m peaks
        
        i_peak = indices_H(mm,pp); 
        w_peak = omega_nat(mm,pp); % w del picco
        csi0i = csi(mm,pp);
        r0i=2*w_peak*csi0i;
        
        % Function reduced to the range of interest
        iini= rf(mm,pp); ifin = rf(mm,pp+1);
        rfHjki = f(iini:ifin);
        Hjkiexp = H(mm, iini:ifin); 
        
        % Constant parameter guessing (initial guess of Aj)
        Aj0i=-imag(H(mm,i_peak))*w_peak*r0i; % all the other constants are set equal to 0

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

        % Reconstruction
        Hjkiid=funHjki(vpar(mm,pp,:), rfHjki);
        
        %FRFreco(mm, iini:ifin) = Hjkiid; %Giustaapposteh

        
        %Plot - plis, do not explode
        ngraph = (mm-1)*C + pp;
        figure(f1)
        subplot(R,C,ngraph);
        plot(rfHjki, abs(Hjkiexp), 'b', rfHjki, abs(Hjkiid), 'r-', 'linewidth',1.2)
        xlabel('Frequency [Hz]')
        ylabel(['|H' num2str(mm) '| [m/N]'])
        
        figure(f2)
        subplot(R,C,ngraph);
        plot(rfHjki, angle(Hjkiexp)*180/pi, 'b', rfHjki, angle(Hjkiid)*180/pi, 'r-', 'linewidth',1.2)
        %plot(rfHjki, unwrap(angle(Hjkiexp)*180/pi), 'b', rfHjki, unwrap(angle(Hjkiid)*180/pi), 'r-', 'linewidth',1.2)
        yticks([-180 -90 0 90 180]); ylim([-200,200]);
        ylabel(['\angleH' num2str(mm) ' [rad]'])
        xlabel('Frequency [Hz]')
        
       
    end
    disp("Measurement " + num2str(mm) + " done." )
end


figure(f1)
title('Magnitude')
legend('Experimental','Identified')
figure(f2)
title('Phase')
legend('Experimental','Identified')

%% Altri plot

% Plots MAX & MIN
figure('Name', 'Experimental FRFs'' modulus - Max & Min');

subplot(2,2,1)
plot(f,abs(H1));
hold on
plot(f(rf(1,:)),abs(H1(rf(1,:))), 'g*');
hold on
plot(f(indices_H(1,:)),abs(H1(indices_H(1,:))), 'r*');
axis([0, f_max, 0, max(abs(H1))*1.1]);
title('FRF between x_1 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_1(f)|','FontSize',12);
grid on

subplot(2,2,2)
plot(f,abs(H2));
hold on
plot(f(rf(2,:)),abs(H2(rf(2,:))), 'g*');
hold on
plot(f(indices_H(2,:)),abs(H2(indices_H(2,:))), 'r*');
axis([0, f_max, 0, max(abs(H2))*1.1]);
title('FRF between x_2 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_2(f)|','FontSize',12);
grid on

subplot(2,2,3)
plot(f,abs(H3));
hold on
plot(f(rf(3,:)),abs(H3(rf(3,:))), 'g*');
hold on
plot(f(indices_H(3,:)),abs(H3(indices_H(3,:))), 'r*');
axis([0, f_max, 0, max(abs(H3))*1.1]);
title('FRF between x_3 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_3(f)|','FontSize',12);
grid on

subplot(2,2,4)
plot(f,abs(H4));
hold on
plot(f(rf(4,:)),abs(H4(rf(4,:))), 'g*');
hold on
plot(f(indices_H(4,:)),abs(H4(indices_H(4,:))), 'r*');
axis([0, f_max, 0, max(abs(H4))*1.1]);
title('FRF between x_4 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_4(f)|','FontSize',12);
grid on

% Plots FRF reco
FRFreco = reco(vpar, f);

figure('Name', 'Reconstructed FRFs'' modulus');

subplot(2,2,1)
plot(f,abs(H1));
hold on
plot(f,abs(FRFreco(1,:)));
axis([0, f_max, 0, max(abs(FRFreco(1,:)))*1.1]);
title('FRF between x_1 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_1(f)|','FontSize',12);
grid on

subplot(2,2,2)
plot(f,abs(H2));
hold on
plot(f,abs(FRFreco(2,:)));
axis([0, f_max, 0, max(abs(FRFreco(2,:)))*1.1]);
title('FRF between x_2 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_2(f)|','FontSize',12);
grid on

subplot(2,2,3)
plot(f,abs(H3));
hold on
plot(f,abs(FRFreco(3,:)));
axis([0, f_max, 0, max(abs(FRFreco(3,:)))*1.1]);
title('FRF between x_3 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_3(f)|','FontSize',12);
grid on

subplot(2,2,4)
plot(f,abs(H4));
hold on
plot(f,abs(FRFreco(4,:)));
axis([0, f_max, 0, max(abs(FRFreco(4,:)))*1.1]);
title('FRF between x_4 and F - modulus','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('|H^{exp}_4(f)|','FontSize',12);
grid on


legend('Experimental', 'Reconstructed')

figure('Name', 'Recontructed FRFs'' phase');

subplot(2,2,1)
plot(f,angle(H1));
hold on
plot(f,angle(FRFreco(1,:)));
axis([0, f_max, -pi, pi]);
title('FRF between x_1 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_1(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

subplot(2,2,2)
plot(f,angle(H2));
hold on
plot(f,angle(FRFreco(2,:)));
axis([0, f_max, -pi, pi]);
title('FRF between x_2 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_2(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

subplot(2,2,3)
plot(f,angle(H3));
hold on
plot(f,angle(FRFreco(3,:)));
axis([0, f_max, -pi, pi]);
title('FRF between x_3 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_3(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

subplot(2,2,4)
plot(f,angle(H4));
hold on
plot(f,angle(FRFreco(4,:)));
axis([0, f_max, -pi, pi]);
title('FRF between x_4 and F - phase','FontSize',8);
xlabel('frequency f [Hz]','FontSize',12); ylabel('$ \angle\bigl(H^\textup{exp}_4(f)\bigr) $ \, [rad]','Interpreter','LaTeX','FontSize',12);
grid on

legend('Experimental', 'Reconstructed')


