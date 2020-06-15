%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 03       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function count = counter(value, matrix)
    count = 0;
    for i = 1:size(matrix,1)
        for j = 1:size(matrix,2)
            if matrix(i,j) == value
                count = count + 1;
            end
        end
    end
end

clc
clear
close all

% Importing experimental data
load("Data.mat");

% xy_bt back plate
% xy_bf top  plate

% [peaks_1, indices_1] = findpeaks(abs(frf_1), 'MinPeakProminence', 15, 'MinPeakWidth', 15);
% plot(freq, abs(frf_1), freq(indices_1), abs(frf_1(indices_1)),'r*');

s_indices = zeros(20, size(FRF, 2));

for i=59:125
    frf_i = abs(FRF(:,i));
    [peaks_i, indices_i] = findpeaks(frf_i, 'MinPeakProminence', 15, 'MinPeakWidth', 4);
    s_indices(1:size(indices_i),i) = indices_i;
end

counter = 0;
verified_modes = 1;

for i=1:size(FRF, 2)
    for j=1:nnz(s_indices(:,i))
        chess = s_indices(j,i);
        
        % se l'elemento non è stato verificato
        if ismember(chess, verified_modes,'legacy') == 0
            % controlla se il valore è ripetuto
            for z=1:size(FRF, 2)
                for k=1:nnz(s_indices(:,i))
                    if(s_indices(k,z) ==  chess)
                        counter = counter + 1;
                    end
                end
            end
            
            if counter >= 4
                verified_modes = [verified_modes chess];
            end
        end
        counter = 0;
    end
end

frf_1 = FRF(:,112);
% plot(freq, abs(frf_1), freq(verified_modes), abs(frf_1(verified_modes)),'r*'); 

counters = zeros(size(freq));

for i=1:length(freq)
    counters(i) = sum(s_indices(:) == i);
end

[peaks, indices] = findpeaks(counters, 'MinPeakProminence', 4, 'MinPeakWidth', 1);

figure;
subplot(1, 2, 1); plot(freq, counters, freq(indices), counters(indices), 'r*');
subplot(1, 2, 2); plot(freq, abs(frf_1), freq(indices), abs(frf_1(indices)), 'r*'); 

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

        % Filling the vector xpar
        xpar0 = [csi0i ; w_peak; Aj0i; zeros(5,1)]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
        
        % Identification: single channel
        option=optimset('fminsearch');
        options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
        xpar=fminsearch(@(xpar) err_i(xpar , rfHjki , Hjkiexp), xpar0, options);
        %xpar=fminsearch(@(xpar) errKjki_cw(xpar , rfHjki , Hjkiexp(:,jj)) , xpar0, options);

        % Plot results of identification
        vpar=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
        %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]

        % Reconstruction
        Hjkiid=funHjki(vpar, rfHjki);
        
        %Plot - plis, do not explode
        ngraph = (mm-1)*C + pp;
        figure(f1)
        subplot(R,C,ngraph);
        plot(rfHjki, abs(Hjkiexp), 'b', rfHjki, abs(Hjkiid), 'r-', 'linewidth',1.2)
        xlabel('Frequency [Hz]')
        ylabel(['|H' num2str(mm) '| [m/N]'])
        
        figure(f2)
        subplot(R,C,ngraph);
        %plot(rfHjki, angle(Hjkiexp)*180/pi, 'b', rfHjki, angle(Hjkiid)*180/pi, 'r-', 'linewidth',1.2)
        plot(rfHjki, unwrap(angle(Hjkiexp)*180/pi), 'b', rfHjki, unwrap(angle(Hjkiid)*180/pi), 'r-', 'linewidth',1.2)
        yticks([-180 -90 0 90 180]); ylim([-200,200]);
        ylabel(['\angleH' num2str(mm) ' [rad]'])
        xlabel('Frequency [Hz]')
        
       
    end
end

figure(f1)
title('Magnitude')
legend('Experimental','Identified')
figure(f2)
title('Phase')
legend('Experimental','Identified')


