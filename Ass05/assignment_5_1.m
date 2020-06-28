%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 03       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Importing experimental data
load("Data.mat");
load("Data2.mat");

% xy_bt back plate
% xy_bf top  plate

%% Filtraggio dati
% Fatto solo una volta e messo in Data2.mat

%{
for ii = 1:size(FRF,2)
    FRF_LP(:,ii) = lowpass(FRF(:,ii), 0.01);
end

save Data2 FRF;
%}
figure('Name', 'Una FRF a caso');
plot(freq, abs(FRF(:,45)), freq, abs(FRF_LP(:,45)));

%% Media
% Media FRF top plate
absFRF_top = mean(abs(FRF(:,1:58)),2);
% Media FRF bottom plate
absFRF_bottom = mean(abs(FRF(:,59:end)),2);
% Pulizia pulizia
f_cut = 0.03;
absFRF_top = lowpass(absFRF_top, f_cut);
absFRF_bottom = lowpass(absFRF_bottom, f_cut);


%{
figure('Name', 'FRF Mediate');
subplot(2,1,1);
plot(freq, abs(FRF_top));
subplot(2,1,2);
plot(freq, abs(FRF_bottom));
%}

%% Picchi - new
prom = 10;
width = 13;
[~, indices_t] = findpeaks(absFRF_top, 'MinPeakProminence', prom, 'MinPeakWidth', width);
[~, indices_b] = findpeaks(absFRF_bottom, 'MinPeakProminence', prom, 'MinPeakWidth', width);
%{
figure('Name', 'FRF Mediate');
subplot(2,1,1);
plot(freq, abs(FRF_top), freq(indices_t), abs(FRF_top(indices_t)), 'r*');

subplot(2,1,2);
plot(freq, abs(FRF_bottom), freq(indices_b), abs(FRF_bottom(indices_b)), 'r*');
%}
%% 2) Finding natural frequencies, damping ratios and mode shapes with simplified methods

% Natural frequencies
f_nat_t = freq(indices_t);
f_nat_b = freq(indices_b);

omega_nat_t = 2*pi.*f_nat_t;
omega_nat_b = 2*pi.*f_nat_b;

% find our frequency ranges (using local minima)
i_max = size(FRF, 1);

% Top Range
wp_dx = 1;
for yy = 1: (length(indices_t) + 1 ) % over the m(+1) peaks
    wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro
    
    if yy> length(indices_t) %(limite destro)
        wp_dx = i_max;
    else
        wp_dx = indices_t(yy); % w del picco successivo (limite destro)
    end
    [~, min_index] = min(absFRF_top(wp_sx:wp_dx));
    rfi_top(1,yy) = min_index + wp_sx - 1;
end
% Bottom range
wp_dx = 1;
for yy = 1: (length(indices_b) + 1 ) % over the m(+1) peaks
    wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro
    
    if yy> length(indices_b) %(limite destro)
        wp_dx = i_max;
    else
        wp_dx = indices_b(yy); % w del picco successivo (limite destro)
    end
    [~, min_index] = min(absFRF_bottom(wp_sx:wp_dx));
    rfi_bottom(1,yy) = min_index + wp_sx - 1;
end

clear yy wp_sx wp_dx min_index

figure('Name', 'FRF Mediate');
subplot(2,1,1);
plot(freq, absFRF_top, freq(indices_t), absFRF_top(indices_t), 'r*', ...
    freq(rfi_top), absFRF_top(rfi_top), 'g*');

subplot(2,1,2);
plot(freq, absFRF_bottom, freq(indices_b), absFRF_bottom(indices_b), 'r*', ...
    freq(rfi_bottom), absFRF_bottom(rfi_bottom), 'g*');

%vpar = zeros(size(indices,1), size(indices,2), 9);

% Minimizzazione
in1 = input('Load Minimization Data (1) or do a new minimization (2) ? \n');

if (in1 ==2)
    FRFreco = zeros(size(FRF));
    maxLen = max( length(indices_t), length(indices_b));
    vpar = NaN(maxLen,125,9);
    C = size(FRF,2);
    for mm = 1:C %over the n measurement
        
        if (C<=58)
            R = length(indices_t);
            for pp = 1:R % over the m peaks
                
                i_peak = indices_t(pp);
                w_peak = omega_nat_t(pp); % w del picco
                
                % Function reduced to the range of interest
                iini= rfi_top(pp); ifin = rfi_top(pp+1);
                rangeFreq = freq(iini:ifin);
                FRFexp = FRF(iini:ifin, mm);
                
                derPh = (angle(FRF(i_peak+1,mm))- angle(FRF(i_peak-1,mm)))/(2*pi*(freq(i_peak+1)-freq(i_peak-1)));
                csi0i = -1/(w_peak*derPh);
                csi(pp,mm) = csi0i;
                r0i = 2*w_peak*csi0i;
                
                % Constant parameter guessing (initial guess of Aj)
                Aj0i=-imag(FRF(i_peak,mm))*w_peak*r0i; % all the other constants are set equal to 0
                
                % Filling the vector xpar
                xpar0 = [csi0i ; w_peak; Aj0i; zeros(5,1)]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
                
                % Identification: single channel
                option=optimset('fminsearch');
                options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
                xpar=fminsearch(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, options);
                %xpar=fminsearch(@(xpar) errKjki_cw(xpar , rfHjki , Hjkiexp(:,jj)) , xpar0, options);
                
                % Plot results of identification
                vpar(pp,mm,:)=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
                %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
                
                
                % Reconstruction
                FRFreco(iini:ifin,mm) = funHjki(vpar(pp,mm,:), rangeFreq);
                %FRFreco(:,mm) = FRFreco(:,mm) + funHjki(vpar(pp,mm,:), freq);
            end
        else
            R = length(indices_b);
            for pp = 1:R % over the m peaks
                
                i_peak = indices_b(pp);
                w_peak = omega_nat_b(pp); % w del picco
                
                % Function reduced to the range of interest
                iini= rfi_bottom(pp); ifin = rfi_bottom(pp+1);
                rangeFreq = freq(iini:ifin);
                FRFexp = FRF(iini:ifin, mm);
                
                derPh = (angle(FRF(i_peak+1,mm))- angle(FRF(i_peak-1,mm)))/(2*pi*(freq(i_peak+1)-freq(i_peak-1)));
                csi0i = -1/(w_peak*derPh);
                csi(pp,mm) = csi0i;
                r0i = 2*w_peak*csi0i;
                
                % Constant parameter guessing (initial guess of Aj)
                Aj0i=-imag(FRF(i_peak,mm))*w_peak*r0i; % all the other constants are set equal to 0
                
                % Filling the vector xpar
                xpar0 = [csi0i ; w_peak; Aj0i; zeros(5,1)]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
                
                % Identification: single channel
                option=optimset('fminsearch');
                options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
                xpar=fminsearch(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, options);
                %xpar=fminsearch(@(xpar) errKjki_cw(xpar , rfHjki , Hjkiexp(:,jj)) , xpar0, options);
                
                % Plot results of identification
                vpar(pp,mm,:)=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
                %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
                
                
                % Reconstruction
                FRFreco(iini:ifin,mm) = funHjki(vpar(pp,mm,:), rangeFreq);
                %FRFreco(:,mm) = FRFreco(:,mm) + funHjki(vpar(pp,mm,:), freq);
            end
        end
    end
    save Data3 freq FRFreco vpar;
    
else
    load("Data3.mat");
end

% PLOT PROVA
absFRF_top2 = mean(abs(FRFreco(:,1:58)),2);
absFRF_bottom2 = mean(abs(FRFreco(:,59:end)),2);

figure('Name', 'FRF Mediate reconstruct');
subplot(2,1,1);
plot(freq, absFRF_top2, freq(indices_t), abs(absFRF_top2(indices_t)), 'r*', ...
    freq(rfi_top), absFRF_top2(rfi_top), 'g*');
ylim([0,500]);

subplot(2,1,2);
plot(freq, absFRF_bottom2, freq(indices_b), abs(absFRF_bottom2(indices_b)), 'r*', ...
    freq(rfi_bottom), absFRF_bottom2(rfi_bottom), 'g*');
ylim([0,500]);

%% Risultati ottenuti
% RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
vpar_t = vpar(1:length(indices_t),1:58,:);
vpar_b = vpar(1:length(indices_b),59:end ,:);

Xi_t = vpar_t(:,:,4);
Xi_b = vpar_b(:,:,4);

Xi_t = rescale(Xi_t, -1 , +1);
Xi_b = rescale(Xi_b, -1 , +1);

%% PARTE 3
modo = input('Which mode do you want to visualize ? ');

figure
%Top
subplot(1,2,1)
hold on
plot(xy_bt(:,1),xy_bt(:,2), '.b') % body %
plot(xy(1:58,1),xy(1:58,2), 'sr') % sensors %

%[X,Y] = meshgrid(-10:0.1:10,0:0.1:35);
[X,Y] = meshgrid(min(xy(1:58,1)):0.1:max(xy(1:58,1)), ...
    min(xy(1:58,2)):0.1:max(xy(1:58,2)));
b = boundary(xy(1:58,1), xy(1:58,2));
inmask = inpolygon(X(:), Y(:), xy(b,1), xy(b,2));
F_t = scatteredInterpolant(xy(1:58,1),xy(1:58,2), Xi_t(modo,:)');
Z_t = F_t(X,Y);
Z_t(~inmask) = NaN;
s = surf(X,Y,Z_t);
s.EdgeColor = 'none';
hold off 

%Bottom
subplot(1,2,2)
hold on
plot(xy_bf(:,1),xy_bf(:,2), '.b') % body %
plot(xy(58:125,1),xy(58:125,2), 'sr')  % sensors %

[X,Y] = meshgrid(min(xy(59:end,1)):0.1:max(xy(59:end,1)), ...
    min(xy(59:end,2)):0.1:max(xy(59:end,2)));
b = boundary(xy(1:58,1), xy(1:58,2));
inmask = inpolygon(X(:), Y(:), xy(b,1), xy(b,2));
F_b = scatteredInterpolant(xy(59:end,1),xy(59:end,2), Xi_b(modo,:)');
Z_b = F_b(X,Y);
Z_b(~inmask) = NaN;
s = surf(X,Y,Z_b);
s.EdgeColor  = 'none';

cb=colorbar;
cb.Position = cb.Position + [0.11, 0, 0, 0];
caxis([-1 +1]);

hold off









