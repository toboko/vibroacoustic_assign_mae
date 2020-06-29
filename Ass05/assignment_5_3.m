%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 05       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Importing experimental data
load("Data.mat");
load("Data2.mat");
disp("Data loaded.");

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
figure('Name', 'Una FRF a caso passabassata');
plot(freq, abs(FRF(:,45)), freq, abs(FRF_LP(:,45)));

%% Media
% Media FRF 
FRFmean = mean(abs(FRF),2);
% Pulizia pulizia
f_cut = 0.03;
FRFmean = lowpass(FRFmean, f_cut);

%{
figure('Name', 'FRF Mediate');
subplot(2,1,1);
plot(freq, abs(FRF_top));
subplot(2,1,2);
plot(freq, abs(FRF_bottom));
%}

%% Picchi - new
prom = 30;
width = 10;
[~, indices] = findpeaks(FRFmean, 'MinPeakProminence', prom, 'MinPeakWidth', width);
Nmodes = length(indices);
disp(num2str(Nmodes) + " modes found.");

%% 2) Finding natural frequencies, damping ratios and mode shapes with simplified methods

% Natural frequencies
f_nat = freq(indices);

omega_nat = 2*pi.*f_nat;

% find our frequency ranges (using local minima)
i_max = size(FRF, 1);

% Range
wp_dx = 1;
for yy = 1: (Nmodes + 1 ) % over the m(+1) peaks
    wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro
    
    if yy> Nmodes %(limite destro)
        wp_dx = i_max;
    else
        wp_dx = indices(yy); % w del picco successivo (limite destro)
    end
    %[~, min_index] = min(FRFmean(wp_sx:wp_dx));
    min_index = round((wp_dx-wp_sx)/2);
    rfi(1,yy) = min_index + wp_sx - 1;
end

clear yy wp_sx wp_dx min_index

figure('Name', 'FRF Mediate');
plot(freq, FRFmean, freq(indices), FRFmean(indices), 'r*', ...
    freq(rfi), FRFmean(rfi), 'g*');

%vpar = zeros(size(indices,1), size(indices,2), 9);

%% Minimizzazione
in1 = input('Load Minimization Data (1) or do a new minimization (2) ? \n');

if (in1 ==2)
    disp("Minimization started."); tic
    FRFreco = zeros(size(FRF));
    maxLen = Nmodes;
    vpar = NaN(9,125,maxLen);
    %C = size(FRF,2);
    
    %for mm = 1:C %over the n measurement
        
        R = Nmodes;
        for pp = 1:R % over the m peaks
            
            i_peak = indices(pp);
            w_peak = omega_nat(pp); % w del picco
            
            % Function reduced to the range of interest
            iini= rfi(pp); ifin = rfi(pp+1);
            rangeFreq = freq(iini:ifin).*ones(1,125);
            FRFexp = FRF(iini:ifin, :);
            
            derPh = (angle(FRF(i_peak+1,:))- angle(FRF(i_peak-1,:)))/(2*pi*(freq(i_peak+1)-freq(i_peak-1)));
            csi0i = -1./(w_peak.*derPh);
            r0i = 2.*w_peak.*csi0i;
            
            % Constant parameter guessing (initial guess of Aj)
            Aj0i=-imag(FRF(i_peak,:)).*w_peak.*r0i; % all the other constants are set equal to 0
            
            % Filling the vector xpar
            xpar0 = [csi0i ; w_peak.*ones(1,125); Aj0i; zeros(5,125)]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
            
            % Identification: single channel
            
%             option=optimset('fminsearch');
%             options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
%             xpar=fminsearch(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, options);
             
            options = optimoptions(@lsqnonlin,'TolFun', 1e-8, 'TolX', 1e-8);
            options.Algorithm = 'levenberg-marquardt';
            options.Display = 'none';
            xpar=lsqnonlin(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, [], [], options);
            % Plot results of identification
            vpar(:,:,pp)= [ones(1,125,1); 2.*xpar(1,:).*xpar(2,:); xpar(2,:).^2; xpar(3:8,:)];
            csi(:,:,pp) = [xpar(1,:); xpar(2,:)];
            %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
            
            % Reconstruction
            FRFreco(iini:ifin,:) = funHjki(vpar(:,:,pp), rangeFreq);
            %FRFreco = FRFreco + funHjki(vpar(:,:,pp), freq.*ones(1,125));
            disp("Peak " + num2str(pp) + " done.");
        end
    %end
    
    disp("Minimization ended.");
    save Data3 freq FRFreco vpar csi;
    toc
    
else
    disp("Loading minimization d a t a . . .");
    load("Data3.mat");
end

% PLOT PROVA
absFRFreco = mean(abs(FRFreco),2);

figure('Name', 'FRF Mediate reconstruct');

plot(freq, absFRFreco, freq(indices), abs(absFRFreco(indices)), 'r*', ...
    freq(rfi), absFRFreco(rfi), 'g*');
ylim([0,500]);

%% Risultati ottenuti
% RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
vpar_t = vpar(:,1:58,:);
vpar_b = vpar(:,59:end ,:);

Xi_t = vpar_t(:,:,4) + 1i.*vpar_t(:,:,5);
Xi_b = vpar_b(:,:,4) + 1i.*vpar_b(:,:,5);

% Xi_t = rescale(abs(Xi_t), -1 , +1);
% Xi_b = rescale(abs(Xi_b), -1 , +1);

%% 3b comparison of experimental and identified FRFs (for a certain reference channel);
ch = input('Which channel do you want to visualize ? \n');
while(ch< 1 || ch > 125)
   disp("Insert a valid channel index (<= 125 )\n");
   ch = input('Which channel do you want to visualize ? \n'); 
end

figure('Name', 'Comparison of experimental and identified FRFs')

plot(freq, abs(FRF(:,ch)));
hold on
plot(freq, abs(FRFreco(:,ch)), 'Linewidth', 2, 'Color', 'g');
hold on
plot(freq(indices), abs(FRF(indices,ch)), 'r*');

%% 3c visualization of the identified modes.

modo = input('Which mode do you want to visualize ? \n');
while(modo< 1 || modo > Nmodes)
   disp("Insert a valid mode index (<= " + num2str(Nmodes) + " ) ");
   modo = input('Which mode do you want to visualize ? \n'); 
end
disp("Natural frequency: " + num2str(f_nat(modo)) + " Hz");
%disp("Damping ratio: " + num2str(csiAv()) + "\n");

figure
%Top
subplot(1,2,1)

[X_t,Y_t] = meshgrid(min(xy(1:58,1)):0.1:max(xy(1:58,1)), ...
    min(xy(1:58,2)):0.1:max(xy(1:58,2)));
b = boundary(xy(1:58,1), xy(1:58,2));
inmask = inpolygon(X_t(:), Y_t(:), xy(b,1), xy(b,2));
F_t = scatteredInterpolant(xy(1:58,1),xy(1:58,2), Xi_t(modo,:)', ...
    'natural', 'none');
Z_t = F_t(X_t,Y_t);
Z_t = rescale(abs(Z_t), -1 , +1);
Z_t(~inmask) = NaN;

%s = surf(X_t,Y_t,Z_t);
s = pcolor(X_t,Y_t,Z_t);
s.LineStyle = 'none';

hold on
plot(xy_bt(:,1),xy_bt(:,2), '.b') % body 
plot(xy(1:58,1),xy(1:58,2), 'sr') % sensors 
xlim([-10, +10]); ylim([0, 35]);

hold off 

%Bottom
subplot(1,2,2)

[X_b,Y_b] = meshgrid(min(xy(59:end,1)):0.1:max(xy(59:end,1)), ...
    min(xy(59:end,2)):0.1:max(xy(59:end,2)));
b = boundary(xy(1:58,1), xy(1:58,2));
inmask = inpolygon(X_b(:), Y_b(:), xy(b,1), xy(b,2));
F_b = scatteredInterpolant(xy(59:end,1),xy(59:end,2), Xi_b(modo,:)', ...
    'natural', 'none');
Z_b = F_b(X_b,Y_b);
Z_b = rescale(abs(Z_b)  , -1 , +1);
Z_b(~inmask) = NaN;

%s = surf(X_b,Y_b,Z_b);
s = pcolor(X_b,Y_b,Z_b);
s.LineStyle = 'none';

hold on
plot(xy_bf(:,1),xy_bf(:,2), '.b') % body %
plot(xy(59:end,1),xy(59:end,2), 'sr')  % sensors %
xlim([-10, +10]); ylim([0, 35]);

cb=colorbar;
cb.Position = cb.Position + [0.11, 0, 0, 0];
caxis([-1 +1]);

hold off









