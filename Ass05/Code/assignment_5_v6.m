%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    FUNDAMENTALS OF VIBRATION ANALYSIS AND VIBROACOUSTICS     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MODULE 2 - VIBROACOUSTICS OF MUSICAL INSTRUMENTS       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ASSIGNMENT 2 - EXPERIMENTAL MODAL ANALYSIS OF A VIOLIN    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all


% Importing experimental data
load("Data.mat");
load("Data2.mat");
disp("Data loaded.");

% xy_bt back plate
% xy_bf top  plate


%% FRF scala
omega = 2*pi.*freq + 1e-10;
FRF = FRF ./ repmat((-(omega.^2)), 1, size(FRF,2));
FRF_LP = FRF_LP ./ repmat((-(omega.^2)), 1, size(FRF_LP,2));

%% Filtraggio dati
% Fatto solo una volta e messo in Data2.mat
%{
for ii = 1:size(FRF,2)
    FRF_LP(:,ii) = lowpass(FRF(:,ii), 0.01);

    asdf = unwrap(angle(FRF(:,ii)));
    minasdf = min(asdf,[],1);
    angleFRF(:,ii) = lowpass(asdf + minasdf, 0.02);
    angleFRF(:,ii) = angleFRF(:,ii) - minasdf;
end

save Data2 freq FRF_LP angleFRF;
%}
figure('Name', 'Una FRF a caso passabassata');
plot(freq, abs(FRF(:,17)), freq, abs(FRF_LP(:,17)));
ylim([0, 1e-04])

%% Media
% Media FRF 
FRFmean = mean(abs(FRF),2);
% Pulizia pulizia
f_cut = 0.4;
FRFmean = lowpass(FRFmean, f_cut);

%% Picchi - new
prom = 0.000003;
width = 8;
tr = 0;
iini = 500;
[~, indices] = findpeaks(FRFmean(iini:end), 'MinPeakProminence', prom, 'MinPeakWidth', width, 'Threshold', tr);
Nmodes = length(indices);
indices = indices + iini;
disp(num2str(Nmodes) + " modes found.");

%% 2) Finding natural frequencies et ranges
% Natural frequencies
f_nat = freq(indices);

omega_nat = 2*pi.*f_nat;

% find our frequency ranges (using local minima)
i_max = size(FRF, 1);

% Range
wp_dx = iini;
for yy = 1: (Nmodes + 1 ) % over the m(+1) peaks
    wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro
    
    if yy> Nmodes %(limite destro)
        wp_dx = i_max;
    else
        wp_dx = indices(yy); % w del picco successivo (limite destro)
    end
    %[~, min_index] = min(FRFmean(wp_sx:wp_dx));
    [~, min_index] = min(abs(FRF_LP(wp_sx:wp_dx, :)), [], 1);
    rfi(:,yy) = min_index + wp_sx - 1;
end

clear yy wp_sx wp_dx min_index

% figure('Name', 'FRF Mediate');
% plot(freq, FRFmean, freq(indices), FRFmean(indices), 'r*', ...
%     freq(rfi), FRFmean(rfi), 'g*');
% ylim([0, 1e-04])

%vpar = zeros(size(indices,1), size(indices,2), 9);

% Plot showing minima and maxima of mean FRF
figure('Name', 'Mean FRF''s modulus')
plot(freq, abs(FRFmean)); hold on
plot(freq(indices), abs(FRFmean(indices)), 'r*');
axis([-inf, inf, 0, 1e-4]);
title('Mean FRF''s modulus', 'FontSize', 16);
xlabel('frequency f [Hz]', 'FontSize', 12);
ylabel('|G^{mean}(f)|   [m/N]', 'FontSize', 12);
grid minor

%% Minimizzazione
in1 = input('Load Minimization Data (1) or do a new minimization (2) ? \n');

if (in1 == 2)
    disp("Minimization started."); tic
    FRFreco = zeros(size(FRF));
    % RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
    vpar = NaN(Nmodes,125,4);
    C = 125;
    i_p = zeros(Nmodes,125);
    derPhmatrix = zeros(Nmodes,125);

for mm = 1:C       % over the n measurements
        R = Nmodes;

        for pp = 1:R % over the p peaks            
            i_peak = indices(pp);
            
            % Function reduced to the range of interest
            iini= rfi(mm,pp); ifin = rfi(mm,pp+1);
            rangeFreq = freq(iini:ifin);
            FRFexp = FRF_LP(iini:ifin, mm);
            omegaRange = 2*pi.*rangeFreq;
            
            % New peak
            safe = 25;
%             i1 = i_peak - safe;
%             i2 = min(i_peak + safe, 5002);
            i1 = iini;
            i2 = ifin;
            [~, i_peak_tmp] = max(abs(FRF_LP(i1:i2,mm)),[],1);

            if (i_peak_tmp == 1 || i_peak_tmp == i2-i1)
                i_peak_tmp = indices(pp);
            else
                i_peak_tmp = i_peak_tmp + i1;
            end

            
            i_peak = i_peak_tmp;
            f_nat_tmp = freq(i_peak);
            w_peak = 2*pi.*f_nat_tmp;
            i_p(pp,mm) = i_peak;
            
            % Phase derivative method (ma con fase smussata)
            safe = 3;
            derPh = (angleFRF(i_peak+safe,mm) - angleFRF(i_peak-safe,mm)) ...
                             ./(2*pi*(freq(i_peak+safe)-freq(i_peak-safe)));
            %derPh = (angle(FRF(i_peak+safe,mm))- angle(FRF(i_peak-safe,mm)))./(2*pi*(freq(i_peak+safe)-freq(i_peak-safe)));
%              if (derPh > -3e-02 )
%                  derPh = -1e03;
%              end
            if (derPh < 1e03)
                derPh = -1e00;
            end
            derPhmatrix(pp,mm) = derPh;

            csi0i = -1./(w_peak.*derPh);
            r0i = 2.*w_peak.*csi0i;
            r0matrix(pp,mm) = r0i;
            
            % Constant parameter guessing (initial guess of Aj)
            Aj0i=-imag(FRF(i_peak,mm)).*w_peak.*r0i; % all the other constants are set equal to 0
            
            % Filling the vector xpar
            %  xpar(par,mis) - 3x125
            xpar0 = [csi0i ; w_peak; Aj0i; zeros(5,1)]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
            
            % Identification
            
%             option=optimset('fminsearch');
%             options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
%             xpar=fminsearch(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, options);
             
            options = optimoptions(@lsqnonlin,'TolFun', 1e-10, 'TolX', 1e-10, 'StepTolerance', 1e-20, ...
                'MaxIter', 1e06, 'MaxFunEval', 1e07);
            %options.Algorithm = 'levenberg-marquardt';
            options.Display = 'none';
            safe = 25;
%             lb = [1e-05; omegaRange(1); -Inf.*ones(6,1)];
%             ub = [0.2; omegaRange(end); Inf.*ones(6,1) ];
            lb = [0; w_peak-safe ; -50*ones(6,1)];
            ub = [0.1; w_peak+safe ; 50*ones(6,1) ];
%             lb = [];
%             ub = [];
            xpar=lsqnonlin(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, lb, ub, options);
            % Plot results of identification
            % RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
            csi(pp,mm,1) = xpar(1); % csi
            csi(pp,mm,2) = xpar(2); % nat freq
            vpar(pp,mm,1) = 1; %m
            vpar(pp,mm,2) = 2*xpar(1).*xpar(2); %c
            vpar(pp,mm,3) = xpar(2).^2; %k
            vpar(pp,mm,4) = xpar(3); %A
            vpar(pp,mm,5) = xpar(4); %B
            vpar(pp,mm,6) = xpar(5); %C
            vpar(pp,mm,7) = xpar(6); %D
            vpar(pp,mm,8) = xpar(7); %E
            vpar(pp,mm,9) = xpar(8); %F
            %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
            
            % Reconstruction
            %FRFreco(iini:ifin,mm) = funHjki(vpar(pp,mm,:), rangeFreq);
            %disp("Peak " + num2str(pp) + " done.");
        end
        disp("Measurement " + num2str(mm) + " done.");
end

    
    disp("Minimization ended.");
    save Data3 freq FRFreco vpar csi i_p;
    toc
    
else
    disp("Loading minimization d a t a . . .");
    load("Data3.mat");
end

% Reconstruction
FRFreco = reco(vpar, freq);

% PLOT PROVA
% absFRFreco = mean(abs(FRFreco),2);
% 
% figure('Name', 'FRF Mediate reconstruct');
% 
% plot(freq, absFRFreco, freq(indices), abs(absFRFreco(indices)), 'r*', ...
%     freq(rfi), absFRFreco(rfi), 'g*');
% ylim([0, 1e-04])

%% Risultati ottenuti
% RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
vpar_t = vpar(:,1:58,:);
vpar_b = vpar(:,59:end ,:);

Xi_t = (vpar_t(:,:,4));
Xi_b = (vpar_b(:,:,4));

%% 3b comparison of experimental and identified FRFs (for a certain reference channel);

% ch = input('Which channel do you want to visualize ? \n');
% while(ch< 1 || ch > 125)
%    disp("Insert a valid channel index (<= 125 )\n");
%    ch = input('Which channel do you want to visualize ? \n'); 
% end

figure('Name', 'Comparison between identified and experimental FRFs')

for mode = 1:Nmodes
    for ch = 1:size(FRF,2)
        hold off
        plot(freq(rfi(ch,mode):rfi(ch,mode+1)), abs(FRF((rfi(ch,mode):rfi(ch,mode+1)),ch)), 'r'); hold on
        plot(freq(rfi(ch,mode):rfi(ch,mode+1)), abs(FRFreco((rfi(ch,mode):rfi(ch,mode+1)),ch)), 'Linewidth', 0.9, 'Color', 'g');
        axis([-inf, inf, 0, max([max(abs(FRF((rfi(ch,mode):rfi(ch,mode+1)),ch))), max(abs(FRFreco((rfi(ch,mode):rfi(ch,mode+1)),ch)))])*1.1]);
        title(sprintf(['Comparison between identified and experimental FRFs' newline 'Mode %d - Channel %d'], mode, ch), 'FontSize', 16);
        xlabel('frequency f [Hz]', 'FontSize', 12);
        ylabel(sprintf('|H_{%d}^{%d}(f)|   [m/N]', ch, mode), 'FontSize', 12);
        legend('Experimental', 'Reconstructed');
        grid minor;

        resp = input('Press ''s'' to skip to next mode, press any other key to go on \n', 's');
        if (strcmp(resp, 's'))
            break
        end
    end
end

figure('Name', 'Comparison between identified and experimental FRFs')

plot(freq, abs(FRF(:,ch)), 'r'); hold on
plot(freq, abs(FRF_LP(:,ch)), 'Color', 'b'); hold on
plot(freq, abs(FRFreco(:,ch)), 'Linewidth', 0.9, 'Color', 'g'); hold on
plot(freq(i_p(:,ch)), abs(FRF(i_p(:,ch),ch)), 'k*');
axis([-inf, inf, 0, 1e-04]);
title(['Comparison between identified and experimental FRFs' newline 'Channel 45'], 'FontSize', 16);
xlabel('frequency f [Hz]', 'FontSize', 12);
ylabel('|H_{45}(f)|   [m/N]', 'FontSize', 12);
legend('Experimental', 'Lowpassed experimental', 'Reconstructed', 'Mode peaks');
grid minor;


%% 3c visualization of the identified modes.

% modo = input('Which mode do you want to visualize ? \n');
% while(modo< 1 || modo > Nmodes)
%    disp("Insert a valid mode index (<= " + num2str(Nmodes) + " ) ");
%    modo = input('Which mode do you want to visualize ? \n'); 
% end
% disp("Natural frequency: " + num2str(f_nat(modo)) + " Hz");

csi_mean = mean(csi(:,:,1),2);

for i = 1:Nmodes

figure
sgtitle(sprintf('Mode %d, natural frequency %.1f Hz, damping ratio %.2f%%', i, f_nat(i), csi_mean(i)*100), 'FontSize', 16);
%Top
subplot(1,2,1)

F_t = scatteredInterpolant(xy(1:58,1),xy(1:58,2), Xi_t(i,:)', ...
    'natural', 'none');
[X_t,Y_t] = meshgrid(min(xy(1:58,1)):0.05:max(xy(1:58,1)), ...
    min(xy(1:58,2)):0.05:max(xy(1:58,2)));
Z_t = F_t(X_t,Y_t);
b = boundary(xy(1:58,1), xy(1:58,2));
inmask = inpolygon(X_t(:), Y_t(:), xy(b,1), xy(b,2));
Z_t(~inmask) = NaN;
Z_t = rescale(Z_t, -1 , +1);

s = pcolor(X_t,Y_t,Z_t);
s.LineStyle = 'none';


hold on
plot(xy_bt(:,1),xy_bt(:,2), '.b') % body 
plot(xy(1:58,1),xy(1:58,2), 'sr') % sensors 
xlim([-10, +10]); ylim([0, 35]);

caxis([-1 +1]);
hold off

title('Front plate', 'FontSize', 12);

%Bottom
subplot(1,2,2)

F_b = scatteredInterpolant(xy(59:end,1),xy(59:end,2), Xi_b(i,:)', ...
    'natural', 'none');
[X_b,Y_b] = meshgrid(min(xy(59:end,1)):0.05:max(xy(59:end,1)), ...
    min(xy(59:end,2)):0.05:max(xy(59:end,2)));
Z_b = F_b(X_b,Y_b);
b = boundary(xy(1:58,1), xy(1:58,2));
inmask = inpolygon(X_b(:), Y_b(:), xy(b,1), xy(b,2));
Z_b(~inmask) = NaN;
Z_b = rescale(Z_b, -1 , +1);

s = pcolor(X_b,Y_b,Z_b);
s.LineStyle = 'none';

hold on
plot(xy_bf(:,1),xy_bf(:,2), '.b') % body %
plot(xy(59:end,1),xy(59:end,2), 'sr')  % sensors %
xlim([-10, +10]); ylim([0, 35]);

caxis([-1 +1]);

colormap jet;
cb=colorbar;
cb.Position = cb.Position + [0.11, 0, 0, 0];
caxis([-1 +1]);

hold off

title('Back plate', 'FontSize', 12);

end

%% Auxiliar plots and other stuff
csis = [0.6, 0.6, 0.8, 0, 0, 1.2, 0; 0.62, 0.58, 0.86, 0.59, 1.26, 1.21, 1.68; [0.00911731052116327, 0.00833520312205319, 0.00873127832851640, 0.00642083648408587, 0.0164981132874108, 0.0130281237120128, 0.0161938724587724] * 100; csi_mean' * 100]';