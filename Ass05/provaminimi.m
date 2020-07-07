clc
clear
close all

% Importing experimental data
load("Data.mat");
load("Data2.mat");
disp("Data loaded.");

%% Media
% Media FRF 
FRFmean = mean(abs(FRF),2);
% Pulizia pulizia
f_cut = 0.03;
FRFmean = lowpass(FRFmean, f_cut);

% Picchi
prom = 30;
width = 10;
[~, indices] = findpeaks(FRFmean, 'MinPeakProminence', prom, 'MinPeakWidth', width);
Nmodes = length(indices);
disp(num2str(Nmodes) + " modes found.");

% Natural frequencies
f_nat = freq(indices);

omega_nat = 2*pi.*f_nat;

% find our frequency ranges (using local minima)
i_max = size(FRF, 1);

%% Range
wp_dx = 1;
for yy = 1: (Nmodes + 1 ) % over the m(+1) peaks
    wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro
    
    if yy> Nmodes %(limite destro)
        wp_dx = i_max;
    else
        wp_dx = indices(yy); % w del picco successivo (limite destro)
    end
    [~, min_index] = min(FRFmean(wp_sx:wp_dx));
    %min_index = round((wp_dx-wp_sx)/2);
    rfi(1,yy) = min_index + wp_sx - 1;
end

clear yy wp_sx wp_dx min_index

figure('Name', 'FRF Mediate');
plot(freq, FRFmean, freq(indices), FRFmean(indices), 'r*', ...
    freq(rfi), FRFmean(rfi), 'g*');

R = Nmodes;
i_p = zeros(Nmodes,125);
for pp = 1:R % over the m peaks
    
    i_peak = indices(pp).*ones(1,125);
    i_peak_tmp = zeros(1,125);
    %             w_peak = omega_nat(pp); % w del picco
    
    % Function reduced to the range of interest
    iini= rfi(pp); ifin = rfi(pp+1);
    rangeFreq = freq(iini:ifin).*ones(1,125);
    FRFexp = FRF(iini:ifin, :);
    
    % New peak
%     i1 = iini + round(abs(i_peak -iini).*0.9);
%     i2 = ifin - round(abs(ifin - i_peak).*0.1);
    safe = 25;
    i1 = i_peak - safe;
    i2 = min(i_peak + safe, 5002);
    [~, i_peak_tmp] = max(abs(FRF_LP(i1:i2,:)),[],1);
    %[~, i_peak_tmp] = max(abs(FRF_LP(iini:ifin,:)),[],1);
    for xx = 1:size(i_peak,2)
        if (i_peak_tmp(xx) ==i1(xx) || i_peak_tmp(xx) ==i2(xx))
            i_peak_tmp(xx) = indices(pp);
        else
            i_peak_tmp(xx) = i_peak_tmp(xx) + i1(xx);
        end
    end
    i_peak = i_peak_tmp;
    f_nat_tmp = freq(i_peak_tmp)';
    w_peak = 2*pi.*f_nat_tmp;
    i_p(pp,:) = i_peak(1,:);
%     derPh = (angle(FRF(i_peak+1))- angle(FRF(i_peak-1)))./(2*pi*(freq(i_peak+1)'-freq(i_peak-1)'));
%     csi0i = -1./(w_peak.*derPh);
%     r0i = 2.*w_peak.*csi0i;
%     
%     % Constant parameter guessing (initial guess of Aj)
%     Aj0i=-imag(FRF(i_peak)).*w_peak.*r0i; % all the other constants are set equal to 0
%     
%     % Filling the vector xpar
%     xpar0 = [csi0i ; w_peak; Aj0i]; % all the unknowns (2 mod parameters + 6 constant (all the other constants are set equal to 0))
%     
%     % Identification: single channel
%     
%     %             option=optimset('fminsearch');
%     %             options=optimset(option, 'TolFun', 1e-8, 'TolX', 1e-8);
%     %             xpar=fminsearch(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, options);
%     
%     options = optimoptions(@lsqnonlin,'TolFun', 1e-8, 'TolX', 1e-8);
%     options.Algorithm = 'levenberg-marquardt';
%     options.Display = 'none';
%     xpar=lsqnonlin(@(xpar) err_i(xpar , rangeFreq , FRFexp), xpar0, [], [], options);
%     % Plot results of identification
%     vpar(:,:,pp)= [ones(1,125,1); 2.*xpar(1,:).*xpar(2,:); xpar(2,:).^2; xpar(3,:)];
%     csi(:,:,pp) = [xpar(1,:); xpar(2,:)];
%     %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
%     
%     % Reconstruction
%     FRFreco(iini:ifin,:) = funHjki(vpar(:,:,pp), rangeFreq);
%     %FRFreco = FRFreco + funHjki(vpar(:,:,pp), freq.*ones(1,125));
%     disp("Peak " + num2str(pp) + " done.");
end

    
%% 3b comparison of experimental and identified FRFs (for a certain reference channel);
ch = input('Which channel do you want to visualize ? \n');
while(ch< 1 || ch > 125)
   disp("Insert a valid channel index (<= 125 )\n");
   ch = input('Which channel do you want to visualize ? \n'); 
end

figure('Name', 'Comparison of experimental and identified FRFs')

plot(freq, abs(FRF(:,ch)), 'y');
hold on
plot(freq, abs(FRF_LP(:,ch)), 'b');
hold on
plot(freq(i_p(:,ch)), abs(FRF(i_p(:,ch),ch)), 'r*');
plot(freq(indices), abs(FRF(indices,ch)), 'k*');
hold on
plot(freq(rfi), abs(FRF(rfi,ch)), 'g*');

