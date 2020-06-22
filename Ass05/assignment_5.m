%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ASSIGNMENT 03       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
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
%}

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
    FRF(:,ii) = lowpass(FRF(:,ii), 0.01);
end

save Data2 FRF;
%}
figure();
plot(freq, abs(FRF(:,45)));

%% Picchi

% [peaks_1, indices_1] = findpeaks(abs(frf_1), 'MinPeakProminence', 15, 'MinPeakWidth', 15);
% plot(freq, abs(frf_1), freq(indices_1), abs(frf_1(indices_1)),'r*');

indices = zeros(30,size(FRF, 2));
for i=1:125
    frf_i = abs(FRF(:,i));
    [~, indices_i] = findpeaks(frf_i, 'MinPeakProminence', 15, 'MinPeakWidth', 4);
    indices(1:size(indices_i),i) = indices_i;
end


%% 2) Finding natural frequencies, damping ratios and mode shapes with simplified methods

% Natural frequencies
f_nat = zeros(size(indices));


for j = 1:size(indices, 2)
    for i = 1:size(indices, 1)
        if(indices(i,j) ~= 0 )
            f_nat(i,j) = freq(indices(i,j));
        else
            break;
        end
    end
end 

omega_nat = 2*pi.*f_nat;

% find our frequency ranges (local minima)
rf_index = zeros(size(indices, 1) + 1, size(indices,2));
i_max = size(FRF, 1);


for xx = 1:size(rf_index,2) %over the nn measurements
    wp_dx = 1;
    
    for yy = 1: (nnz(indices(:,xx)) +1) % over the m(+1) peaks
        wp_sx = wp_dx; % (limite sinistro) is vecchio limite destro
        
        if yy> nnz(indices(:,xx)) %(limite destro)
            wp_dx = i_max; 
        else
            wp_dx = indices(yy,xx); % w del picco successivo (limite destro)
        end
        [~, min_index] = min(abs(FRF((wp_sx:wp_dx), xx)));
        rf_index(yy,xx) = min_index + wp_sx-1;
    end
end

clear wp_sx wp_dx min_index

%vpar = zeros(size(indices,1), size(indices,2), 9);  

C = size(indices,2);
for mm = 1:C %over the n measurement
    R = nnz(indices(:,mm));
    for pp = 1:R % over the m peaks
        
        i_peak = indices(pp,mm); 
        w_peak = omega_nat(pp,mm); % w del picco
               
        % Function reduced to the range of interest
        iini= rf_index(pp,mm); ifin = rf_index(pp+1,mm);
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
        vpar=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
        %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]
      

        % Reconstruction
        FRFreco = funHjki(vpar, rangeFreq);
             
       
    end
end

%% PARTE 3












