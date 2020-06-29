function errore = err_i(xpar , range_H , H_exp) 
%ERR_I Summary of this function goes here
%   Detailed explanation goes here

if (size(xpar,1)~= 8 )
    error('Not valid xpar');
end

if(size(range_H,1) ~= size(H_exp,1))
    error('Not valid lengths');
end

csii= xpar(1);
w0i=xpar(2);
% [Ai,Bi,Ci,Di,Ei,Fi] = xpar(3:end); 
Ai = xpar(3);
Bi = xpar(4);
Ci = xpar(5);
Di = xpar(6);
Ei = xpar(7);
Fi = xpar(8);

% m =1;
% c = 2*xpar(1)*xpar(2);
% k =  xpar(2)^2;

omega = 2*pi.*range_H; 

%vpar=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
%     [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F];

% Controllo su w0i negative
if w0i < 0 
    err0=1000;
else
    err0 = 0;
end


H_anal = (Ai + 1i*Bi)./(-omega.^2 + 2*1i*omega.*(csii*w0i) + w0i^2) + ...
    + (Ci+1i*Di)./ones(size(range_H)) + (Ei+1i*Fi)./(omega.^2+1e-4);

% errore 
c_e = H_anal - H_exp; %complex error
errore = sum(sum(real(c_e).^2) + sum(imag(c_e).^2)) + err0;
end

