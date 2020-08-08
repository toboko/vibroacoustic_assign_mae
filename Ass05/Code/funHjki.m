function H_reco = funHjki(vpar, range_H) 
% RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
% if (size(vpar,3) ~= 9)
%     error('Not valid vpar');
%end

%[m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F] = [1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)]; %non credo sia giusto
m_i = vpar(:,:,1); 
c_i = vpar(:,:,2); 
k_i = vpar(:,:,3); 
Aj = vpar(:,:,4);
Bj = vpar(:,:,5);
Cj = vpar(:,:,6);
Dj = vpar(:,:,7);
Ej = vpar(:,:,8);
Fj = vpar(:,:,9);

w0i = sqrt(k_i./m_i);
csii = c_i ./ (2.*m_i.*w0i); 

omega = 2*pi.*range_H;

if(nnz(ismember(range_H,0)) ~=0)
    omega = omega + 1e-04;
end

%vpar=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
%     [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F];

%H_reco = (Aj)./(- m_i.*omega.^2 + 1i*c_i.*omega + k_i);
% H_reco = (Aj./m_i)./(-omega.^2 + 2*1i*omega.*(csii.*w0i) + w0i.^2);
H_reco = (Aj + 1i*Bj)./(-omega.^2 + 2*1i*omega.*(csii.*w0i) + w0i.^2) ...
          + (Cj + 1i*Dj).*ones(size(omega)) + (Ej + 1i*Fj)./(omega.^2);
end
