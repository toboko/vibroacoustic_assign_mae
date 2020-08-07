function FRF_reco = reco(vparre, reinge)
%REDO Bella
%   RICORDA: vpar(a,b,c) = a:picc b: mis c:parametri
    nPeaks = size(vparre,1);
    meas = size(vparre,2);
    m_i = vparre(:,:,1); 
    c_i = vparre(:,:,2); 
    k_i = vparre(:,:,3); 
    Aj = vparre(:,:,4);
    Bj = vparre(:,:,5);
    omega = 2*pi.*reinge;
    FRF_reco = zeros(length(reinge), meas);
    
    for mm = 1:meas
        for pp = 1:nPeaks
            FRF_reco(:, mm) = FRF_reco(: , mm) + ...
                (Aj(pp,mm) + 1i*Bj(pp,mm))./(- m_i(pp,mm).*omega.^2 + 1i*c_i(pp,mm).*omega + k_i(pp,mm));
        end
    end
    

end

