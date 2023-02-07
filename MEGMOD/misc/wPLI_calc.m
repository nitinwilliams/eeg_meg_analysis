function [PLI,wPLI]=wPLI_calc(DAT)

[M,~]=size(DAT);

PLI=zeros(M);
wPLI=zeros(M);

for ridx=1:M
    
    for cidx=ridx+1:M

    sig1=DAT(ridx,:);
    sig2=DAT(cidx,:);

    % cross-spectral density
    cdd = sig1.*conj(sig2);
    
    % imaginary part of cross-spectrum
    cdi = imag(cdd);
       
    % PLI (Stam et al. 2007)
    PLI(ridx,cidx) = abs(mean(sign(cdi),2));
    
    % weighted phase-lag index (Vinck et al. 2011)
    wPLI(ridx,cidx) = abs(mean(abs(cdi).*sign(cdi),2))./mean(abs(cdi),2);
    
    end
       
end

end