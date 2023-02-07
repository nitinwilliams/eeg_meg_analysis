function [PLVvals,PLVmdir,iPLVvals,iPLVnormvals]=PLVcalc(DAT,mode)

% PLV calculation 
[M,N]=size(DAT);

% Obtaining instantaneous phases
switch mode
    case 'hilbert'
    DAT=unwrap(angle(hilbert(DAT')))';
    case 'morlet'
    DAT=unwrap(angle(DAT'))';
end

PLVvals=zeros(M);
PLVmdir=zeros(M);
iPLVvals=zeros(M);
iPLVnormvals=zeros(M);

for ridx=1:M

    for cidx=ridx:M
           
        % instantaneous phase difference    
        pD=DAT(ridx,:)-DAT(cidx,:);
        
        % PLV values
        PLVvals(ridx,cidx)=abs(sum(exp(1i*pD))/N);
    
        % PLV mean direction
        PLVmdir(ridx,cidx)=angle(sum(exp(1i*pD))/N);
        
        % imaginary PLV values
        iPLVvals(ridx,cidx)=abs(imag(sum(exp(1i*pD))/N));
        
        % normalised imaginary PLV values
        iPLVnormvals(ridx,cidx)=abs(imag(sum(exp(1i*pD))/N))/abs(sum(exp(1i*pD))/N);
        
    end
        
end

end