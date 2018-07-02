function CSR=clustersol_representation(S)

Z=max(S);
M=length(S);
CSR=zeros(M,M,Z);

for zidx=1:Z,
    A=double(S==zidx);
    CSR(:,:,zidx)=A*A';        
end

CSR=sum(CSR,3);

end