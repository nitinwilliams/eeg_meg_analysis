function [Qval,NC]=partitionquality_ss(A,gamma)

numc=size(A,1);
N=length(A{1,1});

BF=zeros(N,N,numc);

for ncidx=1:numc

    Z=A{ncidx,1};   
    
    M  = 1:N;                    
    Q0 = -1; Q1 = 0;            
    while Q1-Q0>1e-5          
    Q0 = Q1;                
    [M, Q1]=community_louvain(Z,gamma,M);
    end
            
    BF(:,:,ncidx)=clustersol_representation(M);
            
end

BF=mean(BF,3);   
BF(1:N+1:end)=0;

M  = 1:N;                    
Q0 = -1; Q1 = 0;            
while Q1-Q0>1e-5         
    Q0 = Q1;                
    [M, Q1] = community_louvain(BF,gamma,M);
end

Qval=Q1;
NC=length(unique(M));

end