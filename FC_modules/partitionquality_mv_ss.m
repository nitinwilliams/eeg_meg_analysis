function [Qval,NC]=partitionquality_mv_ss(A,gamma)

numc=size(A,1);
N=length(A{1,1});
T=size(A,2);

Sbmat=zeros(N,T);      

for tidx=1:T,

BF=zeros(N,N,numc);

for ncidx=1:numc,

    Z=A{ncidx,tidx};   
    
    M  = 1:N;                    
    Q0 = -1; Q1 = 0;            
    while Q1-Q0>1e-5;           
    Q0 = Q1;                
    [M, Q1]=community_louvain(Z,gamma,M);
    end
            
    BF(:,:,ncidx)=clustersol_representation(M);
            
end

BF=mean(BF,3);   
BF(1:N+1:end)=0;

M  = 1:N;                    
Q0 = -1; Q1 = 0;            
while Q1-Q0>1e-5;           
    Q0 = Q1;                
    [M, Q1] = community_louvain(BF,gamma,M);
end

Sbmat(:,tidx)=M;

end

% final (consensus) partition

SSM=zeros(N,N,T);
for pidx=1:T,
    SSM(:,:,pidx)=clustersol_representation(Sbmat(:,pidx));
end
        
SSM=mean(SSM,3);
SSM(1:N+1:end)=0;
    
M  = 1:N;                    
Q0 = -1; Q1 = 0;            
while Q1-Q0>1e-5;           
     Q0 = Q1;                
     [M, Q1] = community_louvain(SSM,gamma,M);
end

Qval=Q1;
NC=length(unique(M));

end