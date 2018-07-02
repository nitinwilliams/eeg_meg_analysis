function [A,CP]=LFR_networks(N,mode)

load('network.dat');

D=zeros(N);
edge_indices=sub2ind([N,N],network(:,1),network(:,2));

mu_wm=0.7;
sigma_wm=0.1;

mu_bm=0.3;
sigma_bm=0.1;

D(edge_indices)=1;

load('community.dat');

CP = community(:,2);

switch mode
    
    case 'binary'
    A = D;
    
    case 'weighted' 
    W=clustersol_representation(CP);
    wm_idxs=find(W(:)==1);
    bm_idxs=find(W(:)==0);
    W(wm_idxs)=randn_positive(length(wm_idxs),mu_wm,sigma_wm);
    W(bm_idxs)=randn_positive(length(bm_idxs),mu_bm,sigma_bm);    
    A = W.*D;
    
end

end