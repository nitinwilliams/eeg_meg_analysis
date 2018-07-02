clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% find mixing range for Louvain

N=74;
load('LFR_binary_networks_mixfac0pt9_SCFCmethods.mat','A','CP');
pctshufvec=[0:10:100];

A_shuffledpairs=cell(size(A,1),2,length(pctshufvec));
CP_shuffledpairs=cell(size(A,1),2,length(pctshufvec));

for pctidx=1:length(pctshufvec),
    
    pctval=pctshufvec(pctidx);
    
    for netidx=1:size(A,1),
          
    D=A{netidx};      
    V=CP{netidx};

    [D_pair,V_pair]=network_shuffle(D,V,pctval);
    
    A_shuffledpairs{netidx,1,pctidx}=(D_pair{1}+D_pair{1}')/2;
    A_shuffledpairs{netidx,2,pctidx}=(D_pair{2}+D_pair{2}')/2;

    CP_shuffledpairs{netidx,1,pctidx}=V_pair{1};
    CP_shuffledpairs{netidx,2,pctidx}=V_pair{2};
    
    end
    
    display(pctidx);
    
end

save('LFR_binary_networks_shuffledpairs_mixfac0pt9_SCFCmethods.mat','A_shuffledpairs','CP_shuffledpairs');

