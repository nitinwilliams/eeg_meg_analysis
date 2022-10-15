clear;
addpath(genpath('B:\Nitin\HBPconnectome\'));

%% Initialisation

numparcs=(148/2);
Nfvec=[3,4,5,7,10,14,20,28,40,57,80,113,135,160,190,226,269,320];
N=numparcs;

load('B:\Nitin\HBPconnectome\GA_data\submission\mSb_Louvain_pct80_CB_lefthem_wFC_67subs_parc2k9_acrossgvals_acrossfreqs.mat','Sb');
Sb(1:N+1:end)=0;
A{2}=Sb;

load('B:\Nitin\HBPconnectome\GA_data\submission\mSb_Louvain_pct80_CB_righthem_wFC_67subs_parc2k9_acrossgvals_acrossfreqs.mat','Sb');
Sb(1:N+1:end)=0;
A{1}=Sb;

%% Multi-slice communities at different levels of gamma & omega

gammavec=[1.1:0.1:2];
omegavec=[0.1:0.1:1];

match_mat=zeros(length(gammavec),length(omegavec));
num_modules=zeros(length(gammavec),length(omegavec));

for gamidx=1:length(gammavec),
    
    for omidx=1:length(omegavec),
        
        gamma=gammavec(gamidx);
        omega=omegavec(omidx);

        T=length(A);
        B=spalloc(N*T,N*T,N*N*T+2*N*T);
        twomu=0;
        for s=1:T
            k=sum(A{s});
            twom=sum(k);
            twomu=twomu+twom;
            indx=[1:N]+(s-1)*N;
            B(indx,indx)=A{s}-gamma*k'*k/twom;
        end
    
        twomu=twomu+2*omega*N*(T-1);
        B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        [S,Q] = genlouvain(B);
        Q = Q/twomu;
        S = reshape(S,N,T);
        
        match_mat(gamidx,omidx)=(sum(S(:,1)-S(:,2))==0);
        num_modules(gamidx,omidx)=max(S(:));

    end

    Sf=S(:,1); save(['B:\Nitin\HBPconnectome\GA_data\submission\Sfss_Louvain_pct80_gamma' num2str(gammavec(gamidx)) '_CB_lefthem_wFC_67subs_parc2k9_acrossgvals_acrossfreqs.mat'],'Sf');
    Sf=S(:,2); save(['B:\Nitin\HBPconnectome\GA_data\submission\Sfss_Louvain_pct80_gamma' num2str(gammavec(gamidx)) '_CB_righthem_wFC_67subs_parc2k9_acrossgvals_acrossfreqs.mat'],'Sf');
    
end

matched_partitions=match_mat.*num_modules;
matched_partitions(matched_partitions==0)=NaN;
