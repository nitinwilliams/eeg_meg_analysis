clear;

%% Initialisation

numparcs=(148/2);
Nfvec=[3,4,5,7,10,14,20,28,40,57,80,113,135,160,190,226,269,320];
N=length(Nfvec);

load('wmPS_zmod_67subs_gamma0.8to5_5Kreps_lefthem.mat','mPS');
mPS(1:N+1:end)=0;
A{1}=mPS;

load('wmPS_zmod_67subs_gamma0.8to5_5Kreps_righthem.mat','mPS');
mPS(1:N+1:end)=0;
A{2}=mPS;

%% Multi-slice communities at different levels of gamma & omega

gammavec=[1:0.05:1.5];
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

end

matched_partitions=match_mat.*num_modules;
matched_partitions(matched_partitions==0)=NaN;

%% testing different combinations

gamma=1.1;
omega=1;
reps=1000;
BF=zeros(N,N,reps);
CH=zeros(reps,1);

for repidx=1:reps,

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

BF(:,:,repidx)=clustersol_representation(S(:,1));
CH(repidx,1)=(sum(S(:,1)-S(:,2))==0);

end

BF=mean(BF,3);

save('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\submission\CS_CP_gamma1pt1.mat','BF');

figure;
imagesc(BF);
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',1:length(Nfvec),'XTickLabel',Nfvec,'YTick',1:length(Nfvec),'YTickLabel',Nfvec);
set(gcf,'Position',[300,300,600,600]);
