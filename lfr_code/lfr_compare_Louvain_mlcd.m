clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% partition similarity w.r.t. ground truth

pctshufvec=[0:10:100];
gamma=1;
reps=0;

N=74;
load('LFR_binary_networks_shuffledpairs_mixfac0pt5_SCFCmethods.mat','A_shuffledpairs','CP_shuffledpairs');
load('optomegavals_binary_mixfac0pt5.mat','optomega_mat1','optomega_mat2'); 

optomega_vec1=mean(optomega_mat1);
optomega_vec2=mean(optomega_mat2);

PMAT=cell(1,4);

for tr_idx=1:size(A_shuffledpairs,1),

EG_A_shuffledpairs=squeeze(A_shuffledpairs(tr_idx,:,:));
EG_CP_shuffledpairs=squeeze(CP_shuffledpairs(tr_idx,:,:));

    for pctidx=1:length(pctshufvec),
         
        A=cell(1,2);
        CP=cell(1,2);
        
        A{1,1}=EG_A_shuffledpairs{1,pctidx};
        A{1,2}=EG_A_shuffledpairs{2,pctidx};
        
        CP{1,1}=EG_CP_shuffledpairs{1,pctidx};
        CP{1,2}=EG_CP_shuffledpairs{2,pctidx};
        
       % find partition similarity w.r.t. ground truth  (Louvain)      
        
        % apply Community Louvain - network 1
               
        M1  = 1:N;                    
        Q0 = -1; Q1 = 0;            
        while Q1-Q0>1e-5;           
           Q0 = Q1;                
           [M1, Q1] = community_louvain(A{1,1},gamma,M1);
        end
        
        % apply Community Louvain - network 2
                       
         M2  = 1:N;                    
         Q0 = -1; Q1 = 0;            
         while Q1-Q0>1e-5;           
           Q0 = Q1;                
           [M2, Q1] = community_louvain(A{1,2},gamma,M2);
         end

         PMAT{1,1}(tr_idx,pctidx)=mean([partsim(M1,CP{1}),partsim(M2,CP{2})]);
         
        % find partition similarity w.r.t. ground truth (omega=1)
         
        omega=1;
        [S,~]=multislice_communitydetection(A,N,gamma,omega,reps,'nonullmodel');
        
        PMAT{1,2}(tr_idx,pctidx)=mean([partsim(S(:,1),CP{1}),partsim(S(:,2),CP{2})]);
         
        % find partition similarity w.r.t. ground truth (MLCD - 1)
         
        omega=optomega_vec1(pctidx);
        [S,~]=multislice_communitydetection(A,N,gamma,omega,reps,'nonullmodel');
        
        PMAT{1,3}(tr_idx,pctidx)=mean([partsim(S(:,1),CP{1}),partsim(S(:,2),CP{2})]);
        
        % find partition similarity w.r.t. ground truth (MLCD - 2)
         
        omega=optomega_vec2(pctidx);
        [S,~]=multislice_communitydetection(A,N,gamma,omega,reps,'nonullmodel');
        
        PMAT{1,4}(tr_idx,pctidx)=mean([partsim(S(:,1),CP{1}),partsim(S(:,2),CP{2})]);
        
    end

    display(tr_idx);
    
end

save('B:\Nitin\SCFC_methods\data\LFR_data\partsim_binary_comparison_Louvain_mlcd_mixfac0pt5.mat','PMAT');

%% plotting

D1=PMAT{1,1};
D2=PMAT{1,2};
D3=PMAT{1,3};
D4=PMAT{1,4};

figure;
hold on;
plot(pctshufvec,mean(D1),'k*');
plot(pctshufvec,mean(D2),'m*');
plot(pctshufvec,mean(D3),'b*');
plot(pctshufvec,mean(D4),'r*');
ylim([0,1]);
xlabel('Percentage of shuffled nodes');
ylabel('Partition similarity');