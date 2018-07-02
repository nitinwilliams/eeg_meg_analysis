clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% find mixing range for Louvain

pctshufvec=[0:10:100];
omegavec=[0:0.1:3];
gamma=1;
num_trials=100;
reps=100;

N=74;
load('LFR_binary_networks_shuffledpairs_mixfac0pt5_SCFCmethods.mat','A_shuffledpairs');

trial_nums=randi(size(A_shuffledpairs,1),[num_trials,1]);

optomega_mat1=zeros(num_trials,length(pctshufvec));
optomega_mat2=zeros(num_trials,length(pctshufvec));

for tr_idx=1:num_trials,
    
    EG_shuffledpairs=squeeze(A_shuffledpairs(trial_nums(tr_idx),:,:));

    parfor pctidx=1:length(pctshufvec),
        
      A=cell(1,2);
      match_mat=zeros(1,length(omegavec));
      
        % apply Community Louvain - network 1

        A{1,1}=EG_shuffledpairs{1,pctidx};
                
        M1  = 1:N;                    
        Q0 = -1; Q1 = 0;            
        while Q1-Q0>1e-5;           
           Q0 = Q1;                
           [M1, Q1] = community_louvain(A{1,1},gamma,M1);
        end
        
        % apply Community Louvain - network 2
        
        A{1,2}=EG_shuffledpairs{2,pctidx};
                
         M2  = 1:N;                    
         Q0 = -1; Q1 = 0;            
         while Q1-Q0>1e-5;           
           Q0 = Q1;                
           [M2, Q1] = community_louvain(A{1,2},gamma,M2);
         end
         
         % similarity between partitions
         
         partsim_val=partsim(M1,M2);
         zQvec=cell(1,1);
         
         for omidx=1:length(omegavec),     
            omega=omegavec(omidx);
            [S,~,zQvec{1,1}(1,omidx)]=multislice_communitydetection(A,N,gamma,omega,reps,'nodalnullmodel');
            match_mat(1,omidx)=(sum(S(:,1)-S(:,2))==0);
         end

          cind=find(match_mat,1,'first');
          optomega_mat1(tr_idx,pctidx)=omegavec(cind)*partsim_val;

          [~,cind]=max(zQvec{1,1});
          optomega_mat2(tr_idx,pctidx)=omegavec(cind);
          
    end
    
    display(tr_idx);

end

save('B:\Nitin\SCFC_methods\data\LFR_data\optomegavals_binary_mixfac0pt5.mat','optomega_mat1','optomega_mat2');

% plotting

figure;
imagesc(optomega_mat);
colorbar;
set(gca,'YDir','Normal');
axis square;

figure;
plot(mean(optomega_mat),'k-*');
ylim([0,2]);
axis square;