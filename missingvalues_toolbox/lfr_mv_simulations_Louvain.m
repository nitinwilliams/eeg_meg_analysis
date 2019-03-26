clear;
addpath(genpath('B:\Nitin\HBPconnectome\'));

%% Partition similarity for networks with missing values

N=1000;
load('LFRbenchmark_networks.mat','CP');
load('LFRbenchmark_networks_wmvals.mat','A_mv');

reps=1000;
num_mix_coeffs=8;
num_prop_coeffs=9;

PS_Louvain_zeros=zeros(num_prop_coeffs,num_mix_coeffs,reps);
PS_Louvain_rcmean=zeros(num_prop_coeffs,num_mix_coeffs,reps);
PS_Louvain_cneigh=zeros(num_prop_coeffs,num_mix_coeffs,reps);
PS_Louvain_conclus=zeros(num_prop_coeffs,num_mix_coeffs,reps);

Nc_Louvain_zeros=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Nc_Louvain_rcmean=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Nc_Louvain_cneigh=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Nc_Louvain_conclus=zeros(num_prop_coeffs,num_mix_coeffs,reps);

Tc_Louvain_zeros=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Tc_Louvain_rcmean=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Tc_Louvain_cneigh=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Tc_Louvain_conclus=zeros(num_prop_coeffs,num_mix_coeffs,reps);

for prop_idx=1:num_prop_coeffs,
    
    parfor mix_idx=1:num_mix_coeffs,
        
        V_ref=CP{mix_idx};
        DAT=A_mv{prop_idx,mix_idx};
        
        for repidx=1:reps,
            
           D=DAT(:,:,repidx);
                
           % Fill zeros
           
           tic;
           Df=fill_zeros(D);
           
           V  = [1:N]';
           Q0 = -1; Q1 = 0;            
           while Q1-Q0>1e-5;           
              Q0 = Q1;                
              [V, Q1] = community_louvain(Df,1,V);
           end  
           
           Tc_Louvain_zeros(prop_idx,mix_idx,repidx)=toc;
           PS_Louvain_zeros(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
           Nc_Louvain_zeros(prop_idx,mix_idx,repidx)=max(V);
           
           % Row-column mean
           
           tic;
           Df=fill_rcmean(D,0.5,'binary');
           
           V  = [1:N]';
           Q0 = -1; Q1 = 0;            
           while Q1-Q0>1e-5;           
              Q0 = Q1;                
              [V, Q1] = community_louvain(Df,1,V);
           end  
           
          Tc_Louvain_rcmean(prop_idx,mix_idx,repidx)=toc;
          PS_Louvain_rcmean(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
          Nc_Louvain_rcmean(prop_idx,mix_idx,repidx)=max(V);
           
          % Common neighbours
          
          tic;
          Df=fill_commonneighbours(D,0.5,'binary');
           
          V  = [1:N]';
          Q0 = -1; Q1 = 0;            
          while Q1-Q0>1e-5;           
            Q0 = Q1;                
            [V, Q1] = community_louvain(Df,1,V);
          end  
           
          Tc_Louvain_cneigh(prop_idx,mix_idx,repidx)=toc;
          PS_Louvain_cneigh(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
          Nc_Louvain_cneigh(prop_idx,mix_idx,repidx)=max(V);
          
          % Consensus clustering
          
          tic;
          Df=fill_consensusclustering(D,1,100,0.5,'Louvain','binary');
           
          V  = [1:N]';
          Q0 = -1; Q1 = 0;            
          while Q1-Q0>1e-5;           
             Q0 = Q1;                
             [V, Q1] = community_louvain(Df,1,V);
          end  
           
          Tc_Louvain_conclus(prop_idx,mix_idx,repidx)=toc;
          PS_Louvain_conclus(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
          Nc_Louvain_conclus(prop_idx,mix_idx,repidx)=max(V);
                     
        end
        
        display(mix_idx);
        
        display(prop_idx);
        
    end    
        
end

save('B:\Nitin\MV_paper\data\LFRnetworks_binary_Louvain_PS.mat','PS_Louvain_zeros','PS_Louvain_rcmean','PS_Louvain_cneigh','PS_Louvain_conclus',...
    'Nc_Louvain_zeros','Nc_Louvain_rcmean','Nc_Louvain_cneigh','Nc_Louvain_conclus',...
    'Tc_Louvain_zeros','Tc_Louvain_rcmean','Tc_Louvain_cneigh','Tc_Louvain_conclus');