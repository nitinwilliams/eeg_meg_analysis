clear;
addpath(genpath('/home/nw04/fcmodules_codedata/Nitin/'));

%% Partition similarity for networks with missing values

N=1000;
load('LFRbenchmark_networks.mat','CP');
load('LFRbenchmark_networks_wmvals.mat','A_mv');

reps=1000;
num_mix_coeffs=8;
num_prop_coeffs=9;

PS_Infomap_zeros=zeros(num_prop_coeffs,num_mix_coeffs,reps);
PS_Infomap_rcmean=zeros(num_prop_coeffs,num_mix_coeffs,reps);
PS_Infomap_cneigh=zeros(num_prop_coeffs,num_mix_coeffs,reps);
PS_Infomap_conclus=zeros(num_prop_coeffs,num_mix_coeffs,reps);

Nc_Infomap_zeros=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Nc_Infomap_rcmean=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Nc_Infomap_cneigh=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Nc_Infomap_conclus=zeros(num_prop_coeffs,num_mix_coeffs,reps);

Tc_Infomap_zeros=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Tc_Infomap_rcmean=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Tc_Infomap_cneigh=zeros(num_prop_coeffs,num_mix_coeffs,reps);
Tc_Infomap_conclus=zeros(num_prop_coeffs,num_mix_coeffs,reps);

for prop_idx=1:num_prop_coeffs,
       
    for mix_idx=1:num_mix_coeffs,
        
        V_ref=CP{mix_idx};
        DAT=A_mv{prop_idx,mix_idx};
        
        for repidx=1:reps,
            
           D=DAT(:,:,repidx);
                       
           % Fill zeros
           
           tic;
           Df=fill_zeros(D);
           V=infomap(Df);
           Tc_Infomap_zeros(prop_idx,mix_idx,repidx)=toc;
           PS_Infomap_zeros(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
           Nc_Infomap_zeros(prop_idx,mix_idx,repidx)=max(V);
           
           % Row-column mean
           
           tic;
           Df=fill_rcmean(D,0.5,'binary');
           V=infomap(Df);
           Tc_Infomap_rcmean(prop_idx,mix_idx,repidx)=toc;
           PS_Infomap_rcmean(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
           Nc_Infomap_rcmean(prop_idx,mix_idx,repidx)=max(V);
           
           % Common neighbours
          
           tic;
           Df=fill_commonneighbours(D,0.5,'binary');
           V=infomap(Df);
           Tc_Infomap_cneigh(prop_idx,mix_idx,repidx)=toc;
           PS_Infomap_cneigh(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
           Nc_Infomap_cneigh(prop_idx,mix_idx,repidx)=max(V);
          
           % Consensus clustering
          
           tic;
           Df=fill_consensusclustering(D,1,100,0.5,'Infomap','binary');
           tval=prctile(Df(:),90);
           Df(Df<tval)=0;
           Df(Df>=tval)=1;
           V=infomap(Df);
           Tc_Infomap_conclus(prop_idx,mix_idx,repidx)=toc;
           PS_Infomap_conclus(prop_idx,mix_idx,repidx)=partsim(V,V_ref);
           Nc_Infomap_conclus(prop_idx,mix_idx,repidx)=max(V);
                                         
        end
    
        display(mix_idx);
        
        display(prop_idx);
        
    end    
        
end

save('/home/nw04/fcmodules_codedata/Nitin/MV_paper/data/LFRnetworks_binary_Infomap_PS.mat','PS_Infomap_zeros','PS_Infomap_rcmean','PS_Infomap_cneigh','PS_Infomap_conclus',...
    'Nc_Infomap_zeros','Nc_Infomap_rcmean','Nc_Infomap_cneigh','Nc_Infomap_conclus',...
    'Tc_Infomap_zeros','Tc_Infomap_rcmean','Tc_Infomap_cneigh','Tc_Infomap_conclus');