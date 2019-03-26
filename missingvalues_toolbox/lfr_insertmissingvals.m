clear;
% addpath(genpath('/home/nw04/fcmodules_codedata/Nitin/'));

%% insert missing vals

load('LFRbenchmark_networks','A');

N=1000;
reps=10;
prop_mv_vec=[0.1:0.1:0.9];
A_mv=cell(length(prop_mv_vec),size(A,2));

for mix_idx=1:size(A,2),

    for prop_idx=1:length(prop_mv_vec),
    
        D=zeros(N,N,reps);
        B=A{1,mix_idx};
        
        for repidx=1:reps,
    
          D(:,:,repidx)=insert_missingvals(B,prop_mv_vec(prop_idx),[],'random');  
    
        end
    
       A_mv{prop_idx,mix_idx}=D;
        
    end
    
    display(mix_idx);

end

save('B:\Nitin\MV_paper\data\LFRbenchmark_networks_wmvals.mat','A_mv','-v7.3');