clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% find mixing range for Louvain

N=74;
load('LFR_binary_networks_mixfac0pt9_SCFCmethods.mat','A','CP');

gamma=1;
PS_Louvain=zeros(length(A),1);

for idx=1:size(A,1),
    
D=A{idx};      
    
% Louvain community detection

V1  = [1:N]';
Q0 = -1; Q1 = 0;            
while Q1-Q0>1e-5;           
   Q0 = Q1;                
   [V1, Q1] = community_louvain(D,gamma,V1);
end  

V_ref=CP{idx};

PS_Louvain(idx,1)=partsim(V1,V_ref); 

display(idx);
    
end

figure;
histogram(PS_Louvain,[0:0.1:1]);

