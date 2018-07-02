clear;
% addpath(genpath('B:\Nitin\SCFC_methods\'));

%% Initialisation

N=74;
gamma=2;
lefthem_idxs=[1:2:N*2];
righthem_idxs=[2:2:N*2];

% Parcel names

parc_acronyms_parc2K9={'mrgF_L','mrgF_R','iO_L','iO_R','paC_L','paC_R','subC_L','subC_R','trF_L','trF_R','Cla_L','Cla_R',...
    'aClm_L','aClm_R','pClm_L','pClm_R','dClp_L','dClp_R','vClp_L','vClp_R','CN_L','CN_R','iFGopc_L','iFGopc_R','iFGorb_L','iFGorbR',...
    'iFGtri_L','iFGtri_R','mFG_L','mFG_R','sFG_L','sFG_R','lgING_L','lgING_R','shING_L','shING_R','mOG_L','mOG_R','sOG_L','sOG_R',...
    'OTGfus_L','OTGfus_R','OTGling_L','OTGling_R','OTGhip_L','OTGhip_R','orbG_L','orbG_R','iPGang_L','iPGang_R','iPGsup_L','iPGsup_R',...
    'sPG_L','sPG_R','poCG_L','poCG_R','prCG_L','prCG_R','prCN_L','prCN_R','rcG_L','rcG_R','subcalG_L','subcalG_R','sTG-trTG_L','sTG-trTG_R',...
    'sTGla_L','sTGla_R','sTGpp_L','sTGpp_R','sTGpt_L','sTGpt_R','iTG_L','iTG_R','mTG_L','mTG_R','laShz_L','laShz_R','laSvr_L','laSvr_R',...
    'laSp_L','laSp_R','Opole_L','Opole_R','Tpole_L','Tpole_R','caS_L','caS_R','CS_L','CS_R','ClSmrg_L','ClSmrg_R','aINS_L','aINS_R',...
    'iINS_L','iINS_R','sINS_L','sINS_R','acolS_L','acolS_R','pcolS_L','pcolS_R','iFS_L','iFS_R','mFS_L','mFS_R','sFS_L','sFS_R','inS_L','inS_R',...
    'intPS_L','intPS_R','mOS_L','mOS_R','sOS_L','sOS_R','aOS_L','aOS_R','laOTS_L','laOTS_R','meOTS_L','meOTS_R','orbSla_L','orbSla_R',...
    'orbSme_L','orbSme_R','orbS_L','orbS_R','POS_L','POS_R','pecalS_L','pecalS_R','poCS_L','poCS_R','iprCS_L','iprCS_R','sprCS_L','sprCS_R',...
    'suborbS_L','suborbS_R','subPS_L','subPS_R','iTS_L','iTS_R','sTS_L','sTS_R','trTS_L','trTS_R'};

%% measuring z-scored modularity

load('B:\Nitin\SCFC_methods\data\sc_data\Sb_weighted_SC_gamma2_pctval0_thr0_lefthem_wTFC_parc2k9.mat','Zmat');

D=Zmat{1,1};

% original modularity

M  = 1:N;
Q0 = -1; Q = 0;
while Q-Q0>1e-5;
    Q0 = Q;
    [M, Q] = community_louvain(D,gamma,M);
end

Q_orig=Q;

% [~,reorder]=sort(M);
% 
% figure;
% imagesc(D(reorder,reorder));
% colorbar;
% axis square;
% 
% set(gca,'YDir','Normal');

% modularity of equivalent random networks

reps=100;
Q_null=zeros(reps,1);

for repidx=1:reps,
    
   D_rand=cell2mat(randomise_und_ensemble(mat2cell(D,N,N)));
   
   M  = 1:N;
   Q0 = -1; Q = 0;
   while Q-Q0>1e-5;
      Q0 = Q;
      [M, Q] = community_louvain(D_rand,gamma,M);
   end
   
%    [~,reorder]=sort(M);
% 
%     figure;
%     imagesc(D_rand(reorder,reorder));
%     colorbar;
%     axis square;
% 
%     set(gca,'YDir','Normal');
   
   Q_null(repidx,1)=Q;
   
end

zQ=(Q_orig-mean(Q_null))/std(Q_null);

display(zQ);