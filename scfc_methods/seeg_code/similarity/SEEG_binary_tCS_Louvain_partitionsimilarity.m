clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% Parcel names

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

%% Initialisation

Nf={'9.99'};

numparcs=74;
lefthem_idxs=[1:2:numparcs*2];
righthem_idxs=[2:2:numparcs*2];

pctval1=80;
pctval2=80;
N=numparcs;
perms=100;
gamma=2;
nct=0;

%% Left hemisphere

load(['B:\Nitin\SCFC_methods\data\seeg_data\Sb_binary_surr_gamma' num2str(gamma) '_pctval' int2str(pctval1) '_thr' int2str(nct) '_' Nf{1,1} 'Hz_lefthem_wTFC_parc2k9.mat'],'Sb');
Sb_surr_SEEG=Sb;

load(['B:\Nitin\SCFC_methods\data\seeg_data\Pstabmat_binary_Louvain_SEEG_FC_' Nf{1,1} 'Hz_lefthem_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval1) '.mat'],'Pstabmat','Pstabmat_rand');  
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK_SEEG=Pstabmat>Pstabmat_thr;

load(['B:\Nitin\SCFC_methods\data\sc_data\Sb_binary_SC_surr_gamma' num2str(gamma) '_pctval' int2str(pctval2) '_thr' int2str(nct) '_lefthem_wTFC_parc2k9.mat'],'Sb');
Sb_surr_SC=Sb;

load(['B:\Nitin\SCFC_methods\data\sc_data\Pstabmat_binary_Louvain_SC_lefthem_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval2) '.mat'],'Pstabmat','Pstabmat_rand');  
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK_SC=Pstabmat>Pstabmat_thr;

MASK=(MASK_SEEG.*MASK_SC);

numreps=size(Sb_surr_SEEG,1);
PS_lefthem=zeros(1,numreps);
PS_nulldist_lefthem=zeros(perms,numreps);

for ridx=1:numreps,
     
   Sb_SEEG=Sb_surr_SEEG{ridx,1};
   Sb_SC=Sb_surr_SC{ridx,1};
    
   val_idxs=find(MASK(:,ridx));
   PS_lefthem(1,ridx)=partsim(Sb_SEEG(val_idxs),Sb_SC(val_idxs));
   PS_nulldist_lefthem(:,ridx)=partsim_nulldist(Sb_SEEG(val_idxs),Sb_SC(val_idxs),perms);
   
   display(ridx); 
   
end

zPS_lefthem=(PS_lefthem-mean(PS_nulldist_lefthem))./std(PS_nulldist_lefthem);

%% Right hemisphere

load(['B:\Nitin\SCFC_methods\data\seeg_data\Sb_binary_surr_gamma' num2str(gamma) '_pctval' int2str(pctval1) '_thr' int2str(nct) '_' Nf{1,1} 'Hz_righthem_wTFC_parc2k9.mat'],'Sb');
Sb_surr_SEEG=Sb;

load(['B:\Nitin\SCFC_methods\data\seeg_data\Pstabmat_binary_Louvain_SEEG_FC_' Nf{1,1} 'Hz_righthem_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval1) '.mat'],'Pstabmat','Pstabmat_rand');  
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK_SEEG=Pstabmat>Pstabmat_thr;

load(['B:\Nitin\SCFC_methods\data\sc_data\Sb_binary_SC_surr_gamma' num2str(gamma) '_pctval' int2str(pctval2) '_thr' int2str(nct) '_righthem_wTFC_parc2k9.mat'],'Sb');
Sb_surr_SC=Sb;

load(['B:\Nitin\SCFC_methods\data\sc_data\Pstabmat_binary_Louvain_SC_righthem_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval2) '.mat'],'Pstabmat','Pstabmat_rand');  
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK_SC=Pstabmat>Pstabmat_thr;

MASK=(MASK_SEEG.*MASK_SC);

numreps=size(Sb_surr_SEEG,1);
PS_righthem=zeros(1,numreps);
PS_nulldist_righthem=zeros(perms,numreps);

for ridx=1:numreps,
     
   Sb_SEEG=Sb_surr_SEEG{ridx,1};
   Sb_SC=Sb_surr_SC{ridx,1};
    
   val_idxs=find(MASK(:,ridx));
   PS_righthem(1,ridx)=partsim(Sb_SEEG(val_idxs),Sb_SC(val_idxs));
   PS_nulldist_righthem(:,ridx)=partsim_nulldist(Sb_SEEG(val_idxs),Sb_SC(val_idxs),perms);
   
   display(ridx); 
   
end

zPS_righthem=(PS_righthem-mean(PS_nulldist_righthem))./std(PS_nulldist_righthem);

PS_bothhems=mean([PS_lefthem;PS_righthem]);
zPS_bothhems=mean([zPS_lefthem;zPS_righthem]);

save(['B:\Nitin\SCFC_methods\data\seeg_data\partsim_binary_Louvain_SC_SEEG_FC_' Nf{1,1} 'Hz_parc2k9_67subs_gamma' num2str(gamma) '.mat'],'PS_lefthem','PS_righthem','PS_bothhems');
save(['B:\Nitin\SCFC_methods\data\seeg_data\zpartsim_binary_Louvain_SC_SEEG_FC_' Nf{1,1} 'Hz_parc2k9_67subs_gamma' num2str(gamma) '.mat'],'zPS_lefthem','zPS_righthem','zPS_bothhems');
