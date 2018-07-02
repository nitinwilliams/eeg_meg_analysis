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

Nf={'9.83'};

numparcs=74;
meg_lefthem_idxs=[1:1*numparcs];
meg_righthem_idxs=[1*numparcs+1:2*numparcs];
sc_lefthem_idxs=[2*numparcs+1:3*numparcs];
sc_righthem_idxs=[3*numparcs+1:4*numparcs];

pctval=0;
N=numparcs;
perms=100;
gamma=2;
nct=0;

load(['B:\Nitin\SCFC_methods\data\sc_meg_data\Pstabmat_weighted_Louvain_SC_FC_MEG_' Nf{1,1} 'Hz_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval) '.mat'],'Pstabmat','Pstabmat_rand','Sb_surr');  

%% Left hemisphere

Sb_surr_MEG=Sb_surr(meg_lefthem_idxs,:);
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));

MASK_MEG=Pstabmat(meg_lefthem_idxs,:)>Pstabmat_thr(meg_lefthem_idxs,:);

Sb_surr_SC=Sb_surr(sc_lefthem_idxs,:);
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK_SC=Pstabmat(sc_lefthem_idxs,:)>Pstabmat_thr(sc_lefthem_idxs,:);

MASK=(MASK_MEG.*MASK_SC);

numreps=size(Sb_surr_MEG,2);
PS_lefthem=zeros(1,numreps);
PS_nulldist_lefthem=zeros(perms,numreps);

for ridx=1:numreps,
     
   Sb_MEG=Sb_surr_MEG(:,ridx);
   Sb_SC=Sb_surr_SC(:,ridx);
    
   val_idxs=find(MASK(:,ridx));
   PS_lefthem(1,ridx)=partsim(Sb_MEG(val_idxs),Sb_SC(val_idxs));
   PS_nulldist_lefthem(:,ridx)=partsim_nulldist(Sb_MEG(val_idxs),Sb_SC(val_idxs),perms);
   
   display(ridx); 
   
end

zPS_lefthem=(PS_lefthem-mean(PS_nulldist_lefthem))./std(PS_nulldist_lefthem);

%% Right hemisphere

Sb_surr_MEG=Sb_surr(meg_righthem_idxs,:);
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));

MASK_MEG=Pstabmat(meg_righthem_idxs,:)>Pstabmat_thr(meg_righthem_idxs,:);

Sb_surr_SC=Sb_surr(sc_righthem_idxs,:);
Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK_SC=Pstabmat(sc_righthem_idxs,:)>Pstabmat_thr(sc_righthem_idxs,:);

MASK=(MASK_MEG.*MASK_SC);

numreps=size(Sb_surr_MEG,2);
PS_righthem=zeros(1,numreps);
PS_nulldist_righthem=zeros(perms,numreps);

for ridx=1:numreps,
     
   Sb_MEG=Sb_surr_MEG(:,ridx);
   Sb_SC=Sb_surr_SC(:,ridx);
    
   val_idxs=find(MASK(:,ridx));
   PS_righthem(1,ridx)=partsim(Sb_MEG(val_idxs),Sb_SC(val_idxs));
   PS_nulldist_righthem(:,ridx)=partsim_nulldist(Sb_MEG(val_idxs),Sb_SC(val_idxs),perms);
   
   display(ridx); 
   
end

zPS_righthem=(PS_righthem-mean(PS_nulldist_righthem))./std(PS_nulldist_righthem);

PS_bothhems=mean([PS_lefthem;PS_righthem]);
zPS_bothhems=mean([zPS_lefthem;zPS_righthem]);

save(['B:\Nitin\SCFC_methods\data\sc_meg_data\partsim_weighted_GenLouvain_SC_MEG_FC_' Nf{1,1} 'Hz_parc2k9_67subs_gamma' num2str(gamma) '.mat'],'PS_lefthem','PS_righthem','PS_bothhems');
save(['B:\Nitin\SCFC_methods\data\sc_meg_data\zpartsim_weighted_GenLouvain_SC_MEG_FC_' Nf{1,1} 'Hz_parc2k9_67subs_gamma' num2str(gamma) '.mat'],'zPS_lefthem','zPS_righthem','zPS_bothhems');
