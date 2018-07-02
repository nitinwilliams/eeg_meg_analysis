clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% Initialisation

numparcs=(148/2);
pctval=80;
gamma=2;

%% Parcel stability summary - right hemisphere

load(['Pstabmat_binary_Louvain_SEEG_FC_9.99Hz_righthem_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval) '.mat'],'Pstabmat','Pstabmat_rand');  

Pstabmat_thr=squeeze(prctile(Pstabmat_rand,95,2));
MASK=Pstabmat>Pstabmat_thr;
numstab_vec=(sum(MASK)/numparcs)*100;

save(['B:\Nitin\SCFC_methods\data\seeg_data\numstableparcs_binary_SEEG_FC_gamma' num2str(gamma) '_parc2k9_Louvain_67subs_pct' int2str(pctval) '_righthem.mat'],'numstab_vec');

display(mean(numstab_vec));
