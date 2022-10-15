clear;

%% Initialisation

numparcs=(148/2);
nct=0;
pctval=80;
Nfvec=[3,4,5,7,10,14,20,28,40,57,80,113,135,160,190,226,269,320];
gammavec=[0.8:0.1:5];

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

%% Parcel stability summary - left hemisphere

Pstab_summary=zeros(length(Nfvec),length(gammavec));
Pstab_chance_summary=zeros(length(Nfvec),length(gammavec),100);

for gamidx=1:length(gammavec),

load(['Pstabmat_accu_x1a_origCS_singleslice_lefthem_parc2k9_67subs_thr' int2str(nct) '_gamma' num2str(gammavec(gamidx)) '_pct' int2str(pctval) '.mat'],'Pstabmat','Pstabmat_rand');  

Pstabmat=mean(Pstabmat,3);
Pstabmat_rand=mean(Pstabmat_rand,4);
Pstabmat_thr=prctile(Pstabmat_rand,95,3);

MASK1=Pstabmat>Pstabmat_thr;
Pstab_summary(:,gamidx)=(sum(MASK1)/numparcs)*100;

MASK2=Pstabmat_rand>repmat(Pstabmat_thr,[1,1,size(Pstabmat_rand,3)]);
Pstab_chance_summary(:,gamidx,:)=ceil((sum(MASK2,1)/numparcs)*100);

disp(gamidx);

end

save('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\submission\parc2k9_67subs_pct80_parcelstability_lefthem.mat','Pstab_summary','Pstab_chance_summary');

% plotting

figure;
imagesc(Pstab_summary,[0,100]);
colorbar;
xlabel('Gamma value');
ylabel('Frequency (Hz)');
title('Left hemisphere - percentage of stable parcels');

set(gca,'YDir','normal','YTick',[1:1:length(Nfvec)],'YTickLabel',Nfvec(1:1:end),'XTick',[1:4:length(gammavec)],'XTickLabel',gammavec(1:4:end));

%% Parcel stability summary - right hemisphere

Pstab_summary=zeros(length(Nfvec),length(gammavec));
Pstab_chance_summary=zeros(length(Nfvec),length(gammavec),100);

for gamidx=1:length(gammavec),

load(['Pstabmat_accu_x1a_origCS_singleslice_righthem_parc2k9_67subs_thr' int2str(nct) '_gamma' num2str(gammavec(gamidx)) '_pct' int2str(pctval) '.mat'],'Pstabmat','Pstabmat_rand');  

Pstabmat=mean(Pstabmat,3);
Pstabmat_rand=mean(Pstabmat_rand,4);
Pstabmat_thr=prctile(Pstabmat_rand,95,3);

MASK1=Pstabmat>Pstabmat_thr;
Pstab_summary(:,gamidx)=(sum(MASK1)/numparcs)*100;

MASK2=Pstabmat_rand>repmat(Pstabmat_thr,[1,1,size(Pstabmat_rand,3)]);
Pstab_chance_summary(:,gamidx,:)=ceil((sum(MASK2,1)/numparcs)*100);

disp(gamidx);

end

save('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\submission\parc2k9_67subs_pct80_parcelstability_righthem.mat','Pstab_summary','Pstab_chance_summary');

% plotting

figure;
imagesc(Pstab_summary,[0,100]);
colorbar;
xlabel('Gamma value');
ylabel('Frequency (Hz)');
title('Right hemisphere - percentage of stable parcels');

set(gca,'YDir','normal','YTick',[1:1:length(Nfvec)],'YTickLabel',Nfvec(1:1:end),'XTick',[1:4:length(gammavec)],'XTickLabel',gammavec(1:4:end));
