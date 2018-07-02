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

%% Fill missing values

Nf={'9.99'};

numparcs=74;
lefthem_idxs=[1:2:numparcs*2];
righthem_idxs=[2:2:numparcs*2];

pctval=80;
N=numparcs;
perms=100;
gamma=2;
nct=0;

Pstabmat=zeros(N,100);
Pstabmat_rand=zeros(N,100,perms);
UTMASK=triu(ones(numparcs),1);
validxs=find(UTMASK>0);

load(['Sb_weighted_gamma' num2str(gamma) '_pctval' int2str(pctval) '_thr0_' Nf{1,1} 'Hz_righthem_wTFC_parc2k9.mat'],'Sb');
Sf=Sb;

load(['Sb_weighted_surr_gamma' num2str(gamma) '_pctval' int2str(pctval) '_thr' int2str(nct) '_' Nf{1,1} 'Hz_righthem_wTFC_parc2k9.mat'],'Sb');
Sb_surr=Sb;

numreps=size(Sb_surr,1);
       
for ridx=1:numreps,
     
   Sbmat=clustersol_representation(Sb_surr{ridx,1});
    
   for parc_idx=1:length(Sf),
            
         parc_indices=setdiff([1:length(Sf)]',parc_idx);
         Sr=(Sf(parc_indices)==Sf(parc_idx));           
         Sm=Sbmat(parc_indices,parc_idx);
                        
            TP=sum(Sr.*Sm);
            TN=sum(~Sr.*~Sm);
            FP=sum(~Sr.*Sm);
            FN=sum(Sr.*~Sm);
                   
            parcel_stability=(TP+TN)/(TP+TN+FP+FN);
            Pstabmat(parc_idx,ridx)=parcel_stability;
            
            parcel_stability_rand=zeros(perms,1);
            
            for permidx=1:perms,
            
                rSm=Sm(randperm(N-1));
                 
                TP=sum(Sr.*rSm);
                TN=sum(~Sr.*~rSm);
                FP=sum(~Sr.*rSm);
                FN=sum(Sr.*~rSm);
        
                parcel_stability_rand(permidx)=(TP+TN)/(TP+TN+FP+FN);
                
            end
        
           Pstabmat_rand(parc_idx,ridx,:)=parcel_stability_rand;
            
   end 

   display(ridx); 
   
end
          
save(['B:\Nitin\SCFC_methods\data\seeg_data\Pstabmat_weighted_Louvain_SEEG_FC_' Nf{1,1} 'Hz_righthem_parc2k9_67subs_gamma' num2str(gamma) '_pct' int2str(pctval) '.mat'],'Pstabmat','Pstabmat_rand');
