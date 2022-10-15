clear;

%% Initialisation

numparcs=(148/2);
Nfvec={'2.50','3.53','5.00','7.07','9.99','14.14','19.99','28.28','39.99','56.56','79.99','113.13','134.54','159.99','190.27','226.27','269.08','319.99'};
gammavec=[0.8:0.1:5];
thr=0;

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

numc=100;
lefthem_idxs=[1:2:numparcs*2];
righthem_idxs=[2:2:numparcs*2];

pctval=80;
N=numparcs;
perms=100;

UTMASK=triu(ones(numparcs),1);
validxs=find(UTMASK>0);
 
for gamidx=1:length(gammavec),

Pstabmat=zeros(N,length(Nfvec),100);
Pstabmat_rand=zeros(N,length(Nfvec),100,perms);

load(['Sb_gamma' num2str(gammavec(gamidx)) '_pctval' int2str(pctval) '_thr0_67subs_lefthem_wTFC_parc2k9.mat'],'Sb');
gamma=gammavec(gamidx);

for fidx=1:length(Nfvec),
    
    Nf=Nfvec(fidx);
    Sf=Sb(:,fidx);
    fname=strcat('PLV_surrmats_rdtmask_',Nf{1,1},'Hz_thr0_67subs_parc2k9_wgroupmat_1Kreps.mat');
    load(fname,'PLVmats_groupmat');
    numreps=size(PLVmats_groupmat,1);
       
    for ridx=1:numreps,
    
        A=cell(numc,1);
        D=squeeze(PLVmats_groupmat{ridx,1}(lefthem_idxs,lefthem_idxs));
        D(1:N+1:end)=0;
    
        % thresholding networks
   
        tval=prctile(D(validxs),pctval);
        D(D<tval)=0;
    
        D_nc=fill_missingvals(D,numc);
        A(1:numc,1)=D_nc;
    
        %% Determining community structure 

        BF=zeros(N,N,numc);

        for ncidx=1:numc

            Z=A{ncidx};   
    
            M  = 1:N;                    
            Q0 = -1; Q1 = 0;            
            while Q1-Q0>1e-5           
                Q0 = Q1;                
                [M, Q1] = community_louvain(Z,gamma,M);
            end
      
            BF(:,:,ncidx)=clustersol_representation(M);
      
        end

        % final partition

        Z=mean(BF,3);
        Z(1:N+1:end)=0;
    
        M  = 1:N;                    
        Q0 = -1; Q1 = 0;            
        while Q1-Q0>1e-5;           
          Q0 = Q1;                
          [M, Q1] = community_louvain(Z,gamma,M);
        end  
    
        Sbmat=clustersol_representation(M);
    
        for parc_idx=1:length(Sf), 
            
            parc_indices=setdiff([1:length(Sf)]',parc_idx);
            Sr=(Sf(parc_indices)==Sf(parc_idx));           
            Sm=Sbmat(parc_indices,parc_idx);
            
            TP=sum(Sr.*Sm);
            TN=sum(~Sr.*~Sm);
            FP=sum(~Sr.*Sm);
            FN=sum(Sr.*~Sm);
                   
            parcel_stability=(TP+TN)/(TP+TN+FP+FN);
            Pstabmat(parc_idx,fidx,ridx)=parcel_stability;
            
            parcel_stability_rand=zeros(perms,1);
            
            for permidx=1:perms,
            
                 rSm=Sm(randperm(N-1));
                 
                 TP=sum(Sr.*rSm);
                 TN=sum(~Sr.*~rSm);
                 FP=sum(~Sr.*rSm);
                 FN=sum(Sr.*~rSm);
        
                 parcel_stability_rand(permidx)=(TP+TN)/(TP+TN+FP+FN);
                
            end
        
           Pstabmat_rand(parc_idx,fidx,ridx,:)=parcel_stability_rand;
            
        end    
          
    end
           
end

save(['C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\submission\Pstabmat_accu_x1a_origCS_singleslice_lefthem_parc2k9_67subs_thr0_gamma' num2str(gammavec(gamidx)) '_pct' int2str(pctval) '.mat'],'Pstabmat','Pstabmat_rand');

disp(gamidx);

end
