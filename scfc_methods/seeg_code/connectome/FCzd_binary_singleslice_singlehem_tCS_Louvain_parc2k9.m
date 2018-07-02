clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% Initialisation

X=74;
lefthem_idxs=[1:2:X*2];
righthem_idxs=[2:2:X*2];

gammavec=[2];
Nfvec={'9.99'};

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

%% Fill missing values

numc=5000;
pctval=80;
nct=0;
A=cell(numc,length(Nfvec));

numparcs=length(righthem_idxs);

UTMASK=triu(ones(numparcs),1);
validxs=find(UTMASK>0);

for gamidx=1:length(gammavec),

for fidx=1:length(Nfvec),
    
    Nf=Nfvec(fidx);
    fname=strcat('PLVmats_rdtmask_',Nf{1,1},'Hz_thr0_67_PARC2k9_wgroupmat.mat');
    load(fname,'PLVmats_groupmat');
    D=PLVmats_groupmat(righthem_idxs,righthem_idxs);
    D(1:numparcs+1:end)=0;
       
    % thresholding networks
   
    tval=prctile(D(validxs),pctval);
    D(D<tval)=0;
    
    D_nc=fill_missingvals(D,numc);
    A(1:numc,fidx)=D_nc;
    
end

%% Determining community structure 

gamma=gammavec(gamidx);
N=numparcs;
T=length(Nfvec);
BF=zeros(N,N,numc,T);

parfor ncidx=1:numc,

Z=A(ncidx,:);   

for s=1:T,
    
      M  = 1:N;                    
      Q0 = -1; Q1 = 0;            
      while Q1-Q0>1e-5;           
          Q0 = Q1;                
          [M, Q1] = community_louvain(Z{s},gamma,M);
      end
      
    BF(:,:,ncidx,s)=clustersol_representation(M);
      
end

end

% final partition

Sb=zeros(N,T);
Zmat=cell(1,T);

for s=1:T,
   
  Z=mean(squeeze(BF(:,:,:,s)),3);
  
  % thresholding
  
  tval=prctile(Z(validxs),pctval);
  Z(Z<tval)=0;
  Z(Z>=tval)=1;
  
  Z(1:N+1:end)=0;
  
  M  = 1:N;                    
  Q0 = -1; Q1 = 0;            
  while Q1-Q0>1e-5;           
     Q0 = Q1;                
     [M, Q1] = community_louvain(Z,gamma,M);
  end   
  
  Sb(:,s)=M;
  Zmat{1,s}=Z;
  
end

display(gamidx);

save(['B:\Nitin\SCFC_methods\data\seeg_data\Sb_binary_gamma' num2str(gamma) '_pctval' int2str(pctval) '_thr' int2str(nct) '_' Nf{1,1} 'Hz_righthem_wTFC_parc2k9.mat'],'Sb','Zmat');

end