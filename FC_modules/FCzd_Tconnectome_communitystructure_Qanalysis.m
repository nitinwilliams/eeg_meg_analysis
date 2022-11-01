clear;
rmpath(genpath('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\'));
addpath(genpath('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\'));
addpath(genpath('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_figs\'));
addpath(genpath('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\BCT\'));
addpath('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_code\GA_code\Connectome_estimation\FC_code\submission\');

%% Initialisation

numparcs=(148/2);
Nfvec={'2.50','3.53','5.00','7.07','9.99','14.14','19.99','28.28','39.99','56.56','79.99','113.13','134.54','159.99','190.27','226.27','269.08','319.99'};

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

%% Estimating strengths

numc=100;
pctval=80;
nct=0;
AF=cell(numc,length(Nfvec));
lefthem_idxs=[1:2:numparcs*2];
righthem_idxs=[2:2:numparcs*2];

UTMASK=triu(ones(numparcs),1);
validxs=find(UTMASK>0);

for fidx=1:length(Nfvec),
    
   Nf=Nfvec(fidx);
   fname=strcat('PLVmats_rdtmask_',Nf{1,1},'Hz_thr',int2str(nct),'_67subs_parc2k9_wgroupmat.mat');
   load(fname,'PLVmats_groupmat');
   D=PLVmats_groupmat(lefthem_idxs,lefthem_idxs);
   D(1:numparcs+1:end)=0;
    
   % thresholding networks
   
   tval=prctile(D(validxs),pctval);
   D(D<tval)=0;
   
   D_nc=fill_missingvals(D,numc);
   AF(1:numc,fidx)=D_nc;
       
end

%% Determining stability - community structure

reps=100;
gammavec=[0.8:0.1:5];
QMAT=zeros(length(Nfvec),length(gammavec));
NMAT=zeros(length(Nfvec),length(gammavec));
zQMAT=zeros(length(Nfvec),length(gammavec));
mQMAT=zeros(length(Nfvec),length(gammavec));
QMATND=zeros(length(Nfvec),length(gammavec),reps);

for freqidx=1:length(Nfvec),

    A=AF(:,freqidx);
    QVEC=zeros(1,length(gammavec));

    parfor gamidx=1:length(gammavec),
        
        gamma=gammavec(gamidx);
        [QVEC(gamidx),NMAT(freqidx,gamidx)]=partitionquality_ss(A,gamma);
          
    end

    %% determining partition quality null distribution 

    for repidx=1:reps
    
        A_ND=randomise_und_ensemble(A); 

        parfor gamidx=1:length(gammavec)
        
            gamma=gammavec(gamidx);
            QMATND(freqidx,gamidx,repidx)=partitionquality_ss(A_ND,gamma);
            
        end
           
    end

    %% z-scoring stability

    ND=squeeze(QMATND(freqidx,:,:))';

    QMAT(freqidx,:) = QVEC;
    zQMAT(freqidx,:)=(QVEC-mean(ND))./std(ND);
    mQMAT(freqidx,:)=(QVEC-mean(ND));

    disp(freqidx);

end

save(['C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\submission\parc2k9_67subs_pct' num2str(pctval) '_Qanalysis_lefthem.mat'],'NMAT','QMAT','zQMAT','mQMAT');

%% Plotting

figure;
imagesc(zQMAT,[0,20]);
colorbar;
xlabel('Gamma');
ylabel('Frequency (Hz)');
set(gca,'YDir','normal','XTick',[1:4:length(gammavec)],'XTickLabel',gammavec(1:4:end),'YTick',[1:length(Nfvec)],'YTickLabel',Nfvec);
set(gcf,'Position',[300,300,1000,750]);
saveas(gcf,'C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_figs\zQMAT_parc2k9_67subs_Qanalysis_lefthem','png');
