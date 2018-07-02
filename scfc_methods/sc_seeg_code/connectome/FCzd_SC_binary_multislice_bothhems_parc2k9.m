clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% Initialisation

N=74;
lefthem_idxs=[1:2:N*2];
righthem_idxs=[2:2:N*2];

omega_vals=cell(1,3);
gamma=2;
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

%% FC matrix

% Partition similarity

A1=cell(1,2);

load(['B:\Nitin\SCFC_methods\data\seeg_data\Sb_binary_gamma' num2str(gamma) '_pctval80_thr0_' Nfvec{1,1} 'Hz_lefthem_wTFC_parc2k9.mat'],'Sb','Zmat'); 
LH=Sb;
A1{1}=Zmat{1,1};

load(['B:\Nitin\SCFC_methods\data\seeg_data\Sb_binary_gamma' num2str(gamma) '_pctval80_thr0_' Nfvec{1,1} 'Hz_righthem_wTFC_parc2k9.mat'],'Sb','Zmat'); 
RH=Sb;
A1{2}=Zmat{1,1};

partsim_val=partsim(LH,RH);

% finding omega value

omegavec=[0:0.1:3];
match_vec=zeros(1,length(omegavec));

for omidx=1:length(omegavec),     
    omval=omegavec(omidx);
    [S,~]=multislice_communitydetection(A1,N,gamma,omval,0,'nonullmodel');
    match_vec(1,omidx)=(sum(S(:,1)-S(:,2))==0);
end

cind=find(match_vec,1,'first');
if isempty(cind), cind=length(omegavec); end
omega_vals{1}=partsim_val*omegavec(cind);

% Multi-slice community detection

[S1,~]=multislice_communitydetection(A1,N,gamma,omega_vals{1},0,'nonullmodel');

%% SC matrix

% Partition similarity

A2=cell(1,2);

load(['B:\Nitin\SCFC_methods\data\sc_data\Sb_binary_SC_gamma' num2str(gamma) '_pctval80_thr0_lefthem_wTFC_parc2k9.mat'],'Sb','Zmat'); 
LH=Sb;
A2{1}=Zmat{1,1};

load(['B:\Nitin\SCFC_methods\data\sc_data\Sb_binary_SC_gamma' num2str(gamma) '_pctval80_thr0_righthem_wTFC_parc2k9.mat'],'Sb','Zmat'); 
RH=Sb;
A2{2}=Zmat{1,1};

partsim_val=partsim(LH,RH);

% finding omega value

omegavec=[0:0.1:3];
match_vec=zeros(1,length(omegavec));

for omidx=1:length(omegavec),     
    omval=omegavec(omidx);
    [S,~]=multislice_communitydetection(A2,N,gamma,omval,0,'nonullmodel');
    match_vec(1,omidx)=(sum(S(:,1)-S(:,2))==0);
end

cind=find(match_vec,1,'first');
if isempty(cind), cind=length(omegavec); end
omega_vals{2}=partsim_val*omegavec(cind);

% Multi-slice community detection

[S2,~]=multislice_communitydetection(A2,N,gamma,omega_vals{2},0,'nonullmodel');

%% SC - FC matrix

FC=zeros(2*N,1);
SC=zeros(2*N,1);

FC(lefthem_idxs)=S1(:,1);
FC(righthem_idxs)=S1(:,2);

SC(lefthem_idxs)=S2(:,1);
SC(righthem_idxs)=S2(:,2);

partsim_val=partsim(SC,FC);

% finding omega value
 
omegavec=[0:0.1:3];
match_vec=zeros(1,length(omegavec));

A3=cell(1,4);
A3{1}=A1{1};
A3{2}=A1{2};
A3{3}=A2{1};
A3{4}=A2{2};

for omidx=1:length(omegavec),     
    omega_vals{3}=omegavec(omidx);
    [S,~]=multislice_communitydetection_ext(A3,N,gamma,omega_vals);
    match_vec(1,omidx)=(sum(S(:,1)-S(:,2))==0);
end

cind=find(match_vec,1,'first');
if isempty(cind), cind=length(omegavec); end
omega_vals{3}=partsim_val*omegavec(cind);

% Multi-slice community detection

[S3,~]=multislice_communitydetection_ext(A3,N,gamma,omega_vals);

save(['B:\Nitin\SCFC_methods\data\sc_seeg_data\Sb_binary_SC_SEEG_FC_parc2k9_mlcd_gamma' num2str(gamma) '_optomegavals.mat'],'S3','A3','omega_vals');
