clear;
addpath(genpath('B:\Nitin\HBPconnectome\'));

%% Initialisation

numparcs=148;
load('parc2K9_xycds');
load('CS_RGBcolormap');
H=zeros(numparcs,3);

lefthem_idxs=[1:2:numparcs];
righthem_idxs=[2:2:numparcs];

load(['B:\Nitin\SCFC_methods\data\sc_meg_data\Sb_binary_SC_MEG_FC_parc2k9_mlcd_gamma2_optomegavals.mat']);

LH=S3(1:(numparcs/2),2);
num_mods=max(LH);

for idx=1:num_mods,
    H(lefthem_idxs((LH==idx)),1)=circ_mean(xycds_fl_av_rc_angle(lefthem_idxs(LH==idx)),[],1);
    H(lefthem_idxs((LH==idx)),2)=1;
    H(lefthem_idxs((LH==idx)),3)=mean(xycds_fl_av_rc_normdist(lefthem_idxs(LH==idx)));
end

RH=S3((numparcs/2)+1:end,2);
num_mods=max(RH);

for idx=1:num_mods,
    H(righthem_idxs((RH==idx)),1)=circ_mean(xycds_fl_av_rc_angle(righthem_idxs(RH==idx)),[],1);
    H(righthem_idxs((RH==idx)),2)=1;
    H(righthem_idxs((RH==idx)),3)=mean(xycds_fl_av_rc_normdist(righthem_idxs(RH==idx)));
end

% masking

load(['B:\Nitin\SCFC_methods\data\sc_meg_data\Pstabmat_binary_Louvain_SC_MEG_FC_9.83Hz_parc2k9_67subs_gamma2_pct80.mat'],'Pstabmat','Pstabmat_rand'); 

Pstabmat=mean(Pstabmat,2);
Pstabmat_rand=mean(Pstabmat_rand,3);
Pstabmat_thr=prctile(Pstabmat_rand,95,2);
MASK_mat=reshape(Pstabmat>Pstabmat_thr,[numparcs,2]);

LH_mask=MASK_mat(1:(numparcs/2),2);
RH_mask=MASK_mat((numparcs/2)+1:end,2);

MASK=zeros(numparcs,1);
MASK(lefthem_idxs)=LH_mask;
MASK(righthem_idxs)=RH_mask;

gray_indices=find(MASK==0)-1;
display(gray_indices);

% normalising angle & distance

H(:,1)=(H(:,1)-(-pi))/(pi-(-pi));

%% Colormap generation

M=hsv2rgb(H);
M(gray_indices+1,:)=0.6;
M_scaled=round(255*M);

figure;
hold on;

for idx=1:size(parc2K9_xycds,1),
    plot(parc2K9_xycds(idx,1),parc2K9_xycds(idx,2),'.','MarkerSize',30,'Color',M(idx,:));        
end

set(gca,'YDir','Reverse','XTick','','YTick','');

save(['B:\Nitin\SCFC_methods\data\sc_meg_data\RGBcolormap_binary_Louvain_SC_modules_mlcd_parc2k9_gamma2_optomega.mat'],'M_scaled');
