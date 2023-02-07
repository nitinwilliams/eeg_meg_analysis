function generateGlobals_allsubs()

subject_sets = {'S0001','S0003','S0005','S0006','S0007','S0008','S0011','S0026','S0034','S0035','S0036','S0038',...
    'S0039','S0049','S0113','S0115','S0116','S0117','S0118','S0119','S0120','S0121','S0122','S0123','S0124','S0125','S0126',...
    'S0127','S0128','S0129','S0130','S0132','S0136','S0137','S0140','S0141','S0142','S0143','S0144','S0145','S0146','S0147',...
    'S0148','S0149','S0150','S0151','S0152','S0153','S0155','S0156','S0157','S0161','S0162','S0166','S0167','S0168','S0169',...
    'S0170','S0171','S0172','S0173','S0174','S0175','S0177','S0179','S0180','S0181','S0184','S0185','S0186','S0188','S0189',...
    'S0193','S0194','S0195'};

global PY2LV_transvec
load('/scratch/work/willian1/MEGMOD/MATLAB/data/PY2LV_transvec.mat','PY2LV_transvec');
global avg_sc
load('/scratch/work/willian1/MEGMOD/MATLAB/data/SC_parc2k9_57HCP.mat','avg_sc');
% global WCnet_noisy_mixeddelays_paramcombinations
% load('/scratch/work/willian1/MEGMOD/MATLAB/data/param_combinations_10ksims.mat','WCnet_noisy_mixeddelays_paramcombinations');
% global WCnet_noisy_isochronousdelays_paramcombinations
% load('/scratch/work/willian1/MEGMOD/MATLAB/data/param_combinations_10ksims.mat','WCnet_noisy_isochronousdelays_paramcombinations');
% global WCnet_noisy_variabledelays_paramcombinations
% load('/scratch/work/willian1/MEGMOD/MATLAB/data/param_combinations_10ksims.mat','WCnet_noisy_variabledelays_paramcombinations');
global WCnet_noisy_isochronousdelays_posteriordistribution
load('/scratch/work/willian1/MEGMOD/MATLAB/data/WCnet_noisy_posteriordistributions_2ksims.mat','WCnet_noisy_isochronousdelays_posteriordistribution');
global WCnet_noisy_distancedependentdelays_posteriordistribution
load('/scratch/work/willian1/MEGMOD/MATLAB/data/WCnet_noisy_posteriordistributions_2ksims.mat','WCnet_noisy_distancedependentdelays_posteriordistribution');
global WCnet_noisy_mixeddelays_posteriordistribution
load('/scratch/work/willian1/MEGMOD/MATLAB/data/WCnet_noisy_posteriordistributions_2ksims.mat','WCnet_noisy_mixeddelays_posteriordistribution')

global source2parc fwd_mat inv_mat

source2parc=cell(length(subject_sets),1);
fwd_mat=cell(length(subject_sets),1);
inv_mat=cell(length(subject_sets),1);

for sub_idx=1:length(subject_sets)
fname = char(strcat('/scratch/work/willian1/MEGMOD/MATLAB/data/MEG_subjectwise_data/sourceIDs/',subject_sets(sub_idx), '_set01__parc2009.csv'));
source2parc{sub_idx}=table2array(readtable(fname,'ReadVariableNames',0,'Format','auto'));    
fname = char(strcat('/scratch/work/willian1/MEGMOD/MATLAB/data/MEG_subjectwise_data/fixed_fwd_ops/',subject_sets(sub_idx), '_set01__py-fwd.csv'));
fwd_mat{sub_idx}=table2array(readtable(fname,'ReadVariableNames',0,'Format','auto'));
fname = char(strcat('/scratch/work/willian1/MEGMOD/MATLAB/data/MEG_subjectwise_data/fidweighted_inv_ops/',subject_sets(sub_idx), '_set01__parc2009.csv'));
inv_mat{sub_idx}=table2array(readtable(fname,'ReadVariableNames',0,'Format','auto'));
end

end

