clear;

%% Loading data

DAT_isochronousdelays=readNPY('Z:\Documents\MEGMOD\Python\data\bolfi_13params_isochronousdelaysmodel_expMEGdata_17102022_9004sims_defthresh_SSIlogdisc_targetprob0pt8_bolfisample.npy');
DAT_distancedependentdelays=readNPY('Z:\Documents\MEGMOD\Python\data\bolfi_12params_distancedependentdelaysmodel_expMEGdata_09092022_9093sims_defthresh_SSIlogdisc_targetprob0pt8_bolfisample.npy');
DAT_mixeddelays=readNPY('Z:\Documents\MEGMOD\Python\data\bolfi_14params_mixeddelaysmodel_expMEGdata_21102022_9063sims_defthresh_SSIlogdisc_targetprob0pt8_bolfisample.npy');

%% Re-ordering and rescaling

N=size(DAT_isochronousdelays,1);

origorder_isochronousdelays=[10,11,12,13,2,3,8,9,6,1,7,4,5];
origorder_distancedependentdelays=[9,10,11,12,2,3,6,7,4,1,5,8];
origorder_mixeddelays=[11,12,13,14,3,4,8,9,6,1,7,5,10,2];

DAT_isochronousdelays_origorder=DAT_isochronousdelays(:,origorder_isochronousdelays);
DAT_distancedependentdelays_origorder=DAT_distancedependentdelays(:,origorder_distancedependentdelays);
DAT_mixeddelays_origorder=DAT_mixeddelays(:,origorder_mixeddelays);

meanvals_isochronousdelays=[20,18,18,1,3,5,18.6,15.1,1.5,2.5,0.15,10,0.2];
stdvals_isochronousdelays=[5,6,6,0.2,1,1,3.8,4.7,0.5,0.5,0.05,3,0.05];

meanvals_distancedependentdelays=[20,18,18,1,3,5,18.6,15.1,1.5,2.5,0.15,8];
stdvals_distancedependentdelays=[5,6,6,0.2,1,1,3.8,4.7,0.5,0.5,0.05,2];

meanvals_mixeddelays=[20,18,18,1,3,5,18.6,15.1,1.5,2.5,0.15,10,8,0.5];
stdvals_mixeddelays=[5,6,6,0.2,1,1,3.8,4.7,0.5,0.5,0.05,3,2,0.15];

WCnet_noisy_isochronousdelays_posteriordistribution=(DAT_isochronousdelays_origorder.*repmat(stdvals_isochronousdelays,[N,1]))+repmat(meanvals_isochronousdelays,[N,1]);
WCnet_noisy_distancedependentdelays_posteriordistribution=(DAT_distancedependentdelays_origorder.*repmat(stdvals_distancedependentdelays,[N,1]))+repmat(meanvals_distancedependentdelays,[N,1]);
WCnet_noisy_mixeddelays_posteriordistribution=(DAT_mixeddelays_origorder.*repmat(stdvals_mixeddelays,[N,1]))+repmat(meanvals_mixeddelays,[N,1]);

%% Saving

save('Z:\Documents\MEGMOD\MATLAB\data\WCnet_noisy_posteriordistributions_2ksims.mat','WCnet_noisy_isochronousdelays_posteriordistribution',...
    'WCnet_noisy_distancedependentdelays_posteriordistribution','WCnet_noisy_mixeddelays_posteriordistribution');

%% Marginal distributions

paramnames_isochronousdelays={'w_{ee}','w_{ei}','w_{ie}','w_{ii}','b_{e}','b_{i}','\tau_{e}','\tau_{i}','k','IH_{scaling}','\psi_{sigma}','delay','coeffvar_{delay}'};
paramnames_distancedependentdelays={'w_{ee}','w_{ei}','w_{ie}','w_{ii}','b_{e}','b_{i}','\tau_{e}','\tau_{i}','k','IH_{scaling}','\psi_{sigma}','velocity'};
paramnames_mixeddelays={'w_{ee}','w_{ei}','w_{ie}','w_{ii}','b_{e}','b_{i}','\tau_{e}','\tau_{i}','k','IH_{scaling}','\psi_{sigma}','delay','velocity','coeff_{balance}'};

disp(median(WCnet_noisy_isochronousdelays_posteriordistribution));
disp(median(WCnet_noisy_distancedependentdelays_posteriordistribution));
disp(median(WCnet_noisy_mixeddelays_posteriordistribution));

% Isochronous delays

numbins_marginals=30;
figure('units','normalized','outerposition',[0 0 1 1])

for col_idx=1:size(WCnet_noisy_isochronousdelays_posteriordistribution,2)

subplot(3,5,col_idx)
histogram(WCnet_noisy_isochronousdelays_posteriordistribution(:,col_idx),numbins_marginals);
xlim([meanvals_isochronousdelays(col_idx)-2*stdvals_isochronousdelays(col_idx),meanvals_isochronousdelays(col_idx)+2*stdvals_isochronousdelays(col_idx)]);
ylim([0,200]);
axis square;
title(paramnames_isochronousdelays{col_idx});

end

saveas(gcf,'Z:\Documents\MEGMOD\MATLAB\figs\marginaldists_13params_WCnet_noisy_isochronousdelays.png');

% Distance-dependent delays

numbins_marginals=30;
figure('units','normalized','outerposition',[0 0 1 1])

for col_idx=1:size(WCnet_noisy_distancedependentdelays_posteriordistribution,2)

subplot(3,5,col_idx)
histogram(WCnet_noisy_distancedependentdelays_posteriordistribution(:,col_idx),numbins_marginals);
xlim([meanvals_distancedependentdelays(col_idx)-2*stdvals_distancedependentdelays(col_idx),meanvals_distancedependentdelays(col_idx)+2*stdvals_distancedependentdelays(col_idx)]);
ylim([0,200]);
axis square;
title(paramnames_distancedependentdelays{col_idx});

end

saveas(gcf,'Z:\Documents\MEGMOD\MATLAB\figs\marginaldists_12params_WCnet_noisy_distancedependentdelays.png');

% Mixed delays

numbins_marginals=30;
figure('units','normalized','outerposition',[0 0 1 1])

for col_idx=1:size(WCnet_noisy_mixeddelays_posteriordistribution,2)

subplot(3,5,col_idx)
histogram(WCnet_noisy_mixeddelays_posteriordistribution(:,col_idx),numbins_marginals);
xlim([meanvals_mixeddelays(col_idx)-2*stdvals_mixeddelays(col_idx),meanvals_mixeddelays(col_idx)+2*stdvals_mixeddelays(col_idx)]);
ylim([0,200]);
axis square;
title(paramnames_mixeddelays{col_idx});

end

saveas(gcf,'Z:\Documents\MEGMOD\MATLAB\figs\marginaldists_14params_WCnet_noisy_mixeddelays.png');

%% Conditional distributions

% Isochronous delays

[RHO,PVALS]=corr(WCnet_noisy_isochronousdelays_posteriordistribution);
bonf_fac=((length(paramnames_isochronousdelays))*(length(paramnames_isochronousdelays)-1))/2;

figure;
imagesc(RHO);
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',[1:length(paramnames_isochronousdelays)],'XTickLabel',paramnames_isochronousdelays,...
    'YTick',[1:length(paramnames_isochronousdelays)],'YTickLabel',paramnames_isochronousdelays);

figure;
imagesc(PVALS<(0.05/bonf_fac));
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',[1:length(paramnames_isochronousdelays)],'XTickLabel',paramnames_isochronousdelays,...
    'YTick',[1:length(paramnames_isochronousdelays)],'YTickLabel',paramnames_isochronousdelays);

% Distance-dependent delays

[RHO,PVALS]=corr(WCnet_noisy_distancedependentdelays_posteriordistribution);
bonf_fac=((length(paramnames_distancedependentdelays))*(length(paramnames_distancedependentdelays)-1))/2;

figure;
imagesc(RHO);
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',[1:length(paramnames_distancedependentdelays)],'XTickLabel',paramnames_distancedependentdelays,...
    'YTick',[1:length(paramnames_distancedependentdelays)],'YTickLabel',paramnames_distancedependentdelays);

figure;
imagesc(PVALS<(0.05/bonf_fac));
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',[1:length(paramnames_distancedependentdelays)],'XTickLabel',paramnames_distancedependentdelays,...
    'YTick',[1:length(paramnames_distancedependentdelays)],'YTickLabel',paramnames_distancedependentdelays);

% Mixed delays

[RHO,PVALS]=corr(WCnet_noisy_mixeddelays_posteriordistribution);
bonf_fac=((length(paramnames_mixeddelays))*(length(paramnames_mixeddelays)-1))/2;

figure;
imagesc(RHO);
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',[1:length(paramnames_mixeddelays)],'XTickLabel',paramnames_mixeddelays,...
    'YTick',[1:length(paramnames_mixeddelays)],'YTickLabel',paramnames_mixeddelays);

figure;
imagesc(PVALS<(0.05/bonf_fac));
colorbar;
axis square;
set(gca,'YDir','Normal','XTick',[1:length(paramnames_mixeddelays)],'XTickLabel',paramnames_mixeddelays,...
    'YTick',[1:length(paramnames_mixeddelays)],'YTickLabel',paramnames_mixeddelays);
