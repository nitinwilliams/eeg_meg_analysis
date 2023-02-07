clear;

%% Load data

N=148;
UTMAT=triu(ones(N),1);
ut_idxs=find(UTMAT(:));

lefthem_idxs=[1:2:N];
righthem_idxs=[2:2:N];

% Load simulation results
reps_per_batch=5;
batch_indices=[1:2000];

wPLI_FC_fullstrengths=cell(length(batch_indices),1);
param_fullmat=zeros(max(batch_indices)*reps_per_batch,14);
for batch_idx=1:length(batch_indices)
    fname=['S:\work\willian1\MEGMOD\MATLAB\data\WCnet_mixeddelays_noisy_60s_simulations_batch',int2str(batch_indices(batch_idx)),'_run1.mat'];
    if isfile(fname)
    load(fname,'wPLImags_acrossreps','mod_params');
    param_fullmat((batch_idx-1)*reps_per_batch+1:(batch_idx*reps_per_batch),:)=mod_params;
    wPLI_FC_fullstrengths{batch_idx}=wPLImags_acrossreps;
    end
    disp(batch_idx);
end

wPLI_FC_fullstrengths_cellarray=cell(max(batch_indices)*reps_per_batch,1);
for batch_idx=1:length(batch_indices)
    if ~isempty(wPLI_FC_fullstrengths{batch_idx})
    wPLI_FC_fullstrengths_cellarray((batch_idx-1)*reps_per_batch+1:(batch_idx*reps_per_batch),1)=num2cell(wPLI_FC_fullstrengths{batch_idx,1},[1,5])';
    end
end

% Cleaning data
val_idxs=[];
for idx=1:size(wPLI_FC_fullstrengths_cellarray,1)
    if ~isempty(wPLI_FC_fullstrengths_cellarray{idx}) && sum(wPLI_FC_fullstrengths_cellarray{idx})>0
        val_idxs(end+1)=idx;
    end
end

wPLI_FC_fullstrengths_cellarray=wPLI_FC_fullstrengths_cellarray(val_idxs);
param_fullmat=param_fullmat(val_idxs,:);

%% Estimate SSI

% Load reference data

fname=strcat('Z:\Documents\MEGMOD\MATLAB\data\wPLImats_9.83Hz_75subs_PARC2k9_148parcs_wgroupmat.mat');
load(fname,'wPLImats_strength_groupmat');
wPLImags_connstrengths_mat=wPLImats_strength_groupmat([lefthem_idxs,righthem_idxs],[lefthem_idxs,righthem_idxs]);
wPLImags_connstrengths_ref=wPLImags_connstrengths_mat(ut_idxs);

save('Z:\Documents\MEGMOD\MATLAB\data\WCnetnoisy_mixeddelays_GPfitting_expMEGdata_75subs_14params_connstrengths_refdata.mat','wPLImags_connstrengths_ref');

% Estimate SSI values

SSI_luminancemap=zeros(size(param_fullmat,1),length(ut_idxs));
SSI_contrastmap=zeros(size(param_fullmat,1),length(ut_idxs));
SSI_structuremap=zeros(size(param_fullmat,1),length(ut_idxs));

for row_idx=1:size(param_fullmat,1)

    wPLImags_connstrengths_est=wPLI_FC_fullstrengths_cellarray{row_idx};
           
    [~,SSI_luminancemap(row_idx,:)]=ssim(wPLImags_connstrengths_ref,wPLImags_connstrengths_est,'Exponents',[1,0,0],'DataFormat',"SS",'DynamicRange',1,'Radius',1.5);
    [~,SSI_contrastmap(row_idx,:)]=ssim(wPLImags_connstrengths_ref,wPLImags_connstrengths_est,'Exponents',[0,1,0],'DataFormat',"SS",'DynamicRange',1,'Radius',1.5);
    [~,SSI_structuremap(row_idx,:)]=ssim(wPLImags_connstrengths_ref,wPLImags_connstrengths_est,'Exponents',[0,0,1],'DataFormat',"SS",'DynamicRange',1,'Radius',1.5);

    disp(row_idx);

end

% Plot images

figure;
imagesc(SSI_luminancemap);
colorbar;
axis square;
set(gca,'YDir','Normal');

figure;
imagesc(SSI_contrastmap);
colorbar;
axis square;
set(gca,'YDir','Normal');

figure;
imagesc(SSI_structuremap);
colorbar;
axis square;
set(gca,'YDir','Normal');

% Plot histogram

figure;
histogram(mean(SSI_luminancemap,2),[0:0.02:1]);
axis square;

figure;
histogram(mean(SSI_contrastmap,2),[0:0.02:1]);
axis square;

figure;
histogram(mean(SSI_structuremap,2),[0:0.02:1]);
axis square;

%% Fit GP model

exponent_vec=[1,1,1];
SSI_combinedmap=((SSI_luminancemap.^exponent_vec(1)).*(SSI_contrastmap.^exponent_vec(2)).*(SSI_structuremap.^exponent_vec(3)));
SSI_combined=mean(SSI_combinedmap,2);

% Plotting combined map

figure;
imagesc(SSI_combinedmap);
colorbar;
axis square;
set(gca,'YDir','Normal');

figure;
histogram(SSI_combined);
axis square;

% GPR fitting

X=param_fullmat;
Y=log(1-SSI_combined);

figure;
histogram(Y);
axis square;

num_reps=1;
LL_best=-Inf;
for rep_idx=1:num_reps
GPMOD_current=fitrgp(X,Y,'KernelFunction','ardsquaredexponential','BasisFunction','constant','Standardize',1);
if GPMOD_current.LogLikelihood>LL_best
    LL_best=GPMOD_current.LogLikelihood;
    GPMOD=GPMOD_current;
end
disp(rep_idx);
end

Y_hat = resubPredict(GPMOD);

figure;
scatter(Y,Y_hat,'o');
xlim([-2,0]);
ylim([-2,0]);
axis square;
h=lsline;
h.Color='r';
h.LineWidth=5;
xlabel('Actual discrepancy');
ylabel('Predicted discrepancy');

set(gca,'FontSize',18);

disp(corr(Y,Y_hat));

% Parameter weights

sigmaL = GPMOD.KernelInformation.KernelParameters(1:end-1);     % Learned length scales
weights = exp(-sigmaL);                                         % Predictor weights
weights = weights/sum(weights);

figure;
bar(weights);

%% Saving data

params_muvec=[20,18,18,1,3,5,18.6,15.1,1.5,2.5,0.15,10,8,0.5];
params_sigmavec=[5,6,6,0.2,1,1,3.8,4.7,0.5,0.5,0.05,3,2,0.15];

N=length(Y);
X=(X-repmat(params_muvec,[N,1]))./(repmat(params_sigmavec,[N,1]));

save('Z:\Documents\MEGMOD\MATLAB\data\WCnetnoisy_mixeddelays_SSI_logdiscrepancies_connstrengths_expMEGdata_75subs_14params_radius1pt5_9062sims.mat','X','Y');