clear;

%% Load (unseen) MEG group-level network

load('Z:\Documents\MEGMOD\MATLAB\data\wPLImats_9.83Hz_30subs_replicates_PARC2k9_148parcs_wgroupmat.mat','wPLImats_strength_groupmat');
MEG_strength_wPLI_FC=wPLImats_strength_groupmat;

%% Estimate SSI values 

N=148;
UTMAT=triu(ones(N),1);
ut_idxs=find(UTMAT(:));

num_reps=5;
num_batches=200;

% Isochronous delays model

SSIvals_isochronousdelaysmodel=zeros(1,num_reps*num_batches);

for batch_idx=1:num_batches
    fname=['Z:\Documents\MEGMOD\MATLAB\data\WCnet_isochronousdelays_noisy_60s_posteriorpredictivechecks_batch',int2str(batch_idx),'.mat'];
    
    if isfile(fname)
    load(fname);
        
        indices=[((batch_idx-1)*num_reps)+1:batch_idx*num_reps];

        for ssi_index=1:length(indices)
            SSIvals_isochronousdelaysmodel(1,indices(ssi_index))=ssim(wPLImags_acrossreps(:,ssi_index),MEG_strength_wPLI_FC(ut_idxs),'Exponents',[1,1,1],'DataFormat',"SS",'DynamicRange',1,'Radius',1.5);
        end

    end

    disp(batch_idx);

end

load('Z:\Documents\MEGMOD\MATLAB\data\WCnetnoisy_isochronousdelays_posteriorpredictivechecks.mat','valid_indices');
SSIvals_isochronousdelaysmodel=SSIvals_isochronousdelaysmodel(valid_indices);

% Mixed delays model

SSIvals_mixeddelaysmodel=zeros(1,num_reps*num_batches);

for batch_idx=1:num_batches
    fname=['Z:\Documents\MEGMOD\MATLAB\data\WCnet_mixeddelays_noisy_60s_posteriorpredictivechecks_batch',int2str(batch_idx),'.mat'];
    
    if isfile(fname)
    load(fname);
        
        indices=[((batch_idx-1)*num_reps)+1:batch_idx*num_reps];

        for ssi_index=1:length(indices)
            SSIvals_mixeddelaysmodel(1,indices(ssi_index))=ssim(wPLImags_acrossreps(:,ssi_index),MEG_strength_wPLI_FC(ut_idxs),'Exponents',[1,1,1],'DataFormat',"SS",'DynamicRange',1,'Radius',1.5);
        end

    end

    disp(batch_idx);

end

load('Z:\Documents\MEGMOD\MATLAB\data\WCnetnoisy_mixeddelays_posteriorpredictivechecks.mat','valid_indices');
SSIvals_mixeddelaysmodel=SSIvals_mixeddelaysmodel(valid_indices);

% Distance-dependent delays model

SSIvals_distancedependentdelaysmodel=zeros(1,num_reps*num_batches);

for batch_idx=1:num_batches
    fname=['Z:\Documents\MEGMOD\MATLAB\data\WCnet_variabledelays_noisy_60s_posteriorpredictivechecks_batch',int2str(batch_idx),'.mat'];
    
    if isfile(fname)
    load(fname);
        
        indices=[((batch_idx-1)*num_reps)+1:batch_idx*num_reps];

        for ssi_index=1:length(indices)
            SSIvals_distancedependentdelaysmodel(1,indices(ssi_index))=ssim(wPLImags_acrossreps(:,ssi_index),MEG_strength_wPLI_FC(ut_idxs),'Exponents',[1,1,1],'DataFormat',"SS",'DynamicRange',1,'Radius',1.5);
        end

    end

    disp(batch_idx);

end

load('Z:\Documents\MEGMOD\MATLAB\data\WCnetnoisy_variabledelays_posteriorpredictivechecks.mat','valid_indices');
SSIvals_distancedependentdelaysmodel=SSIvals_distancedependentdelaysmodel(valid_indices);

%% Thresholding

thresh_vec=[-1:0.02:0];
num_samples=min([length(SSIvals_isochronousdelaysmodel),length(SSIvals_mixeddelaysmodel),length(SSIvals_distancedependentdelaysmodel)]);

SSIlogvals_isochronousdelaysmodel=log(1-SSIvals_isochronousdelaysmodel(1:num_samples));
SSIlogvals_mixeddelaysmodel=log(1-SSIvals_mixeddelaysmodel(1:num_samples));
SSIlogvals_distancedependentdelaysmodel=log(1-SSIvals_distancedependentdelaysmodel(1:num_samples));

figure;
histogram(SSIlogvals_isochronousdelaysmodel,thresh_vec);
axis square;

figure;
histogram(SSIlogvals_mixeddelaysmodel,thresh_vec);
axis square;

figure;
histogram(SSIlogvals_distancedependentdelaysmodel,thresh_vec);
axis square;

%% Estimating model probabilities

num_subthreshidxs=zeros(3,length(thresh_vec));
for thresh_idx=1:length(thresh_vec)
    num_subthreshidxs(1,thresh_idx)=sum(SSIlogvals_isochronousdelaysmodel<thresh_vec(thresh_idx));
    num_subthreshidxs(2,thresh_idx)=sum(SSIlogvals_mixeddelaysmodel<thresh_vec(thresh_idx));
    num_subthreshidxs(3,thresh_idx)=sum(SSIlogvals_distancedependentdelaysmodel<thresh_vec(thresh_idx));
end

model_probabilities=num_subthreshidxs./repmat(sum(num_subthreshidxs),[size(num_subthreshidxs,1),1]);

save('Z:\Documents\MEGMOD\MATLAB\data\model_probabilities.mat','model_probabilities');

%% Line plot

figure;
hold on;
plot(thresh_vec,model_probabilities(1,:),'o','Color',[0,114,178]/255,'LineStyle','-','Linewidth',2);
plot(thresh_vec,model_probabilities(2,:),'o','Color',[213,94,0]/255,'LineStyle','-','Linewidth',2);
plot(thresh_vec,model_probabilities(3,:),'o','Color',[0,158,115]/255,'LineStyle','-','Linewidth',2);
xlim([-1,0]);
ylim([0,1]);
axis square;
xlabel('Threshold');
ylabel('Model probability');
legend('Isochronous delays model','Mixed delays model','Distance-dependent delays model');
