function WCnet_mixeddelays_noisy_simulations(array_idx,run_idx)
%% Initialisation

addpath(genpath('/scratch/work/willian1/MEGMOD/MATLAB/'));

tic;
generateGlobals_allsubs();
global PY2LV_transvec avg_sc WCnet_noisy_mixeddelays_paramcombinations source2parc fwd_mat inv_mat
toc;

subject_sets = {'S0001','S0003','S0005','S0006','S0007','S0008','S0011','S0026','S0034','S0035','S0036','S0038',...
    'S0039','S0049','S0113','S0115','S0116','S0117','S0118','S0119','S0120','S0121','S0122','S0123','S0124','S0125','S0126',...
    'S0127','S0128','S0129','S0130','S0132','S0136','S0137','S0140','S0141','S0142','S0143','S0144','S0145','S0146','S0147',...
    'S0148','S0149','S0150','S0151','S0152','S0153','S0155','S0156','S0157','S0161','S0162','S0166','S0167','S0168','S0169',...
    'S0170','S0171','S0172','S0173','S0174','S0175','S0177','S0179','S0180','S0181','S0184','S0185','S0186','S0188','S0189',...
    'S0193','S0194','S0195'};

m=5;
fc=9.83;
fs=1000;
df=4;

N=148;
fs_ds=(fs/df);

lefthem_idxs=[1:2:N];
righthem_idxs=[2:2:N];

UTMAT=triu(ones(N),1);
ut_idxs=find(UTMAT(:));

% SC network

SC=avg_sc;
SC=SC+SC';
SC=SC([lefthem_idxs,righthem_idxs],[lefthem_idxs,righthem_idxs]);

% Distance matrix

load('/scratch/work/willian1/MEGMOD/MATLAB/data/mnicoords_destrieux.mat','mnicoords_destrieux');

DISTMAT=zeros(N);
for ridx=1:N  
    for cidx=1:N  
       DISTMAT(ridx,cidx)=sqrt(sum((mnicoords_destrieux(ridx,:)-mnicoords_destrieux(cidx,:)).^2));     
    end   
end

% Parameter combinations to simulate

num_reps=5;
mod_params=WCnet_noisy_mixeddelays_paramcombinations(((array_idx-1)*num_reps)+1:(array_idx*num_reps),:);

%% Estimate test statistics for model-generated data

ts=(1/fs_ds)*10^3;

mask_length=5;
mask_samples=mask_length*fs_ds;
data_length=65;
data_samples=data_length*(fs_ds*ts);

rng(run_idx);

% Simulate Wilson-Cowan noisy network model

wPLImags_acrossreps=zeros(length(ut_idxs),num_reps);
PLVmags_acrossreps=zeros(length(ut_idxs),num_reps);
delays_matrix=zeros(length(ut_idxs),num_reps);

for rep_idx=1:num_reps

        tic;

        Vij_1=zeros(N);
        Vij_1(ut_idxs)=mod_params(rep_idx,13);
        Vij_1=Vij_1+Vij_1';
        Vij_2=(DISTMAT./mod_params(rep_idx,12));

        Vij=(mod_params(rep_idx,14)*Vij_1)+((1-mod_params(rep_idx,14))*Vij_2);

        Tij=ceil(DISTMAT./Vij);
        Tij(1:N+1:end)=0;

        delays_matrix(:,rep_idx)=Tij(ut_idxs);

        scaling_mat=ones(N);
        scaling_mat([1:N/2],[(N/2)+1:N])=mod_params(rep_idx,10);
        scaling_mat([(N/2)+1:N],[1:N/2])=mod_params(rep_idx,10);
        SC_scaled=scaling_mat.*SC;
        
        SC_wt=double(SC_scaled);
        Kij = SC_wt./repmat(sum(SC_wt,2),[1,N]);

        Je = 0;
        Ji = 0;
        sys = WilsonCowanNet_variabledelays_noisy(Kij,Tij,Je,Ji);
        sys.pardef=bdSetValue(sys.pardef,'wee',mod_params(rep_idx,1));
        sys.pardef=bdSetValue(sys.pardef,'wei',mod_params(rep_idx,2));
        sys.pardef=bdSetValue(sys.pardef,'wie',mod_params(rep_idx,3));
        sys.pardef=bdSetValue(sys.pardef,'wii',mod_params(rep_idx,4));
        sys.pardef=bdSetValue(sys.pardef,'k',mod_params(rep_idx,9));
        sys.pardef=bdSetValue(sys.pardef,'be',mod_params(rep_idx,5));
        sys.pardef=bdSetValue(sys.pardef,'bi',mod_params(rep_idx,6));
        sys.pardef=bdSetValue(sys.pardef,'taue',mod_params(rep_idx,7));
        sys.pardef=bdSetValue(sys.pardef,'taui',mod_params(rep_idx,8));
        sys.pardef=bdSetValue(sys.pardef,'psi_sigma',mod_params(rep_idx,11));
        
        % Simulate model
    
        tspan = [0 data_samples];                     % Integration time domain
        sol = bdSolve(sys,tspan,@dde23a,'ddesolver'); % Apply the dde23a solver
        tplot = ts:ts:data_samples;                   % Interpolation time points
        Y = bdEval(sol,tplot);                        % Interpolate the solution
    
        D=Y(1:size(Y,1)/2,:);
        
        wPLImags_mat=zeros(N,N,length(subject_sets));
        PLVmags_mat=zeros(N,N,length(subject_sets));
    
        for sub_idx=1:length(subject_sets)
                              
            % Linear mixing
            D_linmix = forwardinverse_modelling_bothems(D,source2parc{sub_idx},PY2LV_transvec,[lefthem_idxs,righthem_idxs],fwd_mat{sub_idx},inv_mat{sub_idx});
            
            % Morlet filtering
            D_linmix_filt = wavelet_filtering(D_linmix,m,fc,fs_ds);           
            
            % Estimate wPLI
            [~,wPLImags]=wPLI_calc(D_linmix_filt(:,mask_samples+1:end));           
            wPLImags(1:N+1:end)=0;
            wPLImags_mat(:,:,sub_idx)=wPLImags+wPLImags';

            % Estimate Phase Locking Value
            [PLVmags,~,~,~]= PLVcalc(D_linmix_filt(:,mask_samples+1:end),'morlet');
            PLVmags(1:N+1:end)=0;
            PLVmags_mat(:,:,sub_idx)=PLVmags+PLVmags';
                        
        end
          
        wPLImags_groupmat=mean(wPLImags_mat,3);
        PLVmags_groupmat=mean(PLVmags_mat,3);
        wPLImags_acrossreps(:,rep_idx)=wPLImags_groupmat(ut_idxs);
        PLVmags_acrossreps(:,rep_idx)=PLVmags_groupmat(ut_idxs);
        
      toc; 
      disp(rep_idx);
        
    save(['/scratch/work/willian1/MEGMOD/MATLAB/data/WCnet_mixeddelays_noisy_60s_simulations_batch',int2str(array_idx),'_run',int2str(run_idx),'.mat'],...
    'wPLImags_acrossreps','PLVmags_acrossreps','mod_params','delays_matrix');

end

end