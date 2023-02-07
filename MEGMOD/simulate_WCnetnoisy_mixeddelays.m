function out = simulate_WCnetnoisy_mixeddelays(wee,wei,wie,wii,be,bi,taue,taui,k,IH,psi_sigma,delay,velocity,balance_coeff)

addpath(genpath('Z:\\Documents\MEGMOD\MATLAB\code\bdtoolkit-2019a\'));

% Subject IDs

subject_sets = {'S0001','S0003','S0005','S0006','S0007','S0008','S0011','S0026','S0034','S0035','S0036','S0038',...
    'S0039','S0049','S0113','S0115','S0116','S0117','S0118','S0119','S0120','S0121','S0122','S0123','S0124','S0125','S0126',...
    'S0127','S0128','S0129','S0130','S0132','S0136','S0137','S0140','S0141','S0142','S0143','S0144','S0145','S0146','S0147',...
    'S0148','S0149','S0150','S0151','S0152','S0153','S0155','S0156','S0157','S0161','S0162','S0166','S0167','S0168','S0169',...
    'S0170','S0171','S0172','S0173','S0174','S0175','S0177','S0179','S0180','S0181','S0184','S0185','S0186','S0188','S0189',...
    'S0193','S0194','S0195'};

% Specify variables

N=148;
lefthem_idxs=[1:2:N];
righthem_idxs=[2:2:N];
    
% Loading required data

load('Z:\\Documents\MEGMOD\MATLAB\data\PY2LV_transvec.mat','PY2LV_transvec');
load('Z:\\Documents\MEGMOD\MATLAB\data\SC_parc2k9_57HCP.mat','avg_sc');

% Loading distance matrix

load('Z:\\Documents\MEGMOD\MATLAB\data\mnicoords_destrieux.mat','mnicoords_destrieux');

DISTMAT=zeros(N);
for ridx=1:N  
    for cidx=1:N  
       DISTMAT(ridx,cidx)=sqrt(sum((mnicoords_destrieux(ridx,:)-mnicoords_destrieux(cidx,:)).^2));     
    end   
end

% SC network

SC=avg_sc;
SC=SC+SC';
SC=SC([lefthem_idxs,righthem_idxs],[lefthem_idxs,righthem_idxs]);

UTMAT=triu(ones(N),1);
ut_idxs=logical(UTMAT(:));

% Parameter values

m=5;
fc=9.83;
fs=250;
ts=(1/fs)*10^3;

mask_length=5;
mask_samples=mask_length*fs;
data_length=65;
data_samples=data_length*(fs*ts);

% Scaling parameters

wee_scaled=5*wee+20;
wei_scaled=6*wei+18;
wie_scaled=6*wie+18;
wii_scaled=0.2*wii+1;
be_scaled=1*be+3;
bi_scaled=1*bi+5;
taue_scaled=3.8*taue+18.6;
taui_scaled=4.7*taui+15.1;
k_scaled=0.5*k+1.5;
IH_scaled=0.5*IH+2.5;
psi_sigma_scaled=0.05*psi_sigma+0.15;
delay_scaled=3*delay+10;
velocity_scaled=2*velocity+8;
balance_coeff_scaled=0.15*balance_coeff+0.5;

% Parameters for Wilson-Cowan network model

Vij_1=zeros(N);
Vij_1(ut_idxs)=velocity_scaled;
Vij_1=Vij_1+Vij_1';
Vij_2=(DISTMAT./delay_scaled);

Vij=(balance_coeff_scaled*Vij_1)+((1-balance_coeff_scaled)*Vij_2);

Tij=ceil(DISTMAT./Vij);
Tij(1:N+1:end)=0;
      
scaling_mat=ones(N);
scaling_mat([1:N/2],[(N/2)+1:N])=IH_scaled;
scaling_mat([(N/2)+1:N],[1:N/2])=IH_scaled;
SC_scaled=scaling_mat.*SC;
        
SC_wt=double(SC_scaled);
Kij = SC_wt./(repmat(sum(SC_wt,2),[1,N]));

Je = 0;
Ji = 0;
sys = WilsonCowanNet_variabledelays_noisy(Kij,Tij,Je,Ji);
sys.pardef=bdSetValue(sys.pardef,'wee',wee_scaled);
sys.pardef=bdSetValue(sys.pardef,'wei',wei_scaled);
sys.pardef=bdSetValue(sys.pardef,'wie',wie_scaled);
sys.pardef=bdSetValue(sys.pardef,'wii',wii_scaled);
sys.pardef=bdSetValue(sys.pardef,'k',k_scaled);
sys.pardef=bdSetValue(sys.pardef,'be',be_scaled);
sys.pardef=bdSetValue(sys.pardef,'bi',bi_scaled);
sys.pardef=bdSetValue(sys.pardef,'taue',taue_scaled);
sys.pardef=bdSetValue(sys.pardef,'taui',taui_scaled);    
sys.pardef=bdSetValue(sys.pardef,'psi_sigma',psi_sigma_scaled);

% Simulate model
 
tspan = [0 data_samples];                              % Integration time domain
sol = bdSolve(sys,tspan,@dde23a,'ddesolver');          % Apply the dde23a solver
tplot = ts:ts:data_samples;                            % Interpolation time points
Y = bdEval(sol,tplot);                                 % Interpolate the solution

D=Y(1:size(Y,1)/2,:);

wPLImags_mat=zeros(N,N,length(subject_sets));

parfor sub_idx=1:length(subject_sets)

    % Source to parcels
    fname = char(strcat('Z:\\Documents\MEGMOD\MATLAB\data\MEG_subjectwise_data\sourceIDs\',subject_sets(sub_idx),'_set01__parc2009.csv'));
    source2parc=table2array(readtable(fname,'ReadVariableNames',0,'Format','auto'));           
            
    % Forward operator
    fname = char(strcat('Z:\\Documents\MEGMOD\MATLAB\data\MEG_subjectwise_data\fixed_fwd_ops\',subject_sets(sub_idx),'_set01__py-fwd.csv'));
    fwd_mat=table2array(readtable(fname,'ReadVariableNames',0,'Format','auto'));
                        
    % Inverse operator
    fname = char(strcat('Z:\\Documents\MEGMOD\MATLAB\data\MEG_subjectwise_data\fidweighted_inv_ops\',subject_sets(sub_idx),'_set01__parc2009.csv'));
    inv_mat=table2array(readtable(fname,'ReadVariableNames',0,'Format','auto'));
            
    % Linear mixing
    D_linmix = forwardinverse_modelling_bothems(D,source2parc,PY2LV_transvec,[lefthem_idxs,righthem_idxs],fwd_mat,inv_mat);
             
    % Morlet filtering      
    D_linmix_filt = wavelet_filtering(D_linmix,m,fc,fs);           
            
    % Estimate weighted Phase Lag Index (wPLI)
    [~,wPLImags]=wPLI_calc(D_linmix_filt(:,mask_samples+1:end));                    
    wPLImags(1:N+1:end)=0;
    wPLImags_mat(:,:,sub_idx)=wPLImags+wPLImags';
                           
end
     
wPLImags_groupmat=mean(wPLImags_mat,3);
out=wPLImags_groupmat(ut_idxs);

end