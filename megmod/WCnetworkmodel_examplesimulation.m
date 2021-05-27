clear;
addpath(genpath('D:\DH7-H916-03\HIIT\bdtoolkit-2019a')); % replace with path to Brain Dynamics Toolbox

%% Initialise 

N=74;       % number of brain regions
pctval=90;  % percentile threshold for generating binary SC

% Generating random structural connectome (SC) to link oscillators - replace with estimated SC

SC_rand=rand(N);
SC=triu(SC_rand,1);
SC=SC+SC';

UTMAT=triu(ones(N),1);
ut_idxs=find(UTMAT(:));

pct_thresh=prctile(SC(ut_idxs),pctval);
SC_bin=double(SC>=pct_thresh);                          % Binary SC of just top 10 percentile connections
Kij = SC_bin./(repmat(sum(SC_bin,2),[1,N]).*ones(N));   % Row-sum normalised version of binary SC

% Parameter values

m=5;            % width parameter of Morlet wavelet
fc=9.83;        % centre frequency for Morlet filtering (Hz)
fs=500;         % sampling frequency of simulated data (Hz)
ts=(1/fs)*10^3; % time interval between each sample of simulated data (ms)

mask_length=5;                      % segment length of simulated data to mask out due to potential transients (s)
mask_samples=mask_length*fs;        % number of samples to mask out due to potential transients
data_length=65;                     % total length of simulated data (s)
data_samples=data_length*(fs*ts);   % number of samples of simulated data

%% Model simulation (requires Brain Dyanmics Toolbox (https://github.com/breakspear/bdtoolkit) on MATLAB path)
     
% Specifying model parameters    
        
Je = 2;
Ji = 0;
sys = WilsonCowanNet(Kij,Je,Ji);
sys.pardef=bdSetValue(sys.pardef,'wee',16);
sys.pardef=bdSetValue(sys.pardef,'wei',12);
sys.pardef=bdSetValue(sys.pardef,'wie',15);
sys.pardef=bdSetValue(sys.pardef,'wii',3);
sys.pardef=bdSetValue(sys.pardef,'k',2);
sys.pardef=bdSetValue(sys.pardef,'be',4);
sys.pardef=bdSetValue(sys.pardef,'bi',3.7);
sys.pardef=bdSetValue(sys.pardef,'taue',23.7);
sys.pardef=bdSetValue(sys.pardef,'taui',23.7);    
    
% Running model simulations
 
tspan = [0 data_samples];               % Integration time domain
sol = bdSolve(sys,tspan,@ode45);        % Apply the ode45 solver
tplot = ts:ts:data_samples;             % Interpolation time points
Y = bdEval(sol,tplot);                  % Interpolate the solution

%% Processing data from excitatory populations of Wilson-Cowan (WC) oscillators

D=Y(1:size(Y,1)/2,:);           % extracting data from just excitatory populations of WC oscillators
           
% Applying forward and inverse operator to simulated model data

num_sensors=306;                % number of MEG sensors
fwd_mat=rand(num_sensors,N);    % random forward matrix relating brain regions to MEG sensors (replace with estimated forward matrix)    
inv_mat=pinv(fwd_mat);          % pseudo inverse of random forward matrix (replace with estimated inverse matrix)
D_linmix = inv_mat*(fwd_mat*D);
             
% Band-pass filtering (Morlet filtering) of simulated data from each brain region      
D_linmix_filt = wavelet_filtering(D_linmix,m,fc,fs);           
            
% Estimate phase-synchronization between simulated data from every pair of regions
[~,wPLImags]=wPLI_calc(D_linmix_filt(:,mask_samples+1:end));               
wPLImags(1:N+1:end)=0;
wPLImags=wPLImags+wPLImags';
