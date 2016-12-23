% below is code to run pipeline on two subjects from BCI oddball task.
%
% please adapt code for your purposes after adding pipeline
% folders/functions to current working directory

clear;
addpath(genpath(pwd));

%% MATLAB pool (start parallel pool)

P=cbupool(16);
P.ResourceTemplate = '-l nodes=^N^,mem=128GB,walltime=150:00:00';
parpool(P);

% parpool close;
% parpool open;

%% Load data

bci_strings={'bcidata_6C_3D','bcidata2_6C_3D'};
k=[2:8];        % samples for baseline window
numcycs=1000;   % samples per trial

CA_nd=repmat({NaN(length(bci_strings),length(k))},[numcycs,1]);
MD_nd=repmat({NaN(length(bci_strings),length(k))},[numcycs,1]);

for bci_idx=1:length(bci_strings),

load(bci_strings{bci_idx},'bcidata_3D'); 

load('chans','chans');
chan_nums=[]; for chidx=1:length(chans), chan_nums(chidx,1)=chans{chidx,2}; end
bcidata_3D=bcidata_3D(chan_nums,:,:); % bcidata_3D is 3D matrix with channels x samples x trials of both conditions

%% Baseline correction

bpts=42;  % samples for baseline window
tpts=240; % samples per trial

bcidata_3D_bc=zeros(size(bcidata_3D,1),size(bcidata_3D,2)-bpts,size(bcidata_3D,3)); 

for idx1=1:size(bcidata_3D,3),
    
    for idx2=1:size(bcidata_3D,1),    
    
         bcidata_3D_bc(idx2,:,idx1)=bcidata_3D(idx2,bpts+[1:tpts],idx1)-mean(bcidata_3D(idx2,1:bpts,idx1));  
        
    end
    
end

%% Pre-processing 1 (z-scoring each trial)

bcidata_3D_bc_z1=zeros(size(bcidata_3D_bc));

for idx1=1:size(bcidata_3D_bc,3),
    
   bcidata_3D_bc_z1(:,:,idx1)=zscore(squeeze(bcidata_3D_bc(:,:,idx1))')';  
    
end

%% Pre-processing 2 (z-scoring across trials)

bcidata_3D_bc_z1_z2=zeros(size(bcidata_3D_bc_z1));

r1=[1:floor(size(bcidata_3D_bc_z1,3)/2)];
r2=[floor(size(bcidata_3D_bc_z1,3)/2)+1:size(bcidata_3D_bc_z1,3)];

% Condition 1

for idx1=1:size(bcidata_3D_bc_z1,1),
    
    bcidata_3D_bc_z1_z2(idx1,:,r1)=zscore(squeeze(bcidata_3D_bc_z1(idx1,:,r1))')';
    
end

% Condition 2

for idx1=1:size(bcidata_3D_bc_z1,1),
    
   bcidata_3D_bc_z1_z2(idx1,:,r2)=zscore(squeeze(bcidata_3D_bc_z1(idx1,:,r2))')';
      
end

%% converting to cell array format

Ys=cell(1,2);

% Condition 1

for trialidx=r1,
       
      Ys{1,1}.trial{1,trialidx}=squeeze(bcidata_3D_bc_z1_z2(:,:,trialidx)); 
       
end

% Condition 2

for trialidx=r2,
       
      Ys{1,2}.trial{1,trialidx-length(r1)}=squeeze(bcidata_3D_bc_z1_z2(:,:,trialidx)); 
       
end

%% Segmentation

srate=240;             % in Hz
M=28;                  % number of electrodes
tpseg=10;              % no. of trials per estimation window
seglen=500;            % width of 1 trial x no. of trials per estimation window (in milliseconds)
T=round(seglen/(1000/srate));
N=round(T/tpseg);

trial_len=size(Ys{1,1}.trial{1,1},2);
segspt=floor(trial_len/N);
segspcond=(size(Ys{1,1}.trial,2)/tpseg)*segspt;
numsegs=segspcond*size(Ys,2);

SEGMAT=cell(numsegs,1);

for Condidx=1:size(Ys,2),
    
    for Tidx=1:tpseg:size(Ys{1,1}.trial,2),
               
        for segidx=1:segspt,
            
           MAR_D=[];
            
           for tsegidx=1:tpseg,
              
            D=Ys{1,Condidx}.trial{1,(Tidx-1)+tsegidx}(:,((segidx-1)*N)+1:segidx*N)';
            MAR_D=vertcat(MAR_D,D);   
         
           end
           
           segid=segidx+floor((Tidx-1)/tpseg)*segspt+(Condidx-1)*segspcond;
           SEGMAT{segid,1}=MAR_D;  
            
        end
        
    end
    
end

%% finding lambda (sMAR)

MOPrange=[1:1:6];
pctval=10;
trunc_len=floor(size(SEGMAT,1)*(pctval/100));
trunc_indices=randperm(size(SEGMAT,1));
trunc_indices=trunc_indices(1,[1:trunc_len])';

lambdaGCV=repmat({NaN(1,length(MOPrange))},[trunc_len 1]);

parfor segidx=1:trunc_len,
    
    DAT=SEGMAT{trunc_indices(segidx,1),1};
    [lambdaGCV{segidx}(1,:)]=findlambda(DAT,MOPrange);
    
    display(segidx);
          
end

lambdaGCV=cell2mat(lambdaGCV);
lambdaGCVvec=nanmean(lambdaGCV);

%% finding model order (sMAR)

pctval_2=10;
trunc_len_2=floor(size(SEGMAT,1)*(pctval_2/100));
trunc_indices_2=randperm(size(SEGMAT,1));
trunc_indices_2=trunc_indices_2(1,[1:trunc_len_2])';

GCV=repmat({NaN(1,length(MOPrange))},[trunc_len_2 1]);
AIC1=repmat({NaN(1,length(MOPrange))},[trunc_len_2 1]);
AIC2=repmat({NaN(1,length(MOPrange))},[trunc_len_2 1]);
BIC1=repmat({NaN(1,length(MOPrange))},[trunc_len_2 1]);
BIC2=repmat({NaN(1,length(MOPrange))},[trunc_len_2 1]);

parfor segidx=1:trunc_len_2,
    
DAT=SEGMAT{trunc_indices_2(segidx,1),1};
[GCV{segidx}(1,:),AIC1{segidx}(1,:),AIC2{segidx}(1,:),BIC1{segidx}(1,:),BIC2{segidx}(1,:)]=findmodelorder(DAT,MOPrange,lambdaGCVvec);

display(segidx);

end

GCV=cell2mat(GCV);
AIC1=cell2mat(AIC1);
BIC1=cell2mat(BIC1);
AIC2=cell2mat(AIC2);
BIC2=cell2mat(BIC2);

[c,ind]=min(real(BIC1),[],2);
MOP=mode(ind);

%% sMAR modelling

lambda(bci_idx)=lambdaGCVvec(1,MOP);
Beta=zeros(numsegs,M*MOP,M);
Bnorm=zeros(numsegs,M,M);

for segidx=1:size(SEGMAT,1),
    
    DAT=SEGMAT{segidx,1};

    [Beta(segidx,:,:)]=mpls(DAT,MOP,lambda(bci_idx));

    Bnorm(segidx,:,:)=mean(sqrt(reshape(squeeze(Beta(segidx,:,:))',[M M MOP]).^2),3);
    
    Bnorm(segidx,:,:)=squeeze(Bnorm(segidx,:,:))-diag(diag(squeeze(Bnorm(segidx,:,:))));
    
    display(segidx);
    
end

% Clustering

CMAT=zeros(size(SEGMAT,1),M^2);

for segidx=1:size(SEGMAT,1),
    
CMAT(segidx,:)=reshape(squeeze(Bnorm(segidx,:,:)),[1 M^2]);   
    
end

CMAT=CMAT(:,sum(CMAT)~=0);

for kidx=1:length(k),

% k-means clustering

[symbols,C,SUMD,D]=kmeans(CMAT,k(kidx),'Distance','correlation','Start','sample','Replicates',1000);

% assembling symbol matrices

% assembling into matrices/condition

SYMBOLMAT=reshape(symbols,[segspt floor(size(Ys{1,1}.trial,2)/tpseg) size(Ys,2)]); 

TOTMAT=cell(1,1);

for Condidx=1:size(Ys,2),
    
  CONDMAT{1,Condidx}=squeeze(SYMBOLMAT(:,:,Condidx))';
  TOTMAT{1,1}=vertcat(TOTMAT{1,1},CONDMAT{1,Condidx});
    
end

Q=k(kidx);

% HMM estimation

% HMM training

% condition 1

prior1=normalise(rand(Q,1));
transmat1=mk_stochastic(rand(Q,Q));
obsmat1=eye(Q); 

[LL1,PR1,TR1,EMIT1]=dhmm_em(CONDMAT{1,1},prior1,transmat1,obsmat1,'max_iter',100,'verbose',0);

% condition 2

prior2=normalise(rand(Q,1));
transmat2=mk_stochastic(rand(Q,Q));
obsmat2=eye(Q); 

[LL2,PR2,TR2,EMIT2]=dhmm_em(CONDMAT{1,2},prior2,transmat2,obsmat2,'max_iter',100,'verbose',0);

% HMM classification

numtrials=size(CONDMAT{1,1},1);         % assumes equal number of trials for both conditions
numfolds=2;
numtperfold=round(numtrials/numfolds);
foldvec=[0:numtperfold:numtrials];
foldmat=[];
for fidx=1:numfolds,
    foldmat(fidx,:)=[foldvec(fidx)+1:foldvec(fidx+1)];
end

% randomising order of trials

CONDMAT{1,1}=CONDMAT{1,1}(randperm(numtrials),:);
CONDMAT{1,2}=CONDMAT{1,2}(randperm(numtrials),:);

TESTMAT=zeros(numtrials,2);

for foldidx=1:numfolds,
    
    FMAT=circshift(foldmat,foldidx-1);
        
    % finding training indices
    trainidxs=FMAT(1:end-1,:);
    testidxs=FMAT(end,:);
        
    % getting training datasets
    trainset1=CONDMAT{1,1}(trainidxs(:),:);
    trainset2=CONDMAT{1,2}(trainidxs(:),:);
    
    % setting priors
    prior=normalise(rand(Q,1));
    transmat=mk_stochastic(rand(Q,Q));
    obsmat=eye(Q); 
    
    % training
    [LL,PRtr1,TRtr1,EMITtr1]=dhmm_em(trainset1,prior,transmat,obsmat,'max_iter',100,'verbose',0);
    [LL,PRtr2,TRtr2,EMITtr2]=dhmm_em(trainset2,prior,transmat,obsmat,'max_iter',100,'verbose',0);
    
    for trialidx=1:size(FMAT,2),    
    
    % testing
    loglik1=dhmm_logprob(CONDMAT{1,1}(testidxs(trialidx),:),PRtr1,TRtr1,EMITtr1);
    loglik2=dhmm_logprob(CONDMAT{1,1}(testidxs(trialidx),:),PRtr2,TRtr2,EMITtr2);
    TESTMAT(testidxs(trialidx),1)=(loglik1>loglik2);    
    
    loglik1=dhmm_logprob(CONDMAT{1,2}(testidxs(trialidx),:),PRtr1,TRtr1,EMITtr1);
    loglik2=dhmm_logprob(CONDMAT{1,2}(testidxs(trialidx),:),PRtr2,TRtr2,EMITtr2);
    TESTMAT(testidxs(trialidx),2)=(loglik2>loglik1);   
    
    end
   
end

CA(bci_idx,kidx)=sum(TESTMAT(:));

%% Permutation testing

% CA

numtrials=size(CONDMAT{1,1},1);   % assumes equal number of trials for both conditions
DMAT=[CONDMAT{1,1};CONDMAT{1,2}];
hp=floor(size(DMAT,1)/2);

parfor cycidx=1:numcycs,
       
     DMAT_nd=DMAT(randperm(size(DMAT,1)),:);
     D1=DMAT_nd(1:hp,:);
     D2=DMAT_nd(hp+1:2*hp,:);
        
     TESTMAT=zeros(numtrials,2);

     for foldidx=1:numfolds,
    
     FMAT=circshift(foldmat,foldidx-1);
        
     % finding training indices
     trainidxs=FMAT(1:end-1,:);
     testidxs=FMAT(end,:);
        
     % getting training datasets
     trainset1=D1(trainidxs(:),:);
     trainset2=D2(trainidxs(:),:);
    
     % setting priors
     prior=normalise(rand(Q,1));
     transmat=mk_stochastic(rand(Q,Q));
     obsmat=eye(Q); 
    
     % training
     [LL,PRtr1,TRtr1,EMITtr1]=dhmm_em(trainset1,prior,transmat,obsmat,'max_iter',100,'verbose',0);
     [LL,PRtr2,TRtr2,EMITtr2]=dhmm_em(trainset2,prior,transmat,obsmat,'max_iter',100,'verbose',0);
     
     for trialidx=1:size(FMAT,2),    
    
        % testing
        loglik1=dhmm_logprob(D1(testidxs(trialidx),:),PRtr1,TRtr1,EMITtr1);
        loglik2=dhmm_logprob(D1(testidxs(trialidx),:),PRtr2,TRtr2,EMITtr2);
        TESTMAT(testidxs(trialidx),1)=(loglik1>loglik2);    
    
        loglik1=dhmm_logprob(D2(testidxs(trialidx),:),PRtr1,TRtr1,EMITtr1);
        loglik2=dhmm_logprob(D2(testidxs(trialidx),:),PRtr2,TRtr2,EMITtr2);
        TESTMAT(testidxs(trialidx),2)=(loglik2>loglik1);   
    
     end
   
    end
    
    CA_nd{cycidx,1}(bci_idx,kidx)=sum(TESTMAT(:));
    
    display(cycidx);
    
end

%% MM distance

Tval=size(CONDMAT{1,1},2);
LL1_D12=dhmm_logprob(CONDMAT{1,2},PR1,TR1,EMIT1);
LL2_D12=dhmm_logprob(CONDMAT{1,2},PR2,TR2,EMIT2);

LL1_D21=dhmm_logprob(CONDMAT{1,1},PR2,TR2,EMIT2);
LL2_D21=dhmm_logprob(CONDMAT{1,1},PR1,TR1,EMIT1);

D12=(LL1_D12-LL2_D12)/Tval;
D21=(LL1_D21-LL2_D21)/Tval;

MD(bci_idx,kidx)=mean([D12,D21]);

% MM-distance null distribution

DMAT=TOTMAT{1,1};
hp=floor(size(DMAT,1)/2);

parfor cycidx=1:numcycs,
       
    DMAT_nd=DMAT(randperm(size(DMAT,1)),:);
    D1=DMAT_nd(1:hp,:);
    D2=DMAT_nd(hp+1:2*hp,:);
    
    % null condition 1

    prior1_nd=normalise(rand(Q,1));
    transmat1_nd=mk_stochastic(rand(Q,Q));
    obsmat1_nd=eye(Q); 

    [LL1_nd,PR1_nd,TR1_nd,EMIT1_nd]=dhmm_em(D1,prior1_nd,transmat1_nd,obsmat1_nd,'max_iter',100,'verbose',0);

    % null condition 2

    prior2_nd=normalise(rand(Q,1));
    transmat2_nd=mk_stochastic(rand(Q,Q));
    obsmat2_nd=eye(Q); 

    [LL2_nd,PR2_nd,TR2_nd,EMIT2_nd]=dhmm_em(D2,prior2_nd,transmat2_nd,obsmat2_nd,'max_iter',100,'verbose',0);
    
    % generating null distribution
    
    LL1_D12=dhmm_logprob(D2,PR1_nd,TR1_nd,EMIT1_nd);
    LL2_D12=dhmm_logprob(D2,PR2_nd,TR2_nd,EMIT2_nd);

    LL1_D21=dhmm_logprob(D1,PR2_nd,TR2_nd,EMIT2_nd);
    LL2_D21=dhmm_logprob(D1,PR1_nd,TR1_nd,EMIT1_nd);

    D12=(LL1_D12-LL2_D12)/Tval;
    D21=(LL1_D21-LL2_D21)/Tval;

    MD_nd{cycidx,1}(bci_idx,kidx)=mean([D12,D21]);
   
end

display(bci_idx);
display(kidx);

end

end

%% obtaining p-values for CA and MD (for up to 2 subjects)

CA_nd=cell2mat(CA_nd); % contains null distribuion for classification accuracy
MD_nd=cell2mat(MD_nd); % contains null distribution for model distance

CA_ndm(:,:,1)=(CA_nd([1:2:size(CA_nd,1)],:))';
CA_ndm(:,:,2)=(CA_nd([2:2:size(CA_nd,1)],:))';

CA_ndm=shiftdim(CA_ndm,2);

MD_ndm(:,:,1)=(MD_nd([1:2:size(MD_nd,1)],:))';
MD_ndm(:,:,2)=(MD_nd([2:2:size(MD_nd,1)],:))';

MD_ndm=shiftdim(MD_ndm,2);

% p-values for CA

pval_permCA=zeros(size(CA));
[M,N]=size(pval_permCA);

for ridx=1:M,
    
   for cidx=1:N, 
           
      t=100;
      while CA(ridx,cidx)<=prctile(CA_ndm(ridx,cidx,:),t),
          t=t-0.1;
      end
          
      pval_permCA(ridx,cidx)=1-(t/100);
      
   end
       
end

% p-values for MD

pval_MD=zeros(size(MD));
[M,N]=size(pval_MD);

for ridx=1:M,
    
   for cidx=1:N, 
    
      try 
       
      t=0;
      while MD(ridx,cidx)>=prctile(MD_ndm(ridx,cidx,:),t),
          t=t+0.1;
      end
      
      catch ME
          
         t=100;
          
      end
    
      pval_MD(ridx,cidx)=(t/100);
      
   end
       
end
