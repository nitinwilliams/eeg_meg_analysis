function LL=markovmodel_loglhood(SMAT,TP)
% finds log-likelihood value of multi-trial sequence of states, for a given
% matrix of transition probabilities
%
% Inputs:
% SMAT - MxN matrix of multi-trial sequence of states (M trials, N states per trial)
% TP - matrix of transition probabilities according to which log-likelihood
% value should estimated
% Outputs:
% LL - log-likelihood value
%
% NOTE 1: assumes, for e.g. states 1 to Q are expressed in STMAT as 1,2,3...Q
%
% NOTE 2: paths with transition probabilities of 0 are set to infinitesimally
% small value (i.e. eps=2.2204e-16)
%

[M,N]=size(SMAT);
LL=zeros(M,1);

% setting zero values to 'eps'

TP(TP(:)==0)=eps;

% Computing log-likelihood

for ridx=1:M
    
    probvec=zeros(1,N-1);
    
   for cidx=1:N-1
    
       rnum=SMAT(ridx,cidx);
       cnum=SMAT(ridx,cidx+1);
       
       probvec(1,cidx)=TP(rnum,cnum);
       
   end
        
   LL(ridx)=sum(log(probvec));
   
end


% LL(LL==-Inf)=NaN;
LL=nansum(LL);

end