function [S,Q,zQ]=multislice_communitydetection(A,N,gamma,omega,reps,mode)

  T=length(A);
  B=spalloc(N*T,N*T,N*N*T+2*N*T);
  twomu=0;
  
  for s=1:T
      k=sum(A{s});
      twom=sum(k);
      twomu=twomu+twom;
      indx=[1:N]+(s-1)*N;
      B(indx,indx)=A{s}-gamma*(k'*k)/twom;
  end
    
  twomu=twomu+2*omega*N*(T-1);
  B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
  [S,Q] = genlouvain(B,10000,0,0);
  Q = Q/twomu;
  S = reshape(S,N,T);
  
  zQ=NaN;
  
  if strcmp(mode,'nodalnullmodel')
  
  % null distribution of Q (nodal null model)
  
  Q_null=zeros(reps,1);
  orig_ord=zeros(N,2);
  
  orig_ord(1:N,1)=[1:N]';
  orig_ord(1:N,2)=(N+[1:N])';
  orig_idxs=sub2ind([2*N,2*N],orig_ord(:,1),orig_ord(:,2));
  
  parfor repidx=1:reps,
      
      B_null= B - omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
      perm_ord=orig_ord;
      perm_ord(:,2)=orig_ord(randperm(N),2);
      perm_idxs_ut=sub2ind([2*N,2*N],perm_ord(:,1),perm_ord(:,2));
      perm_idxs_lt=sub2ind([2*N,2*N],perm_ord(:,2),perm_ord(:,1));
      B_null(perm_idxs_ut)=B(orig_idxs);
      B_null(perm_idxs_lt)=B(orig_idxs);
      [~,Q_null(repidx,1)] = genlouvain(B_null,10000,0,0);
      Q_null(repidx,1) = Q_null(repidx,1)/twomu;
      
  end

  zQ=(Q-mean(Q_null))/std(Q_null);
  
  end
  
end