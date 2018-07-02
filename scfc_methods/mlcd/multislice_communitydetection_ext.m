function [S,Q]=multislice_communitydetection_ext(A,N,gamma,omega_vals)

  T=length(A);
  B=spalloc(N*T,N*T,N*N*T+4*N*T);
  twomu=0;
  for s=1:T
      k=sum(A{s});
      twom=sum(k);
      twomu=twomu+twom;
      indx=[1:N]+(s-1)*N;
      B(indx,indx)=A{s}-gamma*k'*k/twom;
   end
    
   twomu=twomu+2*mean([omega_vals{1},omega_vals{2},omega_vals{3}])*N*(T-1);
   
   % inserting first omega value
   
   for idx=0*N+1:N,
   B(idx,idx+N)=omega_vals{1};
   B(idx+N,idx)=omega_vals{1};
   end
      
   % inserting second omega value
   
   for idx=2*N+1:3*N,
   B(idx,idx+N)=omega_vals{2};
   B(idx+N,idx)=omega_vals{2};
   end
   
   % inserting third omega value
   
   B = B + omega_vals{3}*spdiags(ones(N*T,2),[-2*N,2*N],N*T,N*T);
   
   [S,Q] = genlouvain(B,10000,0,0);
   Q = Q/twomu;
   S = reshape(S,2*N,(T/2));

end