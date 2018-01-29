function [lambdaGCV]=findlambda(DAT,MOPrange)

[T,M]=size(DAT);

for ordidx=1:length(MOPrange),  
 
q = M;    
MOP=MOPrange(1,ordidx);
p=MOP*q;
e = ones(p, T);
L = spdiags([e -e],[0:1],p,p);

PInfo.stype_penalty = 'l1';
PInfo.L = L;

Opts.IC=0;
Opts.verbose=0;
Opts.nlambdas=500;
Opts.hmethod=@MM;

Y=DAT(MOP+1:end,:);

X=zeros(T-MOP,p);

for ridx=1:MOP,
X(:,[1:q]+(ridx-1)*q)=DAT(MOP-(ridx-1):end-ridx,:);
end

[Betas, lambdas, OUTPUT] = penalty4lambdas(Y, X, [], PInfo, Opts);

df=OUTPUT.df;
RSS=OUTPUT.RSS;
sigma2=OUTPUT.sigma2;

Optsgcv.f_graphic = 0;
[beta_GCV,lsel_GCV,critsel_GCV,critvals_GCV,ind_GCV] = lambda4gcv(Betas, lambdas, T-MOP, RSS, df, Optsgcv);

lambdaGCV(1,ordidx)=lsel_GCV;

end

end