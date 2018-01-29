function [GCV,AIC1,AIC2,BIC1,BIC2]=findmodelorder(DAT,MOPrange,lambdaGCVvec)

[T,M]=size(DAT);

for ordidx=1:length(MOPrange), 

q=M;
MOP=MOPrange(1,ordidx);
p=MOP*q;
e = ones(p, T);
L = spdiags([e -e],[0:1],p,p);
PInfo.stype_penalty = 'l1';
PInfo.L = L;

Opts.IC=0;
Opts.verbose=1;
Opts.hmethod=@MM;

Y=DAT(MOP+1:end,:);
X=zeros(T-MOP,p);

for ridx=1:MOP,
X(:,[1:q]+(ridx-1)*q)=DAT(MOP-(ridx-1):end-ridx,:);
end

% applying MM

[Betas, lambdas, OUTPUT] = penalty4lambdas(Y, X, lambdaGCVvec(1,ordidx), PInfo, Opts);

% estimating GCV

RSS=OUTPUT.RSS;
df=OUTPUT.df;
sigma2=OUTPUT.sigma2;

GCV(1,ordidx)=RSS/((T-MOP)-mean(df))^2;
AIC1(1,ordidx)=log(sigma2) + (2-1)*df/(T-MOP);
AIC2(1,ordidx)=log(sigma2) + (2-1)*(df+1)/(T-MOP);
BIC1(1,ordidx)=log(sigma2) + (log(T-MOP)-1)*df/(T-MOP);
BIC2(1,ordidx)=log(sigma2) + (log(T-MOP))*df/(T-MOP);

end

end