function [Beta]=mpls(DAT,MOP,lambda)

% setting parameters

[T,M]=size(DAT);

q = M;
p = MOP*q;
e = ones(p, T);
L = spdiags([e -e],[0:1],p,p);

PInfo.stype_penalty = 'l1';
PInfo.lambda = lambda;
PInfo.L = L;

% preparing Y and X

Y=DAT(MOP+1:end,:);
X=zeros(T-MOP,p);

for ridx=1:MOP,
X(:,[1:q]+(ridx-1)*q)=DAT(MOP-(ridx-1):end-ridx,:);
end

% performing MPLS

Beta = MM(Y, X, PInfo);

end