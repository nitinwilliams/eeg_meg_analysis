function d=randn_positive(N,mn,sd)

d=normrnd(mn,sd,[N,1]);
nz=sum(d<=0);

while nz,
    d(d<=0)=normrnd(mn,sd,[nz,1]);
    nz=sum(d<=0);   
end

end