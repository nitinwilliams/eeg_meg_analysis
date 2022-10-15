function [D_nc]=fill_missingvals(D,numc) 

D_nc={numc,1};
A=triu(ones(size(D,1)),1);

nvalidxs=intersect(find(isnan(D)),find(A>0));
validxs=intersect(find(~isnan(D)),find(A>0));

for ncidx=1:numc,

 Dm=triu(D,1);
 cidxs=randi(length(validxs),[length(nvalidxs),1]);
 Dm(nvalidxs)=D(validxs(cidxs));
 Dm=Dm+Dm';   
 D_nc{ncidx}=Dm;
    
end

end