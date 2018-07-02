function PS_ND=partsim_nulldist(A,B,reps)

N=length(B);
PS_ND=zeros(reps,1);

for repidx=1:reps,
    PS_ND(repidx,1)=partsim(A,B(randperm(N))); 
end

end