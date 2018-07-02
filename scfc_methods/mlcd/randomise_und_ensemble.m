function A_ND=randomise_und_ensemble(A)

% initialise

A_ND=cell(size(A));
N=size(A,1);
T=size(A,2);

M=length(A{1,1});
MASK=triu(ones(M),1);
validxs=find(MASK>0);

% randidxs=validxs(randperm(length(validxs)));

for fidx=1:T,
    
    randidxs=validxs(randperm(length(validxs)));
    
    for ridx=1:N,
                  
     D=A{ridx,fidx}.*MASK;
     D(validxs)=D(randidxs);
     D=D+D';
     
     A_ND{ridx,fidx}=D;
     
    end
        
end

end