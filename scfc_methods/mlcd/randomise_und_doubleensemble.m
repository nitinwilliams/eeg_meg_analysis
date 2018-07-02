function A_ND=randomise_und_doubleensemble(A)

% initialise

A_ND=cell(size(A));
N=size(A,1);
T=size(A,2);

M=length(A{1,1});
MASK=triu(ones(M),1);
validxs=find(MASK>0);

for fidx=1:T,
        
    randidxs=validxs(randperm(length(validxs)));
    
    for ridx=1:N,
                  
     D=A{ridx,fidx,1}.*MASK;
     D(validxs)=D(randidxs);
     D=D+D';
     
     A_ND{ridx,fidx,1}=D;
     
     D=A{ridx,fidx,2}.*MASK;
     D(validxs)=D(randidxs);
     D=D+D';
     
     A_ND{ridx,fidx,2}=D;
     
    end
        
end

end