function A_ND=randomise_und_ensemble_nodeordering_ms(A)

% initialise

A_ND=cell(size(A));
N=size(A,1);
T=size(A,2);

M=length(A{1,1});

for fidx=1:T,
    
    reorder=randperm(M);
           
    for ridx=1:N,
                  
     D=A{ridx,fidx};
     D=D(reorder,reorder);
     
     A_ND{ridx,fidx}=D;
     
    end
        
end

end