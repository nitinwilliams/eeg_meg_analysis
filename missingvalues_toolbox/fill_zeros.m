function D=fill_zeros(A)

D=A;
[M,~]=size(D);

% setting all missing values to 0

D(isnan(A))=0;

% setting all diagonal elements to zero

D(1:M+1:end)=0;

end