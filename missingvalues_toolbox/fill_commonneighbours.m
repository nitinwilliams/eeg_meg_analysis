function D=fill_commonneighbours(A,thresh,mode)

D=A;
[M,N]=size(D);
A(1:M+1:end)=NaN;

% finding indices of missing values

indices=find(isnan(A(:)));
[I,J]=ind2sub([M,N],indices);

switch mode
    
    case 'binary'

    % replacing each missing value by 1 or 0 (depending on proportion of common neighbours)

    for id_num=1:length(indices),
    
        row_indices=find(A(I(id_num),:)==1);
        col_indices=find(A(:,J(id_num))==1)';
        cn_val=length(intersect(row_indices,col_indices))/length(union(row_indices,col_indices));
    
        if cn_val>=thresh,
            D(indices(id_num))=1;
        else
            D(indices(id_num))=0;
        end 
        
    end

    case 'weighted'
        
    % replacing each missing value by proportion of common neighbours

    for id_num=1:length(indices),
    
        row_indices=find(A(I(id_num),:)==1);
        col_indices=find(A(:,J(id_num))==1)';
        D(indices(id_num))=length(intersect(row_indices,col_indices))/length(union(row_indices,col_indices));
    
        % if cn_val>=thresh,
        %    D(indices(id_num))=cn_val;
        % else
        %   D(indices(id_num))=0;
        % end 
        
    end
    
    D(isnan(D))=0;
        
end

% setting all diagonal elements to zero
D(1:M+1:end)=0;

end