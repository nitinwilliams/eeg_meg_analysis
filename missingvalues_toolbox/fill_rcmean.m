function D=fill_rcmean(A,thresh,mode)

D=A;
[M,N]=size(D);
A(1:M+1:end)=NaN;

% finding indices of missing values

indices=find(isnan(A(:)));
[I,J]=ind2sub([M,N],indices);

switch mode
    
    case 'binary'

    % replacing each missing value by 1 or 0

    for id_num=1:length(indices),
    
        rowcol_vec=[A(I(id_num),:),A(:,J(id_num))'];
        mn_val=nanmean(rowcol_vec);
    
        if mn_val>=thresh,
            D(indices(id_num))=1;
        else
            D(indices(id_num))=0;
        end 
        
    end
    
    case 'weighted'
        
    % replacing each missing value by mean or 0

    for id_num=1:length(indices),
    
        rowcol_vec=[A(I(id_num),:),A(:,J(id_num))'];
        D(indices(id_num))=nanmean(rowcol_vec);
        
        % if mn_val>=thresh,
        %    D(indices(id_num))=mn_val;
        % else
        %    D(indices(id_num))=0;
        % end 
        
    end
        
end

% setting all diagonal elements to zero

D(1:M+1:end)=0;

end