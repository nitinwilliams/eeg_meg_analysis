function A_mv=insert_missingvals(A,prop,sMAT,mode) 

% NOTE: assumes prop.>(prop. of initial NaN values)

A=double(A);
UTmat=triu(ones(size(A,1)),1);
validxs=intersect(find(UTmat>0),find(~isnan(A)));
init_prop=sum(isnan(A(:)))/numel(A);
prop=((prop-init_prop)/(1-init_prop));
num_mvals=floor(prop*length(validxs));

switch mode
    
    case 'random'

    % cidxs=randi(length(validxs),[num_mvals,1]);
    perm_indices=randperm(length(validxs));
    cidxs=perm_indices(1:num_mvals);

    case 'samples'
        
    [~,numsamp_indices]=sort(sMAT(validxs),'ascend');
    cidxs=numsamp_indices(1:num_mvals);          
        
end

A=triu(A,1);
A(validxs(cidxs))=NaN;
A_mv=A+A';   

end