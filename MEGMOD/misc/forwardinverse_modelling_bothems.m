function D_linmix=forwardinverse_modelling_bothems(D,source2parc,py2lv,hem_idxs,fwd_mat,inv_mat)

% translating source to parcel vector

source2parc_trans=zeros(size(source2parc));
source2parc_trans(source2parc==-1)=-1;
for parc_idx=1:length(py2lv)
    source2parc_trans(source2parc==(parc_idx-1))=find(py2lv==parc_idx);       
end

% duplicating original dataset

num_sources=sum(source2parc<=73)+sum(source2parc>73)-sum(source2parc==-1);
D_sources=zeros(num_sources,size(D,2));
fwd_mat_reordered=zeros(size(fwd_mat,1),num_sources);
inv_mat_reordered=zeros(num_sources,size(inv_mat,2));

curr_idx=1;
for idx=1:length(hem_idxs)
    M=sum(source2parc_trans==hem_idxs(idx));
    D_sources([curr_idx:curr_idx+M-1],:)=repmat(D(idx,:),[M,1]);
    fwd_mat_reordered(:,[curr_idx:curr_idx+M-1])=fwd_mat(:,source2parc_trans==hem_idxs(idx));
    inv_mat_reordered([curr_idx:curr_idx+M-1],:)=inv_mat(source2parc_trans==hem_idxs(idx),:);
    curr_idx=curr_idx+M;
end

% forward-inverse matrix multiplication

D_linmix_sources=inv_mat_reordered*(fwd_mat_reordered*D_sources);

% averaging within parcels

row_num=1;
D_linmix=zeros(size(D));
for parc_idx=1:size(D,1)
    num_rows=sum(source2parc_trans==hem_idxs(parc_idx));
    D_linmix(parc_idx,:)=mean(D_linmix_sources([row_num:(row_num+num_rows-1)],:));
    row_num=row_num+num_rows;
end

end