function  [D_pair,V_pair]=network_shuffle(D,V,pctval)

N=length(V);

% find number of nodes to shuffle
num_nodes=round((pctval/100)*N);

% obtain random indices
rand_order=randperm(N);
[rand_nodes,~]=sort(rand_order(1:num_nodes),'ascend');

shuf_order1=randperm(num_nodes);
shuf_order2=randperm(num_nodes);

orig_order=[1:N];

reorder1=orig_order;
reorder1(rand_nodes)=rand_nodes(shuf_order1);

reorder2=orig_order;
reorder2(rand_nodes)=rand_nodes(shuf_order2);

% perform shuffling

D_pair=cell(1,2);
V_pair=cell(1,2);

D_pair{1}=D(reorder1,reorder1);
D_pair{2}=D(reorder2,reorder2);

V_pair{1}=V(reorder1);
V_pair{2}=V(reorder2);

end
