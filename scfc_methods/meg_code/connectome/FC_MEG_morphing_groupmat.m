clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% initialsing

load('PLVmats_9.83Hz_19_PARC2k9_200parcs_wgroupmat.mat','PLVmats_groupmat');
load('parclist_Destrieux_200to148','parclist');

numparcs=148;
PLVmats_groupmat_morphed=zeros(numparcs);

for ridx=1:numparcs,
   
    for cidx=ridx+1:numparcs,
                                               
            A=zeros(length(parclist),1); B=zeros(length(parclist),1);
                        
            RN=ridx-1;
            CN=cidx-1;
                        
            cindices1=(parclist==(RN));
            cindices2=(parclist==(CN));
            
            A(cindices1,1)=1;
            B(cindices2,1)=1;
            
            C=A*B';
            findices=(C(:)==1);
            
            PLVmats_groupmat_morphed(ridx,cidx)=mean(PLVmats_groupmat(findices));
            PLVmats_groupmat_morphed(cidx,ridx)=mean(PLVmats_groupmat(findices));
                                
    end
    
    display(ridx);
    
end

% saving morphed group matrix

PLVmats_groupmat=PLVmats_groupmat_morphed;
save('B:\Nitin\SCFC_methods\data\meg_data\PLVmats_9.83Hz_19_PARC2k9_148parcs_wgroupmat.mat','PLVmats_groupmat');