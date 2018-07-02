clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% initialising

load('PLVmats_9.83Hz_19_PARC2k9_200parcs_wsubmat.mat','PLVmats_submats');
load('parclist_Destrieux_200to148','parclist');

numparcs=148;
numsubs=size(PLVmats_submats,3);
PLVmats_submats_morphed=zeros(numparcs,numparcs,numsubs);

for subidx=1:numsubs,

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
            
            D=squeeze(PLVmats_submats(:,:,subidx));
            
            PLVmats_submats_morphed(ridx,cidx,subidx)=mean(D(findices));
            PLVmats_submats_morphed(cidx,ridx,subidx)=mean(D(findices));
                                
            end
        
     end

    display(subidx);

end

% saving morphed group matrix

PLVmats_submats=PLVmats_submats_morphed;
save('B:\Nitin\SCFC_methods\data\meg_data\PLVmats_9.83Hz_19_PARC2k9_148parcs_wsubmat.mat','PLVmats_submats');