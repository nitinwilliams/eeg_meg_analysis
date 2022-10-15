clear;

%% Initialisation

% 67 subjects

subnames={'subject1','subject2','subject3','subject4','subject5','subject6','subject7',...
    'subject8','subject9','subject10','subject11','subject12','subject13','subject14',...
    'subject15','subject16','subject17','subject18','subject19','subject20','subject21',...
    'subject22','subject23','subject24','subject25','subject26','subject27','subject28',...
    'subject29','subject30','subject31','subject32','subject33','subject34','subject35',...
    'subject36','subject37','subject38','subject39','subject40','subject41','subject42',...
    'subject43','subject44','subject45','subject46','subject47','subject48','subject49',...
    'subject50','subject51','subject52','subject53','subject54','subject55','subject56',...
    'subject57','subject58','subject59','subject60','subject61','subject62','subject63',...
    'subject64','subject65','subject66','subject67'};

thr=0;
numparcs=(148/2);
lefthem_idxs=[1:2:numparcs*2];
righthem_idxs=[2:2:numparcs*2];

Nfvec={'2.50','3.53','5.00','7.07','9.99','14.14','19.99','28.28','39.99','56.56','79.99','113.13','134.54','159.99','190.27','226.27','269.08','319.99'};

%% PLV calculation

for Nfidx=1:length(Nfvec),

    for subidx=1:length(subnames);
    
        Nf=Nfvec(Nfidx);
        fname=strcat(Nf{1,1},'Hz_',subnames{subidx},'.mat');
        load(fname,'contactPLV','finalMask');
   
        [nr,~]=size(contactPLV);
   
        contactPLV=contactPLV+tril(contactPLV)';
        contactPLV([1:nr+1:end])=NaN;
        PLVmats{subidx}=contactPLV;
   
        FMASK{subidx}=finalMask; 
        
     end

%% Parcel numbers

    PARClist=cell(1,length(subnames));

    for subidx=1:length(subnames),
        load(['Parc2009_mOp_' subnames{subidx}],'Ax');
        Ax(Ax==-1)=NaN;
        PARClist{1,subidx}=diag(Ax);
    end

%% Morphing to parcellation

M=148;
PLVmats_groupmat=NaN(M);
PLV_numsamps=NaN(M);
       
for ridx=1:M,
        
   for cidx=ridx:M,
       
       PLV_vec=[];
            
        for subidx=1:length(subnames),
            
            MASK=FMASK{subidx};
            
            parclist = PARClist{1,subidx};
            A=zeros(length(parclist),1); B=zeros(length(parclist),1);
            
            % PARCELS
            
            RN=ridx-1;
            CN=cidx-1;
                        
            cindices1=(parclist==(RN));
            cindices2=(parclist==(CN));
            
            A(cindices1,1)=1;
            B(cindices2,1)=1;
            
            C=A*B';
            C=C.*MASK;
            findices=(C(:)==1);
            
            PLV_vec=vertcat(PLV_vec,PLVmats{1,subidx}(findices));
            
        end
        
        if sum(~isnan(PLV_vec))>thr,
        
            PLVmats_groupmat(ridx,cidx)=nanmean(PLV_vec);
            PLVmats_groupmat(cidx,ridx)=nanmean(PLV_vec);
            
            PLV_numsamps(ridx,cidx)=sum(~isnan(PLV_vec));
            PLV_numsamps(cidx,ridx)=sum(~isnan(PLV_vec));
            
        end
            
    end
        
end
   
plvfname=strcat('C:\Users\willian1\OneDrive - Aalto University\Previous_projects\Nitin\HBPconnectome\GA_data\submission\PLVmats_rdtmask_',Nf{1,1},'Hz_thr',int2str(thr),'_',int2str(length(subnames)),'subs_parc2k9_wgroupmat.mat');
save(plvfname,'PLVmats_groupmat','PLV_numsamps');

display(Nfidx);

end