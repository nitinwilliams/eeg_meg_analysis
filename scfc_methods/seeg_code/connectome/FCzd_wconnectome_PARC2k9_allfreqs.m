clear;
rmpath(genpath('B:\Nitin\HBPconnectome\'));
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% Initialisation

% 67 subjects

subnames={'andenna','barrell','bartali','benedetti','berberi','bianchini','bonso','borella','bosi','caravona','catapano',...
    'cioffi','cipriani','cortinovis','depropis','dobrota','farci','fazekas','federico','ferloni','ferranti','ghidini','giordano','ioje',...
    'laface','lori','marchesi','marchesini','marino','martini','menghini','menjivar','meringolo','merlo','micci','migliaccio','minolfi','moni',...
    'morandini','nannarone','nodari','ottoboni','perotti','petriceanu','pintaldi','prezioso','privitera','rocchini','rolla','romano','romanoAlessandra',...
    'roncelli','russo','salimbeni','sava','scotti','selvabonino','sorenti','tarau','testa','tirico','tonelli','ventrone',...
    'veugelers','vivian','zanardo','zanoni'};

thr=0;
numparcs=(148/2);
lefthem_idxs=[1:2:numparcs*2];
righthem_idxs=[2:2:numparcs*2];

Nfvec={'2.50','3.53','5.00','7.07','9.99','14.14','19.99','28.28','39.99','56.56','79.99','113.13','134.54','159.99','190.27','226.27','269.08','319.99'};

%% PLV calculation

for Nfidx=1:6, %length(Nfvec),

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
   
plvfname=strcat('B:\Nitin\SCFC_methods\data\seeg_data\PLVmats_rdtmask_',Nf{1,1},'Hz_thr',int2str(thr),'_',int2str(length(subnames)),'_PARC2k9_wgroupmat.mat');
save(plvfname,'PLVmats_groupmat','PLV_numsamps');

display(Nfidx);

end