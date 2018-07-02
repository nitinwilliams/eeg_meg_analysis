clear;
addpath(genpath('B:\Nitin\HBPconnectome\'));

%% Initialisation

% 67 subjects

subnames={'andenna','barrell','bartali','benedetti','berberi','bianchini','bonso','borella','bosi','caravona','catapano',...
    'cioffi','cipriani','cortinovis','depropis','dobrota','farci','fazekas','federico','ferloni','ferranti','ghidini','giordano','ioje',...
    'laface','lori','marchesi','marchesini','marino','martini','menghini','menjivar','meringolo','merlo','micci','migliaccio','minolfi','moni',...
    'morandini','nannarone','nodari','ottoboni','perotti','petriceanu','pintaldi','prezioso','privitera','rocchini','rolla','romano','romanoAlessandra',...
    'roncelli','russo','salimbeni','sava','scotti','selvabonino','sorenti','tarau','testa','tirico','tonelli','ventrone',...
    'veugelers','vivian','zanardo','zanoni'};

Nfvec={'2.50','3.53','5.00','7.07','9.99','14.14','19.99','28.28','39.99','56.56','79.99','113.13','134.54','159.99','190.27','226.27','269.08','319.99'};
reps=100;
thr=0;
M=148;

%% PLV calculation

for Nfidx=5:5, %length(Nfvec),
    
    PLVmats_groupmat=NaN(M,M,reps);
    PLV_numsamps=NaN(M,M,reps);
    
    for repidx=1:reps,
        
        sublist=randi(length(subnames),[length(subnames),1]);
        subnames_surr=subnames(sublist);
        
        PLVmats=cell(1,length(subnames));
        
        for subidx=1:length(subnames);
    
            Nf=Nfvec(Nfidx);
            fname=strcat(Nf{1,1},'Hz_',subnames_surr{subidx},'.mat');
            load(fname,'contactPLV','finalMask');
   
            [nr,~]=size(contactPLV);
   
            contactPLV=contactPLV+tril(contactPLV)';
            contactPLV([1:nr+1:end])=NaN;
            % contactPLV(contactPLV>0.5)=0;
            PLVmats{subidx}=contactPLV;
   
            FMASK{subidx}=finalMask; 
        
        end
        
        %% Parcel numbers
        
        PARClist=cell(1,length(subnames));

        for subidx=1:length(subnames),
        load(['Parc2009_mOp_' subnames_surr{subidx}],'Ax');
            Ax(Ax==-1)=NaN;
            PARClist{1,subidx}=diag(Ax);
        end
        
        %% Morphing to parcellation
        
        for ridx=1:M,
            
            parfor cidx=ridx:M,
                
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
                    
                    PLV_vec=vertcat(PLV_vec,PLVmats{subidx}(findices));
                    
                end
                
                if sum(~isnan(PLV_vec))>thr,
                    
                    PLVmats_groupmat(ridx,cidx,repidx)=nanmean(PLV_vec);
                    PLV_numsamps(ridx,cidx,repidx)=sum(~isnan(PLV_vec));                    
                    
                end
                
            end
            
        end
        
        D1=triu(squeeze(PLVmats_groupmat(:,:,repidx)));
        PLVmats_groupmat(:,:,repidx)=D1+triu(D1,1)';
        
        D2=triu(squeeze(PLV_numsamps(:,:,repidx)));
        PLV_numsamps(:,:,repidx)=D2+triu(D2,1)';
        
        display(repidx);
        
    end
    
    plvfname=strcat('B:\Nitin\SCFC_methods\data\seeg_data\PLV_surrmats_rdtmask_',Nf{1,1},'Hz_thr',int2str(thr),'_',int2str(length(subnames)),'_PARC2k9_wgroupmat.mat');
    save(plvfname,'PLVmats_groupmat','PLV_numsamps');
      
    display(Nfidx);
    
end