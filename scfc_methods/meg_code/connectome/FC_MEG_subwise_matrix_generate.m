clear;
addpath(genpath('B:\Nitin\SCFC_methods\'));

%% import group average file

numparcs=200;

subject_sets = {'S0006 set01 PLV 9.83','S0008 set01 PLV 9.83','S0035 set01 PLV 9.83','S0038 set02 PLV 9.83','S0039 set01 PLV 9.83','S0049 set01 PLV 9.83',...
    'S0113 set02 PLV 9.83','S0116 set01 PLV 9.83','S0117 set01 PLV 9.83','S0118 set01 PLV 9.83','S0119 set01 PLV 9.83','S0120 set01 PLV 9.83','S0121 set01 PLV 9.83',...
    'S0123 set01 PLV 9.83','S0124 set01 PLV 9.83','S0126 set01 PLV 9.83','S0127 set01 PLV 9.83','S0128 set01 PLV 9.83','S0130 set01 PLV 9.83'};

PLVmats_submats=zeros(numparcs,numparcs,length(subject_sets));

for subidx=1:length(subject_sets),

fname=strcat('B:\Nitin\SCFC_methods\data\meg_data\single sets\PLV\',subject_sets{subidx},'.csv');
D=abs(dlmread(fname));
D(1:numparcs+1:end)=0;
PLVmats_submats(:,:,subidx)=D;

end

save('B:\Nitin\SCFC_methods\data\meg_data\PLVmats_9.83Hz_19_PARC2k9_200parcs_wsubmat.mat','PLVmats_submats');