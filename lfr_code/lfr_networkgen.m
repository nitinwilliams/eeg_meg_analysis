clear;
addpath(genpath('/home/nw04/fcmodules_codedata/Nitin/'));

%% generating LFR benchmark network

N=74;
k=10;		
maxk=20;
t1=2;		
t2=1;		
minc=6;	
maxc=12;	
on=0;		
om=0;		

muvec=[0.9];
reps=1000;

A=cell(reps,length(muvec));
CP=cell(reps,length(muvec));

for muidx=1:length(muvec),
    
mu=muvec(muidx);

for repidx=1:reps,

unix(['./benchmark -N ' num2str(N) ' -k ' num2str(k) ' -maxk ' num2str(maxk) ' -mu ' num2str(mu) ' -t1 ' num2str(t1) ' -t2 ' num2str(t2) ' -minc ' num2str(minc) ' -maxc ' num2str(maxc) ' -on ' int2str(on) ' -om ' int2str(om)]);

[A{repidx,muidx},CP{repidx,muidx}]=LFR_networks(N,'binary');

display(repidx);

end

display(muidx);

end

save('/home/nw04/fcmodules_codedata/Nitin/SCFC_methods/data/LFR_data/LFR_binary_networks_mixfac0pt9_SCFCmethods.mat','A','CP');

% Visualising modules

D=A{1,1};
V=CP{1,1};

[xynew,~,~]=FRKK(sparse(D),V,0.01,9999,2,4); 

figure;
graphplot2d(xynew,D,2,V);
axis square;
set(gca,'Visible','off');