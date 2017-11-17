function pTE=phaseTE(Xf,lag)
% Inputs:
% Xf is 2D matrix, i.e. channels x samples
% lag is the expected delay, expressed in samples
% Outputs:
% pTE(1,2) is phase TE from 1 to 2
% pTE(2,1) is phase TE from 2 to 1
%
% Function written by Nitin Williams, with help from Felix Siebenhuehner and Muriel Lobier 
%
%
%% Inst. phase

Xfh=angle(hilbert(Xf'))';

%% Number of bins

for idx=1:size(Xfh,1),
    [~,dirsd(idx,1)]=circ_std(Xfh(idx,:),[],[],2);
end

N=size(Xf,2);
bw=(((3.49).*dirsd)./(N^(1/3)));
numbins=round((2*pi)./bw);

%% Computing normalised histograms

tX=Xfh(1,1+lag:end);
tY=Xfh(2,1+lag:end);
tXd=Xfh(1,1:end-lag);
tYd=Xfh(2,1:end-lag);

tXh=discretize(tX,linspace(-pi,pi,numbins(1,1)));
tYh=discretize(tY,linspace(-pi,pi,numbins(2,1)));
tXdh=discretize(tXd,linspace(-pi,pi,numbins(1,1)));
tYdh=discretize(tYd,linspace(-pi,pi,numbins(2,1)));

[P_YYd,~,~]=crosstab(tYh,tYdh); P_YYd=P_YYd./sum(P_YYd(:));
[P_YdXd,~,~]=crosstab(tYdh,tXdh); P_YdXd=P_YdXd./sum(P_YdXd(:));
[P_Yd,~,~]=crosstab(tYdh); P_Yd=P_Yd./sum(P_Yd(:));
[P_YYdXd,~,~]=crosstab(tYh,tYdh,tXdh); P_YYdXd=P_YYdXd./sum(P_YYdXd(:));

[P_XXd,~,~]=crosstab(tXh,tXdh); P_XXd=P_XXd./sum(P_XXd(:));
[P_XdYd,~,~]=crosstab(tXdh,tYdh); P_XdYd=P_XdYd./sum(P_XdYd(:));
[P_Xd,~,~]=crosstab(tXdh); P_Xd=P_Xd./sum(P_Xd(:));
[P_XXdYd,~,~]=crosstab(tXh,tXdh,tYdh); P_XXdYd=P_XXdYd./sum(P_XXdYd(:));

%% Entropy & pTE calculation

H_YYd = sum(-(P_YYd(P_YYd>0).*(log2(P_YYd(P_YYd>0)))));
H_YdXd = sum(-(P_YdXd(P_YdXd>0).*(log2(P_YdXd(P_YdXd>0)))));
H_Yd = sum(-(P_Yd(P_Yd>0).*(log2(P_Yd(P_Yd>0)))));
H_YYdXd = sum(-(P_YYdXd(P_YYdXd>0).*(log2(P_YYdXd(P_YYdXd>0)))));

pTE(1,1)=NaN;
pTE(2,2)=NaN;
pTE(1,2)=H_YYd+H_YdXd-H_Yd-H_YYdXd;

H_XXd = sum(-(P_XXd(P_XXd>0).*(log2(P_XXd(P_XXd>0)))));
H_XdYd = sum(-(P_XdYd(P_XdYd>0).*(log2(P_XdYd(P_XdYd>0)))));
H_Xd = sum(-(P_Xd(P_Xd>0).*(log2(P_Xd(P_Xd>0)))));
H_XXdYd = sum(-(P_XXdYd(P_XXdYd>0).*(log2(P_XXdYd(P_XXdYd>0)))));

pTE(2,1)=H_XXd+H_XdYd-H_Xd-H_XXdYd;

end