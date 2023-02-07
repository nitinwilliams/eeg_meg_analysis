% WilsonCowanNet A network of Wilson-Cowan neural population models.
% The Wilson-Cowan equations describe the mean firing rates of reciprocally
% coupled populations of excitatory and inhibitory neurons. Here we model a
% a network of Wilson-Cowan neural populations where the excitatory populations
% are coupled according to the connectivity matrix Kij.
%
% The differential equations are
%    Ue' = (-Ue + F(wee*Ue - wei*Ui - be + Je + k*Kij*Ue) )./taue;
%    Ui' = (-Ui + F(wie*Ue - wii*Ui - bi + Ji) )./taui;
% where
%    Ue is the mean firing rate of the excitatory cells (nx1)
%    Ui is the mean firing rate of the inhibitory cells (nx1)
%    wei is the weight of the connection to E from I
%    Kij is the network coupling matrix (nxn)
%    k is a scaling constant
%    be and bi are threshold constants
%    Je and Ji are injection currents (either 1x1 or nx1)
%    taue and taui are time constants
%    F(v)=1/(1+\exp(-v)) is a sigmoid function
%
% Returns a sys structure where the number of equations is determined from
% the size of the coupling matrix Kij (which must be nxn).
%
% USAGE:
%    sys = WilsonCowanNet(Kij,Je,Ji)
%
% EXAMPLE:
%    % load the 47x47 CoCoMac connectivity matrix
%    load cocomac047 CIJ
%
%    % Normalise the row-sums of the connectivity matrix
%    Kij = CIJ ./ (sum(CIJ,2) .* ones(47,47));
%
%    % Scalar injection currents
%    Je = -2.5;
%    Ji = -8.5;
%
%    % Construct and run the model
%    sys = WilsonCowanNet(Kij,Je,Ji);
%    gui = bdGUI(sys);
%
% SEE ALSO:
%   WilsonCowan
%   WilsonCowanRing
%
% AUTHOR
%   Stewart Heitmann (2018b)

% Copyright (C) 2016-2019 QIMR Berghofer Medical Research Institute
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in
%    the documentation and/or other materials provided with the
%    distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function sys = WilsonCowanNet_variabledelays_noisy(Kij,Tij,Je,Ji)
    
    % number of nodes
    n = size(Kij,1);
    
    % building lag vector
    lags = [1:max(Tij(:))];
    
    % sanity checks
    assert(numel(Kij)==(n*n), 'Kij must be a square matrix');
    assert(numel(Tij)==(n*n), 'Tij must be a square matrix');
    assert(numel(Je)==1 || numel(Je)==n, 'Je must be 1x1 or nx1');
    assert(numel(Ji)==1 || numel(Ji)==n, 'Ji must be 1x1 or nx1');
    
    % Handle to the DDE function
    sys.ddefun = @ddefun;
    
    % Determine the best limits for Kij    
    % KijLim = bdPanel.RoundLim(min(Kij(:)),max(Kij(:)));
    
    % DDE parameters
    sys.pardef = [ struct('name','wee', 'value', 0,  'lim',[0 30]);
                   struct('name','wei', 'value', 0,  'lim',[0 30]);
                   struct('name','wie', 'value', 0,  'lim',[0 30]);
                   struct('name','wii', 'value', 0,   'lim',[0 30]);
                   struct('name','k',  'value', 0,   'lim',[0 5]);
                   struct('name','Kij', 'value', Kij, 'lim',[0,1]);
                   struct('name','Tij', 'value', Tij, 'lim',[0,100]);
                   struct('name','be',  'value', 0, 'lim',[0 10]);
                   struct('name','bi',  'value', 0,   'lim',[0 10]);
                   struct('name','Je',  'value', 0,  'lim',[0 5]); 
                   struct('name','Ji',  'value', 0,  'lim',[0 5]);
                   struct('name','taue', 'value', 0,   'lim',[1 20]);
                   struct('name','taui', 'value', 0,   'lim',[1 20]);
                   struct('name','psi_sigma', 'value', 0,   'lim',[0,1])];
              
    % DDE variables
    sys.vardef = [ struct('name','Ue', 'value',rand(n,1), 'lim',[0 1]);
                   struct('name','Ui', 'value',rand(n,1), 'lim',[0 1])];
                              
    % DDE lag parameters
    sys.lagdef = [ struct('name','lag', 'value',lags) ];
 
    % Default time span
    sys.tspan = [0 100];
    
    % Default DDE options
    sys.ddeoption.RelTol = 1e-3;
    
    % Latex Panel
    sys.panels.bdLatexPanel.latex = {
        '\textbf{WilsonCowanNet}'
        ''
        'A network of Wilson-Cowan equations where the nodes of the network'
        'are local populations of excitatory and inhibitory neurons. Only the'
        'excitatory cells are connected by the network. The inhibitory interact-'
        'ions are local only. The dynamical equations are'
        ''
        '\qquad $\tau_e \; \dot U_e = -U_e + F\big(w_{ee} U_e - w_{ei} U_i - b_e + J_e + k \sum_j K_{ij} U_e\big)$'
        '\qquad $\tau_i \; \dot U_i\; = -U_i \; + F\big(w_{ie} U_e - w_{ii} U_i - b_i + J_i \big)$'
        ''
        'where'
        '\qquad $U_e$ is the firing rate of the \textit{excitatory} populations (nx1),'
        '\qquad $U_i$ is the firing rate of the \textit{inhibitory} populations (nx1),'
        '\qquad $w_{ei}$ is the weight of the connection to $e$ from $i$,'
        '\qquad $K_{ij}$ in an nxn connectivity matrix,'
        '\qquad $T_{ij}$ in an nxn delay matrix,'
        '\qquad $k$ is a scaling constant,'
        '\qquad $b_{e}$ and $b_{i}$ are threshold constants,'
        '\qquad $J_{e}$ and $J_i$ are injection currents (1x1 or nx1),'
        '\qquad $\tau_{e}$ and $\tau_{i}$ are time constants,'
        '\qquad $F(v)=1/(1+\exp(-v))$ is a sigmoidal firing-rate function,'
        ''
        '\textbf{References}'
        'Wilson \& Cowan (1972) Biophysics Journal 12(1):1-24.'
        'Hlinka \& Coombes (2012) Euro J Neurosci 36:2137-2145.'
        };
    
    % Other Panels
    sys.panels.bdTimePortrait = [];
    sys.panels.bdSpaceTime = [];
    sys.panels.bdSolverPanel = [];
end

function dU = ddefun(~,U,Z,wee,wei,wie,wii,k,Kij,Tij,be,bi,Je,Ji,taue,taui,psi_sigma)
    
    % Extract incoming data
    n = numel(U)/2;
    Ue = U(1:n);
    Ui = U([1:n]+n);
    Ze = Z(1:n,:)';
    
    % Building set of time-delayed values
    Ue_delaymat=zeros(n);
    for nidx=1:n
        source_idxs=setdiff(1:n,nidx);
        delays=Tij(nidx,source_idxs);
        delay_idxs=((source_idxs-1)*size(Ze,1))+delays;
        Ue_delaymat(source_idxs,nidx)=Ze(delay_idxs);      
    end
    
    % Wilson-Cowan dynamics
    dUe = (-Ue + F(wee*Ue - wei*Ui + Je - be + psi_sigma*randn + k*diag(Kij*Ue_delaymat)) )./taue;
    dUi = (-Ui + F(wie*Ue - wii*Ui + Ji - bi + psi_sigma*randn) )./taui;
    
    % concatenate results
    dU = [dUe;dUi];
    
end

% Sigmoidal firing-rate function
function y = F(x)
    y = 1./(1+exp(-x));
end
