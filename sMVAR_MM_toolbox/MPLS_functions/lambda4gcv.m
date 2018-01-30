function [beta, lsel, critsel, critvals, ind] = lambda4gcv(Betas, lambdas, N, RSS, df, Opts)
% [beta, lsel, critsel, critvals] = lambda4gcv(Betas, lambdas, Opts);
%
% This function implement the GCV-function
%                      ||y - Xb||^2         RSS
%              GCV = ---------------- = ----------- 
%                      (tr(I - H))^2     (N - df)^2
% where H is the 'hat' matrix, y_hat = H*y.
%
% INPUT
% Betas:    Cell array of solutions from which user wish to select the best
%           one, attending to GCV criterion. These solutions correspond to
%           different models (nlambdas) applied to the regression problem
%           with multiples combined a-prioris. Cell dimension is
%           dim(1, nlambdas). If the user is only interested in getting the
%           'best' lambda, he can input Betas==[].
% lambdas:  Weight assigned to each combined a-priori for giving strength
%           to this assumption. It is a 2D matrix or a column vector:
%           rows corresponding to the different models tested (nlambdas)
%           and columns to the a-prioris combined. Dimension is
%           dim(nlambdas, npenalties).
% N         Number of independent observations for the data. More exactly,
%           the rank of the covariate matrix X or the highest degree of
%           freedom reachable for any model.
% RSS       Residual sum of squares for each adjusted model. Dimension is
%           dim(nlambdas, 1).
% df        Measure of 'degree of freedom' reached by each adjusted
%           model. Dimension is dim(nlambdas, 1).
% Opts:     Structure whose fields are options used by lambda4gcv for
%           computing GCV criterion. These are:
% .f_graphic
%           - equals 1 to plot curve, 0 instead.
%
% OUTPUT
% beta:     If Betas~=[], then beta is the solution selected by the user,
%           else it is empty.
% lsel:     The lambda selected to minimize the GCV criterion or one
%           selected attending to an user subjective criterion.
% critsel:  The criterion value that correspond to lsel.
% critvals: The criterion values estimated for GCV selection criterion.
% ind:      index selected: lsel = lambdas(ind).

%% PARAMETERS & CHECKING...
if ~exist('Opts', 'var')
    Opts = [];
end
% local lambda selected?
flocal = 0; if isfield(Opts, 'flocal'), flocal = Opts.flocal; end
% Graphic option.
f_graphic = 0; if isfield(Opts, 'f_graphic'), f_graphic = Opts.f_graphic; end

%% critvals for GCV.
critvals = RSS./(N-df).^2;

%% Select lsel and critsel according to user selection.
if flocal
    ind = localminima(critvals);
else
    [critsel, ind] = min(critvals);
end
if isempty(ind)
    if (critvals(1) < critvals(end))
        ind = 1;
    else
        ind = length(critvals);
    end
else
    ind = find(critvals == min(critvals(ind)));
    ind = ind(1);
end
critsel = critvals(ind);
lsel = lambdas(ind, :);

if isempty(Betas)
    beta = [];
else
    beta = Betas{ind};
end

%% Plot cross-validation curve and lambda selected.
if (f_graphic == 1)
    if (size(lambdas, 2) > 1)
        lambdas = sum(lambdas, 2);
    end
    loglog(lambdas, critvals, sum(lsel), critsel, '*r');
    xlabel('\lambda');
    ylabel(['GCV' '(\lambda)']);
    title(['GCV' ' function, selection at \lambda = ' num2str(sum(lsel))]);
end
return;