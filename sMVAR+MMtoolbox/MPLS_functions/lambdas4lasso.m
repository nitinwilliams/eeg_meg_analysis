function lambdas = lambdas4lasso(X, nlambdas, dfmin);
% lambdas = lambdas4lasso(X, nlambdas, dfmin);
%
% Find lambdas for a lasso regression problem., in case of multiple linear
% regression problem using a penalization operator L, the user must first
% standardize X (X = X*pinv(L)). It is assumed then that L, of dimensions
% dim(k, p) is full rank p and k >= p.
%
% see lambdas4ridge function.
% >> help lambdas4ridge
%
% INPUT
% X:        The covariate matrix.
% nlambdas: Number of lambda estimated in the interval [lmin, lmax]
%           including the bounds: nlambdas >= 3.
% dfmin:    The superior bound of lambdas is selected such that the
%           solution had at least dfmin degree of freedom.
%
% OUTPUT
% lambdas:  Lambdas estimated for a lasso regression problem.
%
% NeuroStatistic Group, CNC, Cuba.
% February 20, 2006.

%% PARAMETERS & CHECKING...
if (nargin < 3)
    dfmin = 1;
end
% Data dimensions.
if iscell(X)
    N = X{1}(1);
    p = X{1}(2);
else
    [N, p] = size(X);
end

lambdas = lambdas4ridge(X, nlambdas, dfmin);
lambdas = lambdas.^2/sqrt(N);

return;