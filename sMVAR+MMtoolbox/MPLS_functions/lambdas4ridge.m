function lambdas = lambdas4ridge(varargin);
% lambdas = lambdas4ridge(X, npoints, dfmin);
% lambdas = lambdas4ridge(s, N, p, npoints, dfmin);
%
% Find lambdas for a ridge regression problem, in case of multiple linear
% regression problem using a penalization operator L, the user must first
% standardize X (X = X*pinv(L)). It is assumed then that L, of dimensions
% dim(k, p) is full rank p and k >= p.
% This function operate first finding the superior and inferior bounds of
% lambdas. This is based on the closed form of the ridge solution b.
%                   b = inv(X'*X + lambda^2*eye(p))*X'*y
% With this the minimum lambda (lmin) must make X'*X + lmin^2 invertible.
% For estimating the superior bound (lmax) we based on another known
% measure measure: degree of freedom. 
%                df = trace(X*inv(X'*X + lambda^2*eye(p))*X')
% Thus we set lmax such that guarantee a minimun value for df measure.
%                trace(X*inv(X'*X + lmax^2*eye(p))*X') > dfmin
%
% INPUT
% X:        The covariate matrix.
% npoints:  Number of lambda estimated in the interval [lmin, lmax]
%           including the bounds: npoints >= 3.
% dfmin:    The superior bound of lambdas is selected such that the
%           solution had at least dfmin degree of freedom.
%
% OUTPUT
% lambdas:  Lambdas estimated for a ridge regression problem.
%
% We base the lambdas interval estimation in searching first for l2min
% (lmin^2) and l2max (lmax^2).
%
% NeuroStatistic Group, CNC, Cuba.
% May 1st, 2007.
% 
% Reference page in Help browser

%% PARAMETERS & CHECKING
if (nargin == 3)
    X = varargin{1};
    npoints = varargin{2};
    dfmin = varargin{3};
    [N, p] = size(X);
    if (N >= p)
        s2 = svd(full(X'*X));
    else
        s2 = svd(full(X*X'));
    end
    clear X;
elseif (nargin == 5)
    s2 = varargin{1}.^2;
    N = varargin{2};
    p = varargin{3};
    npoints = varargin{4};
    dfmin = varargin{5};
end
% Maximum degree of freedom.
dfmax = min(N, p);
% Minimum singular value
if (N < p)
    minsv = 0;
else
    minsv = min(s2);
end

% X assumed rank.
r = length(s2);
if (dfmin <= 0) || (dfmin >= dfmax)
    error(['The minimum degree of freedom (dfmin) must be 0 < dfmin < ' num2str(dfmax) '.']);
end
if (npoints < 3)
    error('The number of lambdas must be >= 3');
end
lambdas = zeros(npoints, 1);
% Number of bounds interval needed for searching l2min and l2max.
nbounds = 100;
% tolerance value of matrix X'*X
tol = p*eps*s2(1);
% Rank of X
Xrank = sum(s2 > tol);

%% Find bounds for searching l2min: according to make X'*X of full rank
%% if X'*X is full rank already then set to 0 the minimum regularization lambda.
if (Xrank == p)
    l2min = 0;
else
    % search for l2min
    infbound = 0;
    supbound = ((1+eps)*tol)/(1-p*eps);
    for (it = 1:nbounds)
        nextbound = (infbound + supbound)/2;
        newtol = tol + p*eps*nextbound;
        % if we take l2min = nextbound this make X'*X full rank?
        if (minsv+nextbound > newtol)
            supbound = nextbound;
        else
            infbound = nextbound;
        end
    end
    l2min = supbound;
end

%% Find bounds for searching l2max.
infbound = dfmax*s2(1);
infdf = sum(s2./(s2+infbound));
supbound = 0;
supdf = r;
% search for l2max into (supbound < l2max < infbound) guaranteeing that
% infdf < dfmin < supdf
nbounds = max(nbounds, ceil(log2(infbound))+10);
for (it = 1:nbounds)
    nextbound = (infbound + supbound)/2;
    newdf = sum(s2./(s2+nextbound));
    
    if (newdf > dfmin)
        % infdf < dfmin < newdf
        supbound = nextbound;
        supdf = newdf;
    else
        % newdf <= dfmin < supdf
        infbound = nextbound;
        infdf = newdf;
    end
end
l2max = supbound;

%% CHECH LAMBDAS CONSISTENCY
if (l2min > l2max)
    error(['The maximum degree of freedom reachable for this problem is ' num2str(sum(s2./(s2+l2min)))]);
end

%% CONSTRUCT LAMBDAS VECTOR
if (l2min == 0)
    lambdas(1) = l2min;
    lambdas(2) = min(eps*l2max, s2(r));
    lambdas(2:npoints) = logspace(log10(lambdas(2)), log10(l2max), npoints-1);
else
    lambdas(1:npoints) = logspace(log10(l2min), log10(l2max), npoints);
end

lambdas = sqrt(lambdas);
% lambdas = sqrt(lambdas)*sqrt(N);
return;