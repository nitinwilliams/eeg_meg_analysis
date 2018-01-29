function [Beta, df, RSS, Niter, W] = MM(Y, X, PInfo, Opts)
% [Beta, df, RSS, Niter, W] = MM(Y, X, PInfo, Opts);
%
% This function solves the multivariate least squares model with penalty
% combination:
%        ||y - Xb||^2 + l1*pen1(b) + l2*pen2(b) + ... + lk*penk(b)
% where pen1, pen2, ..., penk is a combination of penalties applied on the
% solution to guarantee convenient properties or satisfy a-priors.
%
% INPUT
% Y:        Data matrix of dimensions dim(N,q). N correspond to data's
%           number of observations and q the number of regressions to do.
% X:        The covariate matrix, dimension is dim(N,p) where p is the
%           number of covariates (coefficients).
% PInfo:    Structure with information about penalties (a-priors).
%           length(PInfo) define the number of penalties combined.
% .L        - the laplace operator applied on coefficients for this prior.
%           Dimension is dim(Nr, p), number of restrictions (Nr) on the p
%           coefficients.
% .lambda   - scalar that impose weight on this prior.
% .stype_penalty
%           - this can be one of {'Ridge', 'l2', 'LASSO', 'l1', 'SCAD',
%           'HardThresh', 'Lq', 'GroupPenalty'}. Only in the case of
%           'GroupPenalty' user must provide PInfo.arr_group and
%           PInfo.sfunc_deriv fields.
% .arr_group
%           - integer vector that define groups to which belong the
%           coefficients. Dimension is dim(p, 1).
% .sfunc_deriv
%           - this can be one of {'Ridge', 'l2', 'LASSO', 'l1', 'SCAD',
%           'HardThresh', 'Lq'}. In this case penalties are applied
%           according to GroupPenalty strategy, see refs.
% Opts:     Structure with information related to the application of the
%           algorithm.
% .b0       - The algorithm here implemented is iterative. This is the
%           initial solution to begin the iteration. Dimension is dim(p, q).
%
% OUTPUT
% Beta:     Coefficient matrix that is the solution of the regression
%           problem defined above, dimension is dim(p, q).
% df:       Degree of freedom. Dimension is dim(1, q).
% RSS:      Is a scalar that represent the residual sum of aquares for all
%           regressions together.
% Niter:    Number of iterations needed for convergence. Dimension is
%           dim(1, q).
% W:        The pseudo-inverse of X operator: Beta = W*Y. It is used to
%           estimate the covariance matrix of coefficients: 
%           var(Beta) = W*var(Y)*W'. For the case of non-linear least
%           squares where Y is a matrix, W is different for each column of
%           Y. Dimension is dim(q, N).
%
% Refs:
% The algorithm implemented is the Minorization-Majorization (MM)
% algorithm, see Hunter & Lange (2004). A tutorial on MM algorithms. Am. Stat. 58, 30�37.
%
% Reference page in Help browser
% <a href="matlab: web('MPLS\Testing\html\html_MM.html', '-helpbrowser')">doc MM</a>

%% Initialization

[N,p] = size(X);
q = size(Y, 2);
Maxit = 200;
tau = 1e-3;
critbeta = 0;
eta = 1e-8;
theta = 0;

%% PARAMETERS & CHECKING

if (N >= p)
    XtY = X'*Y;
    XtX = X'*X;
else
    IN = speye(N);
    Xt = X';
end
% Input options
if ~exist('Opts', 'var')
    Opts = [];
end

%% Initialize output and auxialiar variables for iterations

theta = 0;
Beta = zeros(p, q);
df = zeros(1, q);
Niter = zeros(1, q);
%% Ridge iteration for each column of Y.
% in order to catch the last warning.
lastwarn('');
LASTMSG = '';
for (itq = 1:q)
    % For each regression determined for each column of Y.
    if (N >= p)
        Xy = XtY(:, itq);
    else
        y = Y(:, itq);
    end
    
    % Get or estimate the initial solution.
    if isfield(Opts, 'b0')
        beta = Opts.b0(:, itq);
    else
        Sigma = sparse(p, p);
        for (itp = 1:length(PInfo))
            L = PInfo(itp).L;
            Sigma = Sigma + PInfo(itp).lambda*L'*L;
        end
        if (N >= p)
            %beta =(XtX + Sigma)\Xy;
            beta = pinv(XtX + Sigma)*Xy;
        else
            x = Sigma\Xt;
            % beta = x*((X*x + IN)\y);
            % x=pinv(Sigma)*Xt;
            beta = x*(pinv(X*x + IN)*y);            
        end
    end

    % Iterative procedure.
    if (N >= p)
        for (niter = 1:Maxit)
            oldbeta = beta;
            Sigma = sparse(p, p);
            for (itp = 1:length(PInfo))
                Penalty = PInfo(itp);
                lambda = PInfo(itp).lambda;
                d = eval(PInfo(itp).stype_penalty);
                % Condition d for build Sigma.
                tol = p*eps*max(d);
                d(d<=tol) = (1+eps)*tol;
                Sigma = Sigma + PInfo(itp).L'*spdiags(d, 0, length(d), length(d))*PInfo(itp).L;
            end
            % beta = (XtX + Sigma)\Xy;
            beta = pinv(XtX + Sigma)*Xy;
            
            if (max(abs(oldbeta-beta))./max(max(abs(beta))+eps) < tau), break; end
            [LASTMSG, LASTID] = lastwarn;
            if ~isempty(LASTMSG)
                break;
            end
        end
        
        % Use the equality Tr(A*B) = Tr(B*A) because P is only use to compute the trace.
        % P = (XtX + Sigma)\XtX;
        P = pinv(XtX + Sigma)*XtX;
        
        if (nargout > 4) && (q == 1)
            %W = (XtX + Sigma)\Xt;
            W = pinv(XtX + Sigma)*Xt;
        end
    else % p > N
        for (niter = 1:Maxit)
            oldbeta = beta;
            Sigma = sparse(p, p);
            for (itp = 1:length(PInfo))
                Penalty = PInfo(itp);
                lambda = PInfo(itp).lambda;
                d = eval(PInfo(itp).stype_penalty);
                % Condition d for build Sigma.
                tol = p*eps*max(d);
                d(d<=tol) = (1+eps)*tol;
                Sigma = Sigma + PInfo(itp).L'*spdiags(d, 0, length(d), length(d))*PInfo(itp).L;
            end
            x = Sigma\Xt;
            % beta = x*((X*x + IN)\y);
            % x = pinv(Sigma)*Xt;
            beta = x*(pinv(X*x + IN)*y);
            
            if (max(abs(oldbeta-beta))./max(max(abs(beta))+eps) < tau), break; end
            [LASTMSG, LASTID] = lastwarn;
            if ~isempty(LASTMSG)
                break;
            end
        end
        % P = (X*x)/(X*x+IN);
        P = pinv(X*x)*(X*x+IN);
        if (nargout > 4) && (q == 1)
            % W = x/(X*x + IN);
            W = pinv(x)*(X*x + IN);
        end
    end
    
    Beta(:, itq) = beta;
    df(itq) = trace(P);
    Niter(itq) = niter;
    if ~isempty(LASTMSG)
        break;
    end
end
RSS = mean(sum(abs(Y-X*Beta).^2, 1));

%% Derivative of principal penalties
    function pderiv = penalty_deriv(stype_penalty)
        switch (stype_penalty)
            case 'Ridge'
                pderiv = 2*lambda^2*theta;
            case 'l2'
                pderiv = 2*lambda*theta;
            case {'LASSO', 'l1'}
                pderiv = lambda;
            case 'SCAD'
                a = 3.7; % according to Fan and Li ()
                temp = (((a*lambda-theta)>0).*(a*lambda-theta))*pinv((a-1)*lambda);
                % temp = (((a*lambda-theta)>0).*(a*lambda-theta))/((a-1)*lambda);
                pderiv = lambda*((theta<=lambda)+(temp.*(theta>lambda)));
            case 'HardThresh'
                pderiv = -2*(theta-lambda).*(theta<lambda);
            case 'Lq'
                pderiv = qcoeff*lambda*theta.^(qcoeff - 1);
            otherwise
                error('Penalty unknown.');
        end
    end

%% --------------deriv of penalty functions(theta) ./ theta (theta positive)
% Estimate:
%                      0.5*pderiv(theta)
%               d =  ---------------------
%                         eta + theta

    function result = Ridge
        result = lambda^2*ones(size(Penalty.L, 1), 1);
    end

    function result = l2
        result = lambda*ones(size(Penalty.L, 1), 1);
    end

    function result = LASSO
        theta = abs(Penalty.L*oldbeta);
        pderiv = penalty_deriv('LASSO');
        result = (0.5*pderiv)./(theta+eta);
    end

    function result = l1
        theta = abs(Penalty.L*oldbeta);
        pderiv = penalty_deriv('l1');
        result = (0.5*pderiv)./(theta+eta);
    end

    function result = SCAD
        theta = abs(Penalty.L*oldbeta);
        pderiv = penalty_deriv('SCAD');
        result = (0.5*pderiv)./(theta+eta);
    end

    function result = HardThresh
        theta = abs(Penalty.L*oldbeta);
        pderiv = penalty_deriv('HardThresh');
        result = pderiv./(theta+eta);
    end

    function result = Lq
        theta = abs(Penalty.L*oldbeta);
        pderiv = penalty_deriv('Lq');
        result = pderiv./(theta+eta);
    end

    function result = GroupPenalty
        theta = Penalty.L*oldbeta;
        theta = sqrt(accumarray(Penalty.arr_group, theta.^2));
        theta = theta(Penalty.arr_group);
        pderiv = penalty_deriv(Penalty.sfunc_deriv);
        result = pderiv./(theta+eta);
    end
end