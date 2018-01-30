function [Betas, lambdas, OUTPUT] = penalty4lambdas(Y, X, lambdas, PInfo, Opts)
% [Betas, lambdas, OUTPUT] = penalty4lambdas(Y, X, lambdas, PInfo, Opts);
%
% This function resolves automatically the solution to the multivariate
% linear regression problem with combination of penalties (>> help MM).
%
% INPUT
% Y:        Data matrix of dimensions dim(N, m). N correspond to data's
%           number of observations and m to number of regressions to do.
%
% X:        The covariate matrix, dimension is dim(N, p) where p is the
%           number of covariates.
%
% lambdas:  If user wish to input its own lambdas (USER_LAMBDAS option) he
%           must set this parameter. lambdas must be a 2D matrix with
%           dimension dim(nlambdas, npenalties) where nlambdas correspond
%           to the number of lambdas that user wish to compute and
%           npenalties to the number of penalties combined. If lambdas==[]
%           then only automatic (AUTOMATIC_SELECTION option) or proportion
%           selection (LAMBDA_PROPORTION) are possible.
%
% PInfo:    Structure with information about penalties (a-priors).
%           length(PInfo) define the number of penalties combined, see
%           >> help MM.
%
% Opts:     Structure with information related to the application of the
%           MM's algorithm and options to execute at each iteration, see
%           >> help MM for some field definitions. Other fields are:
% .IC       - equals 1 iff coefficient initialization is defined, this
%           refer to use an initial solution for compute solution at a
%           particular lambda, this initial solution is the final solution
%           computed for previous lambda. Note that this initialization
%           step could be better or worst in dependece with the grid of
%           lambdas, if it is more fine or thick, respectively.
%           equals 0 for use Ridge aproximation as initial solution.
% .verbose  - equals 1 to display iteration number for impatient users or
%           huge data that take long time to execute.
% .mu       - this parameter set LAMBDA_PROPORTION option and it is used
%           when multiple penalties (npenalties > 1) are combined.
%           Dimension is dim(1, npenalties). Proportions are in (0, 1)
%           interval and add 1.
% .nlambdas - number of lambdas defined by user, this parameter is only
%           admissible when lambdas==[]. If lambdas is not empty
%           nlambdas==length(lambdas), default is 50.
% .hmethod   - the algorithm used to compute the model, one of {@MM @LQA}.
%
% OUTPUT
% Betas:    Cell array which contains the solution of the multivariate
%           regression for each lambdas. Cell dimension is dim(1, nlambdas).
% lambdas:  The lambdas entered by user or selected automatically by this
%           program attending to input parameters. Dimension is
%           dim(nlambdas, npenalties).
% OUTPUT:   Structure whose fields are needed for compute statistical
%           measures as AIC, BIC, GCV, ...
% .df       - degree of freedom. Dimension is dim(nlambdas, 1).
% .RSS      - Residual sum squares of the regression problem. Dimension is dim(nlambdas, 1).
% .sigma2   - error variance estimation for a model whose error term has 
%           gaussian distribution with parameters (0, sigma2*eye(N)).
%           Dimension is dim(nlambdas, 1)
% 
% Reference page in Help browser
%   <a href="matlab: web('MPLS\Testing\html\html_penalty4lambdas.html', '-helpbrowser')">doc penalty4lambdas</a>

%% PARAMETERS & CHECKING...
% Data dimensions.
[N, p] = size(X);
m = size(Y, 2);
% Number of penalties.
npenalty = length(PInfo);
hmethod = [];
verbose = 0;
check_input;

%% Set lambda option according to user input.
% Constants.
USER_LAMBDAS = 1;
LAMBDA_PROPORTION = 2;
AUTOMATIC_SELECTION = 3;
[option, nlambdas, mu] = lambda_option;

%% lambda selection automatically.
% if isempty(lambdas)
%     lambdas = lambdas_automatically;
% end
lambdas = lambdas_automatically;

%% Compute solutions for each lambda.
% beta, degree of freedom and RSS are recovered from the engine for each lambda.
Betas = {};
RSS = zeros(nlambdas, 1);
for (it_pinfo = 1:npenalty)
    PInfo(it_pinfo).lambda = lambdas(nlambdas, it_pinfo);
end
if INIT_COEFF
    % dfraw has dimensions dim(nboots, m, nlambdas).
    [Betas{nlambdas}, dfraw(:, :, nlambdas), RSS(nlambdas), Niter] = feval(hmethod, Y, X, PInfo, OptsAlg);
else
    % dfraw has dimensions dim(1, m, nlambdas).
    [beta, dfraw(:, :, nlambdas), RSS(nlambdas), Niter] = feval(hmethod, Y, X, PInfo, OptsAlg);
end
if (verbose == 2)
    %h = waitbar(0,'Please wait...');
    guidlg_time('inicializar');
end
for (it = nlambdas-1:-1:1)
    if (verbose == 1)
        disp(it);
    elseif (verbose == 2)
        %waitbar(, h);
        guidlg_time((nlambdas-it)/nlambdas);
        drawnow;
    end
    for (it_pinfo = 1:npenalty)
        PInfo(it_pinfo).lambda = lambdas(it, it_pinfo);
    end
    if INIT_COEFF
        OptsAlg.b0 = Betas{it+1};
        [Betas{it}, dfraw(:, :, it), RSS(it), Niter] = feval(hmethod, Y, X, PInfo, OptsAlg);
    else
        [beta, dfraw(:, :, it), RSS(it), Niter] = feval(hmethod, Y, X, PInfo);
    end
    % Check if the engine had some warning or kind of error, break if TRUE.
    [LASTMSG, LASTID] = lastwarn;
	if ~isempty(LASTMSG)
        if CONJUGATED_GRADIENT
            error('This warning with CONJUGATE_GRADIENT option is unexpected!');
        end
        dfraw(:, :, 1:it) = repmat(dfraw(:, :, it+1), [1 1 it]);
        RSS(1:it) = 1.1*RSS(it+1);
        if INIT_COEFF
            vecNaN = ones(size(Betas{it+1}))*NaN;
            for (itB = 1:it), Betas{itB} = vecNaN; end
        end
        break;
    end
end
if (verbose == 2), guidlg_time('finalizar'); end

% Postprocessing of dfraw.
if CONJUGATED_GRADIENT
    for (it = m:-1:1)
        dfest = squeeze(mean(dfraw(:, it, :)));
        dfest = monofit([log(lambdas) log(dfest)], [0 2 1 0]);
        df(:, it) = exp(dfest(:, 2));
    end
    %df = sum(df, 2);
    df = mean(df, 2);
    if (df(1) > min(N, p))
        dfnum = min(N, p)-1e-8 - df(end);
        dfden = df(1)-df(end);
        df = df(end) + (dfnum/dfden)*(df - df(end));
    end
    OUTPUT.df = df;
else
    %OUTPUT.df = squeeze(sum(dfraw, 2));
    OUTPUT.df = squeeze(mean(dfraw, 2));
end
if isreal(Y) && isreal(X)
    OUTPUT.RSS = RSS;
else
    OUTPUT.RSS = RSS/2;
end
%OUTPUT.sigma2 = RSS./(m*N-OUTPUT.df);
OUTPUT.sigma2 = OUTPUT.RSS./(N-OUTPUT.df);

%% Check input parameters
    function check_input
        % Check PInfo struct
        if ~isfield(PInfo, 'stype_penalty')
            error('The struct PInfo must define the ''stype_penalty'' field.');
        end
        if ~isfield(PInfo, 'L')
            PInfo(1).L = [];
        end
        for (it_pinfo = 1:length(PInfo))
            if isempty(PInfo(it_pinfo).L), PInfo(it_pinfo).L = speye(p); end
        end
        % check Opts.
        if ~exist('Opts')
            Opts = [];
        end
        % Is Regression method defined?
        if isfield(Opts, 'hmethod')
            hmethod = Opts.hmethod;
        end
        % Is CONJUGATED_GRADIENT defined?
        if isfield(Opts, 'CG')
            CONJUGATED_GRADIENT = Opts.CG;
            OptsAlg.CG = Opts.CG;
        else
            CONJUGATED_GRADIENT = 0;
            OptsAlg = [];
        end
        % is coefficient initialization defined?
        if isfield(Opts, 'IC')
            INIT_COEFF = Opts.IC;
        else
            INIT_COEFF = 0;
        end
        % Verbose option.
        if isfield(Opts, 'verbose'), verbose = Opts.verbose; end
    end

%% Lambda option according to user input.
    function [option, nlambdas, mu] = lambda_option
        % What option for selecting lambdas to choice?
        option = [];
        % Number of lambdas.
        nlambdas = [];
        % Number of proportion.
        mu = [];
        if ~isempty(lambdas)
            % User provide its own lambdas.
            option = USER_LAMBDAS;
            % In this case the user must provide a square matrix lambdas, with so
            % many columns as penalties are tested in combination. The number of
            % rows correspond to the number of critvals estimated values for each
            % combination of lambdas and penalties.
            nlambdas = size(lambdas, 1);
            if (size(lambdas, 2) ~= npenalty)
                error(['The number of columns of ''lambdas'' must coincide with the ' ...
                    'number of penalties.']);
            end
        end
        if isfield(Opts, 'mu')
            % User provide proportion relations between lambdas.
            if (option == USER_LAMBDAS)
                error(['Decide selection between Opts.mu (lambda proportion version) ' ...
                    'or input your own lambdas (user lambdas version).']);
            else
                option = LAMBDA_PROPORTION;
            end
            % mu proportion between penalties.
            mu = Opts.mu;
            % mu must be a row vector.
            if (size(mu, 2) ~= npenalty)
                error('Opts.mu must be a row vector.');
            end
            % There are more than one penalty and many penalties as length(mu).
            if (npenalty <= 1) || (length(mu) ~= npenalty)
                error(['This option is for more than one penalty. Also provide ' ...
                    'many proportions (length(Opts.mu)) as penalties (length(PInfo)).']);
            end
            % mu must be in the range (0, 1) and must add 1.
            if ~isempty(find(mu <= 0)) || (sum(abs(mu)) ~= 1)
                error('mu must be in the range (0, 1) and must add 1.');
            end
        end
        if isfield(Opts, 'nlambdas')
            if (option == USER_LAMBDAS)
                error('The parameter Opts.nlambdas is incompatible with ''user lambdas'' version.');
            end
            nlambdas = Opts.nlambdas;
        end
        if isempty(option)
            % In this case lambdas are selected automatically
            if (length(PInfo) ~= 1)
                error(['This version is only available for the case of single penalty ' ...
                    'applications (Ridge, LASSO, SCAD, ...). No admit combinations.']);
            end
            option = AUTOMATIC_SELECTION;
        end
        if isempty(nlambdas)
            % Default number of lambdas.
            nlambdas = 50;
        end
    end

%% Select lambdas semi-automatically.
    function automatic_lambdas = lambdas_automatically
        switch (option)
            case USER_LAMBDAS
                % Do nothing.
                automatic_lambdas = lambdas;

            case LAMBDA_PROPORTION
                L = [];
                for (it_pinfo = 1:length(PInfo))
                    L = [L; mu(it_pinfo)*PInfo(it_pinfo).L];
                end
                automatic_lambdas = lambdas4lasso(X/L, nlambdas, 1);
                automatic_lambdas = repmat(automatic_lambdas, 1, npenalty).*repmat(mu, nlambdas, 1);

            case AUTOMATIC_SELECTION
                if strcmp(PInfo(1).stype_penalty, 'Ridge')
                    automatic_lambdas = lambdas4ridge(X/PInfo(1).L, nlambdas, 1);
                else
                    automatic_lambdas = lambdas4lasso(X/PInfo(1).L, nlambdas, 1);
                end
        end
    end
end