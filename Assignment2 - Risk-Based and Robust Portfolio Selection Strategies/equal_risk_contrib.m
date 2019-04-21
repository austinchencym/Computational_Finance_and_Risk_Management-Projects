clear all;
clc;
format long;

global Q A_ineq A_eq

%% Add PATH to CPLEX and other solvers
addpath('C:\Austin\Applications\CPLEX_Optimization_Studio\cplex\matlab\x64_win64');

% Random data for 10 stocks
n = 10;
if(exist('portf_data.mat','file')) % read data from file if it exists
    load('portf_data', 'Q', 'mu')
else
    Q = randn(n); Q = Q*Q'/1000; % covariance matrix
    mu  = rand(1,n)'/100;        % expected return
end

% Equality constraints
A_eq = ones(1,n);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];
           
% Define initial portfolio ("equally weighted" or "1/n portfolio")
w0 = repmat(1.0/n, n, 1);

options.lb = zeros(1,n);       % lower bounds on variables
options.lu = ones (1,n);       % upper bounds on variables
options.cl = [b_eq' b_ineql']; % lower bounds on constraints
options.cu = [b_eq' b_inequ']; % upper bounds on constraints

% Set the IPOPT options
options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-10;

% The callback functions
funcs.objective         = @computeObjERC;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradERC;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @computeJacobian;

% !!!! Function "computeGradERC" is just the placeholder
% !!!! You need to compute the gradient yourself
  
%% Run IPOPT
[wsol info] = ipopt(w0',funcs,options);

% Make solution a column vector
if(size(wsol,1)==1)
    w_erc = wsol';
else
    w_erc = wsol;
end

% Compute return, variance and risk contribution for the ERC portfolio
ret_ERC = dot(mu, w_erc);
var_ERC = w_erc'*Q*w_erc;
 RC_ERC = (w_erc .* ( Q*w_erc )) / sqrt(w_erc'*Q*w_erc);

%% 1/n portfolio return
ret_init = dot(mu, w0);
%% 1/n portfolio variance
var_init = w0' * Q * w0;

% Bounds on variables
lb_rMV = zeros(n,1);
ub_rMV = inf*ones(n,1);

%% Compute minimum variance portfolio
cplex_minVar = Cplex('MinVar');
cplex_minVar.addCols(zeros(1,n)', [], lb_rMV, ub_rMV);
cplex_minVar.addRows(1, ones(1,n), 1);
cplex_minVar.Model.Q = 2*Q;
cplex_minVar.Param.qpmethod.Cur = 6;
cplex_minVar.solve();
cplex_minVar.Solution
w_minVar = cplex_minVar.Solution.x; % asset weights
ret_minVar = dot(mu, w_minVar);
var_minVar = w_minVar' * Q * w_minVar;
RC_minVar  = (w_minVar .* ( Q*w_minVar )) / sqrt(w_minVar' * Q * w_minVar);
RC_init  = (w0 .* ( Q*w0 )) / sqrt(w0' * Q * w0);

disp(' ')
disp(['       Portfolio ERC return = ' num2str(ret_ERC,'%9.5f')])
disp(['    Portfolio minVar return = ' num2str(ret_minVar,'%9.5f')])
disp(['       Portfolio 1/n return = ' num2str(ret_init,'%9.5f')])   
disp(['      Portfolio ERC st.dev. = ' num2str(sqrt(var_ERC),'%9.5f')])
disp(['   Portfolio minVar st.dev. = ' num2str(sqrt(var_minVar),'%9.5f')])
disp(['      Portfolio 1/n st.dev. = ' num2str(sqrt(var_init),'%9.5f')])
disp(' ')

fprintf('\n\nPortfolio weights for ERC, minVar and 1/n portfolios:\n')
[w_erc w_minVar w0]

fprintf('\n\nAsset risk contributions for ERC, minVar and 1/n portfolios:\n')
[RC_ERC RC_minVar RC_init]

fprintf('\n\nSum of asset risk contributions for ERC, minVar and 1/n portfolios:\n')
[sum(RC_ERC) sum(RC_minVar) sum(RC_init)]

fprintf('\n\nStandard deviation for ERC, minVar and 1/n portfolios:\n')
[sqrt(var_ERC) sqrt(var_minVar) sqrt(var_init)]
