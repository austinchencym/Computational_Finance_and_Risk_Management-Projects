function  [x_optimal cash_optimal weight] = strat_equal_risk_contr(x_init, cash_init, mu, Q, cur_prices)
    global Q A_ineq A_eq
    portfolio_value = cur_prices*x_init+cash_init;
    n=20;
    

    % Equality constraints
    A_eq = ones(1,n);
    b_eq = 1;

    % Inequality constraints
    A_ineq = [];
    b_ineql = [];
    b_inequ = [];

    % Define initial portfolio ("equally weighted" or "1/n portfolio")
    w0 = cur_prices' .* x_init / portfolio_value;

    options.lb = zeros(1,n);       % lower bounds on variables
    options.lu = ones (1,n);       % upper bounds on variables
    options.cl = [b_eq' b_ineql']; % lower bounds on constraints
    options.cu = [b_eq' b_inequ']; % upper bounds on constraints
    

    % Set the IPOPT options
    options.ipopt.jac_c_constant        = 'yes';
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.tol                   = 1e-10;
    options.ipopt.print_level           = 0;

    % The callback functions
    funcs.objective         = @computeObjERC;
    funcs.constraints       = @computeConstraints;
    funcs.gradient          = @computeGradERC;
    funcs.jacobian          = @computeJacobian;
    funcs.jacobianstructure = @computeJacobian;

    % !!!! Function "computeGradERC" is just the placeholder
    % !!!! You need to compute the gradient yourself

    %Run IPOPT
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

    w=w_erc;
    
    %generate portfolio
    x_optimal=floor(w*portfolio_value./cur_prices');

    %generate transaction cost
    trans = cur_prices*abs(x_optimal-x_init)*0.005;
    cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
    
     if cash_optimal < 0
        % we take the amount of money that still needed         
        % from investment budget, and recalculate x_optimal
        portfolio_value_negative = portfolio_value + cash_optimal*1.2;
        % newly calculated x_optimal since a portion of money been
        % taked to compensate transaction cost
        x_optimal=floor(w*portfolio_value_negative./cur_prices');
        % new transaction cost (would probably lower than needed amount)
        % therefore the solution is feasible but not optimal
        trans = cur_prices*abs(x_optimal-x_init)*0.005;
        % new cash_left
        cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
     end
    
    weight = w;
    warning('off','all')


end