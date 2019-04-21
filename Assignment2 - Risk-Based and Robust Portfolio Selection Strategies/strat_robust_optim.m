function  [x_optimal cash_optimal weight] = strat_robust_optim(x_init, cash_init, mu, Q, cur_prices)
    % 20 stocks
    n = 20;
    global year 
    if year == 2008
        r_rf = 0.045;
    else
        r_rf = 0.025;
    end
    % compute the total asset 
    portfolio_value = cur_prices * x_init + cash_init;

    % Define initial portfolio
    w0 = cur_prices' .* x_init / portfolio_value;

    % Bounds on variables
    lb_rMV = zeros(n,1);
    ub_rMV = inf*ones(n,1);

    % Required portfolio robustness
    var_matr = diag(diag(Q));

    % Target portfolio return estimation error is 
    % return estimation error of 1/n portfolio
    % return estimation error of initial portfolio
    rob_init = w0' * var_matr * w0; 
    % target return estimation error
    rob_bnd = rob_init; 

    % Target portfolio return we try risk free rate
    Portf_Retn = r_rf/252;
    
%     % Compute minimum variance portfolio
%     cplex_minVar = Cplex('MinVar');
%     cplex_minVar.addCols(zeros(1,n)', [], lb_rMV, ub_rMV);
%     cplex_minVar.addRows(1, ones(1,n), 1);
%     cplex_minVar.Model.Q = 2*Q;
%     cplex_minVar.Param.qpmethod.Cur = 6;
%     cplex_minVar.DisplayFunc = []; % disable output to screen 
%     cplex_minVar.solve();
%     cplex_minVar.Solution;
%     w_minVar = cplex_minVar.Solution.x; % asset weights
%     ret_minVar = dot(mu, w_minVar);
%     var_minVar = w_minVar' * Q * w_minVar;
%     rob_minVar = w_minVar' * var_matr * w_minVar;

    

    % Formulate and solve robust mean-variance problem
    % Objective function
    f_rMV  = zeros(n,1);
    % Constraints
    A_rMV  = sparse([  mu';
                     ones(1,n)]);
    lhs_rMV = [Portf_Retn; 1];
    rhs_rMV = [inf; 1];
    % Initialize CPLEX environment
    cplex_rMV = Cplex('Robust_MV');
    % Add objective function and variable bounds
    cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
    % Add constraints
    cplex_rMV.addRows(lhs_rMV, A_rMV, rhs_rMV);
    % Add quadratic objective
    cplex_rMV.Model.Q = 2*Q;
    % Add quadratic constraint on return estimation error (robustness constraint)
    Qq_rMV = var_matr;
    cplex_rMV.addQCs(zeros(size(f_rMV)), Qq_rMV, 'L', rob_bnd, {'qc_robust'});
    % Set CPLEX parameters
    cplex_rMV.Param.threads.Cur = 4;
    cplex_rMV.Param.timelimit.Cur = 60;
    cplex_rMV.Param.barrier.qcpconvergetol.Cur = 1e-12; % solution tolerance
    cplex_rMV.DisplayFunc = []; % disable output to screen 
    cplex_rMV.solve();   
    cplex_rMV.Solution;
    
    % if no solution can be obtained during financial crisis --> infeasible
    % keep the portfolio unchanged --> applying "buy and hold strategy"
    if (strcmp(cplex_rMV.Solution.statusstring,'infeasible') == 1)
         x_optimal = x_init;
         cash_optimal = cash_init;
         w_rMV = cur_prices' .* x_init / portfolio_value;
    else        
        w_rMV = cplex_rMV.Solution.x;
    end
   

    % Round near-zero portfolio weights
    w_rMV_nonrnd = w_rMV;
    w_rMV(find(w_rMV<=1e-6)) = 0;
    w_rMV = w_rMV / sum(w_rMV);

    %generate portfolio
    invest_money = w_rMV * portfolio_value;
    x_optimal=floor(invest_money'./cur_prices)';

    %generate transaction cost
    trans = cur_prices*abs(x_optimal-x_init)*0.005;
    cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
    
     if cash_optimal < 0
        % we take the amount of money that still needed         
        % from investment budget, and recalculate x_optimal
        portfolio_value_negative = portfolio_value + cash_optimal*1.2;
        % newly calculated x_optimal since a portion of money been
        % taked to compensate transaction cost
        invest_money = w_rMV * portfolio_value_negative;
        x_optimal=floor(invest_money'./cur_prices)';
        % new transaction cost (would probably lower than needed amount)
        % therefore the solution is feasible but not optimal
        trans = cur_prices*abs(x_optimal-x_init)*0.005;
        % new cash_left
        cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
     end
    
    weight = w_rMV;
 
end