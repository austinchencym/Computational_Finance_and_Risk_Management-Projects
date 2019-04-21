function  [x_optimal cash_optimal weight cash_reserve] = start_min_variance(x_init, cash_init, mu, Q, cur_prices, cash_reserve)
    portfolio_value = cur_prices*x_init+cash_init;
    
    % Optimization problem data
    lb=zeros(20,1);
    ub=inf*ones(20,1);
    A=ones(20,1)';
    b=1;
    
    % Compute minimum variance portfolio
    cplex1 = Cplex('min_Variance');
    cplex1.addCols(zeros(20,1), [], lb, ub);
    cplex1.addRows(b, A, b);
    cplex1.Model.Q = 2*Q;
    cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
    cplex1.Param.barrier.crossover.Cur = 1; % enable crossover
    cplex1.DisplayFunc = []; % disable output to screen
    cplex1.solve();
    w=cplex1.Solution.x;
    
    %generate portfolio
    x_optimal=floor(w*portfolio_value./cur_prices');

    %generate transaction cost
    trans = cur_prices*abs(x_optimal-x_init)*0.005;
    cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
    
     if cash_optimal < 0
        % we take the amount of money that still needed         
        % from investment budget, and recalculate x_optimal
        portfolio_value_negative = portfolio_value + cash_optimal;
        % newly calculated x_optimal since a portion of money been
        % taked to compensate transaction cost
        x_optimal=floor(w*portfolio_value_negative./cur_prices');
        % new transaction cost (would probably lower than needed amount)
        % therefore the solution is feasible but not optimal
        trans = cur_prices*abs(x_optimal-x_init)*0.005;
        % new cash_left
        cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
    
    weight = w;
    cash_reserve=cash_reserve;

end