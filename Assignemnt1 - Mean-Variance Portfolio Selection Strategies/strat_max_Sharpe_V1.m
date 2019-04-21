function  [x_optimal cash_optimal weight cash_reserve] = strat_max_Sharpe_V1(x_init, cash_init, mu, Q, cur_prices, cash_reserve)

    %Q2 = [Q zeros(size(Q,1),1)];
    portfolio_value = cur_prices*x_init+cash_init;
%     invest_asset=portfolio_value*0.99;
%     trans_asset=portfolio_value*0.01;
    
    % Reshape Q tobe 21*21 to allow computing K
    Q=[Q,zeros(20,1)]; % add a row of zero since K does not appear in objective function
    Q=[Q;zeros(21,1)']; % add a column of zero as well
    
    % Optimization problem data
    lb=zeros(21,1);
    ub=inf*ones(21,1);
    A1=[mu-0.025/252, ones(20,1)];% first constraint and second constraint coefficient for y
    A = [A1;[0,-1]]'; % add coefficient for k so that in the second constraint: lhs - k = 0
    b=[1;0];
    
    % Compute minimum variance portfolio
    cplex2 = Cplex('max_Sharpe');
    cplex2.addCols(zeros(21,1), [], lb, ub); % Add objective function and bounds on variables 
    cplex2.addRows(b, A, b);% Add constraints to CPLEX model
    cplex2.Model.Q = 2*Q; % Add quadratic part of objective function to CPLEX model
    cplex2.Param.qpmethod.Cur = 6; % concurrent algorithm
    cplex2.Param.barrier.crossover.Cur = 1; % enable crossover
    cplex2.DisplayFunc = []; % disable output to screen
    cplex2.solve();
    
    y=cplex2.Solution.x(1:20);
    k=cplex2.Solution.x(21);
    
    %generate weights w=y/k
    w=y/k;
    
    %generate portfolio
    x_optimal=floor(w*portfolio_value./cur_prices');

    %generate transaction cost
    trans = cur_prices*abs(x_optimal-x_init)*0.005;
    cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
    
    
    % if we have a negative cash_optimal
    % it means our current cash are not sufficient to afford transaction cost 
    if cash_optimal < 0
        % if the cash_reserve can pay the negative cash_optimal
        % then we don't need to take money from investment budget and
        % recalculate optimal postitions
        if cash_reserve >= abs(cash_optimal)
            % subtract cash reserve with corresponding amount of needed
            % money, and set the cash optimal as 0, it equal to take needed
            % amount of money from cash reserve and make sure the
            % transaction can be processed.
            cash_reserve = cash_reserve + cash_optimal;
            cash_optimal = 0;
        else
            % if the reserved cash cannot pay for the needed negative cash
            % optimal, we take the rest of all reserved cash and add it to
            % portfolio, then we take the amount of money that still needed
            % from investment budget, and recalculate x_optimal
            portfolio_value_negative = portfolio_value + cash_reserve +cash_optimal;
            % newly calculated x_optimal since a portion of money been
            % taked to compensate transaction cost
            x_optimal=floor(w*portfolio_value_negative./cur_prices');
            % new transaction cost (would probably lower than needed amount)
            % therefore the solution is feasible but not optimal
            trans = cur_prices*abs(x_optimal-x_init)*0.005;
            % new cash_left
            cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
            % used up all cash reserve
            cash_reserve = 0;
        end 
    else
        cash_reserve=cash_reserve; % if current cash is enough, keep cash reserve
    end
    
    
    weight = w;
    


end