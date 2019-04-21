function  [x_optimal cash_optimal weight] = strat_equally_weighted(x_init, cash_init, mu, Q, cur_prices)
    
    %portfolio value
    portfolio_value = cur_prices*x_init+cash_init;
    %divide value into equally 20 subset and get x_optinal corresponding to
    %stock prices
    x_optimal = floor((portfolio_value./cur_prices)/20)';
    
    % transaction cost
    trans = cur_prices*abs(x_optimal-x_init)*0.005;
    
    cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
    
    if cash_optimal < 0
        portfolio_value_negative = portfolio_value + cash_optimal;
        x_optimal = floor((portfolio_value_negative./cur_prices)/20)';
        trans = cur_prices*abs(x_optimal-x_init)*0.005;
        cash_optimal = portfolio_value-cur_prices*x_optimal-trans;
        
    end
   % portfolio_value = cur_prices*x_init+cash_init;
   
    weight = x_optimal*cur_prices/portfolio_value;
end

