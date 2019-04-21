function  [x_optimal cash_optimal weight cash_reserve] = strat_buy_and_hold(x_init, cash_init, mu, Q, cur_prices, cash_reserve)

   x_optimal = x_init;
   cash_optimal = cash_init;
   portfolio_value = cur_prices*x_init+cash_init;
   
   % newly introduced variables  to help plot and compare variations
   weight = x_optimal*cur_prices/portfolio_value;
   cash_reserve = cash_reserve;

end