clc;
clear all;
format long
global period Q
% Input files
input_file_prices  = 'Daily_closing_prices.csv';

% Add path to CPLEX
addpath('C:/Austin/Applications/CPLEX_Optimization_Studio/cplex/matlab/x64_win64');

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);

% Find the number of trading days in Nov-Dec 2014 and
% compute expected return and covariance matrix for period 1
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
cur_returns0 = data_prices(day_ind_start0+1:day_ind_end0,:) ./ data_prices(day_ind_start0:day_ind_end0-1,:) - 1;
mu = mean(cur_returns0)';
Q = cov(cur_returns0);

% Remove datapoints for year 2014
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Initial positions in the portfolio
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';

% Initial value of the portfolio
init_value = data_prices(1,:) * init_positions;
fprintf('\nInitial portfolio value = $ %10.2f\n\n', init_value);

% Initial portfolio weights
w_init = (data_prices(1,:) .* init_positions')' / init_value;

% Number of periods, assets, trading days
N_periods = 6*length(unique(dates_array(:,1))); % 6 periods per year
N = length(tickers);
N_days = length(dates);

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;

% Number of strategies
strategy_functions = {'strat_buy_and_hold' 'strat_equally_weighted' 'strat_min_variance' 'strat_max_Sharpe' 'strat_equal_risk_contr' 'strat_lever_equal_risk_contr' 'strat_robust_optim'};
strategy_names     = {'Buy and Hold' 'Equally Weighted Portfolio' 'Mininum Variance Portfolio' 'Maximum Sharpe Ratio Portfolio' 'Equal risk contributions' 'Leverage equal risk contributions' 'Robust Optimization'};
N_strat = 7; % comment this in your code
%N_strat = length(strategy_functions); % uncomment this in your code
fh_array = cellfun(@str2func, strategy_functions, 'UniformOutput', false);


for (period = 1:N_periods)
   % Compute current year and month, first and last day of the period
   if(dates_array(1,1)==15)
       cur_year  = 15 + floor(period/7);
   else
       cur_year  = 2015 + floor(period/7);
   end
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   fprintf('\nPeriod %d: start date %s, end date %s\n', period, char(dates(day_ind_start)), char(dates(day_ind_end)));

   % Prices for the current day
   cur_prices = data_prices(day_ind_start,:);

   % Execute portfolio selection strategies
   for(strategy = 1:N_strat)

       
      
      % Get current portfolio positions
      if(period==1) %initiation
         curr_positions = init_positions;
         curr_cash = 0;
         portf_value{strategy} = zeros(N_days,1);
      else
         curr_positions = x{strategy,period-1};
         curr_cash = cash{strategy,period-1};
      end

      % Compute strategy
      [x{strategy,period} cash{strategy,period} weight{strategy,period}] = fh_array{strategy}(curr_positions, curr_cash, mu, Q, cur_prices);

      % Verify that strategy is feasible (you have enough budget to re-balance portfolio)
      % Check that cash account is >= 0
      % Check that we can buy new portfolio subject to transaction costs

      %%%%%%%%%%% Insert your code here %%%%%%%%%%%%

       % if the cash for specfic strategy in specific period is negative, 
       % it probably means the transaction cost for that period is great
       % enough that we could not ignore, therefore,  witin in each function, 
       % few new lines are introduced to help compensate the cash by properly 
       % reducing the amount of money to invest(rebalance).
       
         
      % Compute portfolio value in Q2
      if strategy == 6
        portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period} - 1000000; 
      else
        portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period};
      end 
      %portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period};

      fprintf('   Strategy "%s", value begin = $ %10.2f, value end = $ %10.2f\n', char(strategy_names{strategy}), portf_value{strategy}(day_ind_start), portf_value{strategy}(day_ind_end));

   end
      
   % Compute expected returns and covariances for the next period
   cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
   mu = mean(cur_returns)';
   Q = cov(cur_returns);
   
   
   

   
   
end



% Plot results
% figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%

% assign result values
strat_result_1 = portf_value{1};
strat_result_2 = portf_value{2};
strat_result_3 = portf_value{3};
strat_result_4 = portf_value{4};
strat_result_5 = portf_value{5};
strat_result_6 = portf_value{6};
strat_result_7 = portf_value{7};


% figure 1
% plot results for all 4 strategies 
figure(1);
plot(strat_result_1,"b")
hold on
plot(strat_result_2,"k")
hold on
plot(strat_result_3,"r")
hold on
plot(strat_result_4,"g")
hold on
plot(strat_result_5,"y")
hold on
plot(strat_result_6,"c")
hold on
plot(strat_result_7,"m")
xlim([0 504])
legend('Buy and Hold','Equally Weighted', 'Min Variance', 'Max Sharpe','Equal Risk', 'Leveraged Equal Risk', 'Robust Optimization');
xlabel('Trading Days');
ylabel('Portfolio Value ($)'); 
title('Portfolio Value VS Trading day');
hold off;

% figure 2 dynamic change for strategy 7
x_7 = [weight{7,:}];
x_7 = [w_init x_7];
figure(2);
plot([0:12],x_7');
ylim([0 1]);
legend('Stock 1','Stock 2', 'Stock 3', 'Stock 4','Stock 5','Stock 6','Stock 7','Stock 8','Stock 9','Stock 10','Stock 11','Stock 12','Stock 13','Stock 14','Stock 15','Stock 16','Stock 17','Stock 18','Stock 19','Stock 20');
xlabel('Period'); 
ylabel('Stock Weight'); 
title('Stock Weight VS Period for strategy 7');
% 
% figure 3 dynamic change for strategy 3
figure(3);
x_3 = [weight{3,:}];
x_3 = [w_init x_3];
plot([0:12],x_3');
ylim([0 1]);
legend('Stock 1','Stock 2', 'Stock 3', 'Stock 4','Stock 5','Stock 6','Stock 7','Stock 8','Stock 9','Stock 10','Stock 11','Stock 12','Stock 13','Stock 14','Stock 15','Stock 16','Stock 17','Stock 18','Stock 19','Stock 20');
xlabel('Period'); 
ylabel('Stock Weight'); 
title('Stock Weight VS Period for strategy 3');

% figure 4 dynamic change for strategy 4
figure(4);
x_4 = [weight{4,:}];
x_4 = [w_init x_4];
plot([0:12],x_4');
ylim([0 1]);
legend('Stock 1','Stock 2', 'Stock 3', 'Stock 4','Stock 5','Stock 6','Stock 7','Stock 8','Stock 9','Stock 10','Stock 11','Stock 12','Stock 13','Stock 14','Stock 15','Stock 16','Stock 17','Stock 18','Stock 19','Stock 20');
xlabel('Period'); 
ylabel('Stock Weight'); 
title('Stock Weight VS Period for strategy 4');

figure(5);
plot(strat_result_5,"b")
hold on
plot(strat_result_6,"c")
hold on
plot(strat_result_7,"m")
xlim([0 504])
legend('Equal Risk', 'Leveraged Equal Risk', 'Robust Optimization');
xlabel('Trading Days');
ylabel('Portfolio Value ($)'); 
title('Portfolio Value for strategy 5,6,7 VS Trading days');
hold off;
% 
% % figure 4 Q3
% figure(4);
% plot(strat_result_4,"g")
% hold on
% plot(strat_result_5,"b")
% xlim([0 504])
% legend('Max Sharpe','Max Sharpe V1');
% xlabel('Trading Days');
% ylabel('Portfolio Value ($)'); 
% title('Portfolio Value VS Trading day');
% hold off;
%  
