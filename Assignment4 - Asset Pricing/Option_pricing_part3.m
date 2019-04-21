clc;
clear all;
format long

% Pricing a European option using Black-Scholes formula and Monte Carlo simulations
% Pricing a Barrier option using Monte Carlo simulations

S0 = 100;     % spot price of the underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.2;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
Sb = 110;     % barrier

% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations

% Implement your Black-Scholes pricing formula
[call_BS_European_Price, put_BS_European_Price] = BS_european_price(S0, K, T, r, sigma);
disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(put_BS_European_Price)])
fprintf('\n');

%% Part 3
% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations
numPaths_list = [10000 20000 30000 40000 50000 60000 70000 80000 90000 100000];
numSteps_list = [2 12 24 52 252];
%% one-step MC
% we need to compare  call_BS_European_Price, put_BS_European_Price
% minimize difference
residual_call = zeros(1,length(numPaths_list));
residual_put = zeros(1,length(numPaths_list));
for i = 1:length(numPaths_list)
    % Implement your one-step Monte Carlo pricing procedure for European option
    [callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, 1, numPaths_list(i));
    residual_call(i) = abs(call_BS_European_Price-callMC_European_Price_1_step);
    residual_put(i) = abs(put_BS_European_Price-putMC_European_Price_1_step);
end
[minval_call, minidx_call] = min(residual_call);
[minval_put, minidx_put] = min(residual_put);
optimal_num_path_call = numPaths_list(minidx_call);
optimal_num_path_put = numPaths_list(minidx_put);
disp(['1-step MC Call Option, Optimal numPaths is ',num2str(optimal_num_path_call)])
disp(['                       With a Minimum difference of ',num2str(abs(call_BS_European_Price-callMC_European_Price_1_step))])
disp(['1-step MC Put Option, Optimal numPaths is ',num2str(optimal_num_path_put)])
disp(['                      With a Minimum difference of ',num2str(abs(put_BS_European_Price-putMC_European_Price_1_step))])
fprintf('\n');
%% multi-step
residual_mul_call = zeros(length(numSteps_list),length(numPaths_list));
residual_mul_put = zeros(length(numSteps_list),length(numPaths_list));
for i = 1:length(numSteps_list)
    for j = 1:length(numPaths_list)
        [callMC_European_Price_multi_step, putMC_European_Price_multi_step] = MC_european_price(S0, K, T, r, mu, sigma, numSteps_list(i), numPaths_list(j));
        residual_mul_call(i,j) = abs(call_BS_European_Price-callMC_European_Price_multi_step);
        residual_mul_put(i,j) = abs(put_BS_European_Price-putMC_European_Price_multi_step);
    end
end
[minval_mul_row_call, minidx_mul_row_call] = min(residual_mul_call);
[minval_mul_final_call, minidx_mul_final_call] = min(min(residual_mul_call));
[minval_mul_row_put, minidx_mul_row_put] = min(residual_mul_put);
[minval_mul_final_put, minidx_mul_final_put] = min(min(residual_mul_put));

optimal_num_mul_step_call = numSteps_list(minidx_mul_row_call(1));
optimal_num_mul_path_call = numPaths_list(minidx_mul_final_call);
optimal_num_mul_step_put = numSteps_list(minidx_mul_row_put(1));
optimal_num_mul_path_put = numPaths_list(minidx_mul_final_put);

disp(['Multi-step MC Call Option, Optimal numStep is ',num2str(optimal_num_mul_step_call)])
disp(['                           Optimal numPaths is ',num2str(optimal_num_mul_path_call)])
disp(['                           With a Minimum difference of ',num2str(abs(call_BS_European_Price-callMC_European_Price_multi_step))])
disp(['Multi-step MC Put Option, Optimal numStep is ',num2str(optimal_num_mul_step_put)])
disp(['                          Optimal numPaths is ',num2str(optimal_num_mul_path_put)])
disp(['                          With a Minimum difference of ',num2str(abs(put_BS_European_Price-putMC_European_Price_multi_step))])
fprintf('\n');

%% functions Black-Scholes pricing model
function [call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma)
    t = 0;
    d1=(log(S0/K)+(r+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
    d2=d1-sigma*sqrt(T-t);
    call_BS_European_Price = normcdf(d1)*S0-normcdf(d2)*K*exp(-r*(T-t));
    putBS_European_Price = normcdf(-d2)*K*exp(-r*(T-t))-normcdf(-d1)*S0;
end

%% Monte Carlo Geometric Random Walk Path
function paths = GRWPaths(initPrice, mu, sigma, T, numSteps, numPaths)
    % Computes numPaths random paths for a geometric random walk
    % mu is the annual drift, sigma the annual volatility
    % T is the total length of time for the path (in years)
    % dT is the time increment (in years)
       
    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    
    % Vector of paths will store realizations of the asset price
    % First asset price is the initial price
    paths(1,:) = initPrice;
 
    % Generate paths
    for iPath = 1:numPaths
        for iStep = 1:numSteps
            paths(iStep+1, iPath) = paths(iStep, iPath) * exp( (mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1) );
        end
    end 
end

%% Monte Carlo pricing model
function [callMC_European_Price_multi_step, putMC_European_Price_multi_step, S] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths)
    
    % Simulate asset paths for the geometric random walk
    S = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
    % Calculate the payoff for each path for a Put
    PutPayoff = max(K-S(numSteps+1,:),0);
    % Calculate the payoff for each path for a Call
    CallPayoff = max(S(numSteps+1,:)-K,0);
    % Discount back
    putMC_European_Price_multi_step = mean(PutPayoff)*exp(-r*T);
    callMC_European_Price_multi_step = mean(CallPayoff)*exp(-r*T);
end

