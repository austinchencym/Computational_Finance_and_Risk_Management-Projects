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
numPaths = 50000;
numSteps = 12;

% Implement your Black-Scholes pricing formula
[call_BS_European_Price, put_BS_European_Price] = BS_european_price(S0, K, T, r, sigma);

% Implement your one-step Monte Carlo pricing procedure for European option
% numSteps = 1;
[callMC_European_Price_1_step, putMC_European_Price_1_step,Path_1] = MC_european_price(S0, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for European option
[callMC_European_Price_multi_step, putMC_European_Price_multi_step,Path_mul] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths);

% Implement your one-step Monte Carlo pricing procedure for Barrier option
% numSteps = 1;
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step, Path_1_B] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step, Path_mul_B] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);

% Implement your one-step Monte Carlo pricing procedure for Barrier option with increased volatility 10%
% numSteps = 1;
[callMC_Barrier_IncreaseV_Knockin_Price_1_step, putMC_Barrier_IncreaseV_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, 1.1*sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option with increased volatility 10
[callMC_Barrier_IncreaseV_Knockin_Price_multi_step, putMC_Barrier_IncreaseV_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, 1.1*sigma, numSteps, numPaths);

% Implement your one-step Monte Carlo pricing procedure for Barrier option with decreased volatility 10%
% numSteps = 1;
[callMC_Barrier_DecreaseV_Knockin_Price_1_step, putMC_Barrier_DecreaseV_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, 0.9*sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option with decreased volatility 10
[callMC_Barrier_DecreaseV_Knockin_Price_multi_step, putMC_Barrier_DecreaseV_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, 0.9*sigma, numSteps, numPaths);

disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(put_BS_European_Price)])
fprintf('\n');
disp(['One-step MC price of an European call option is ',num2str(callMC_European_Price_1_step)])
disp(['One-step MC price of an European put option is ',num2str(putMC_European_Price_1_step)])
disp(['Multi-step MC price of an European call option is ',num2str(callMC_European_Price_multi_step)])
disp(['Multi-step MC price of an European put option is ',num2str(putMC_European_Price_multi_step)])
fprintf('\n');
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])
fprintf('\n');
disp(['One-step MC price of an Barrier call option with 10% increased Volality is ',num2str(callMC_Barrier_IncreaseV_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option with 10% increased Volality is ',num2str(putMC_Barrier_IncreaseV_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option with 10% increased Volality is ',num2str(callMC_Barrier_IncreaseV_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option with 10% increased Volality is ',num2str(putMC_Barrier_IncreaseV_Knockin_Price_multi_step)])
fprintf('\n');
disp(['One-step MC price of an Barrier call option with 10% decreased Volality is ',num2str(callMC_Barrier_DecreaseV_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option with 10% decreased Volality is ',num2str(putMC_Barrier_DecreaseV_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option with 10% decreased Volality is ',num2str(callMC_Barrier_DecreaseV_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option with 10% decreased Volality is ',num2str(putMC_Barrier_DecreaseV_Knockin_Price_multi_step)])

%%
% Plot results
figure(1);  
numPaths = 50000;
numSteps = 12;
for i=1:numPaths 
    plot(0:numSteps,Path_mul(:,i));
    hold on;
end
title('Multi-Step MC - Underlying Stock Price Simulations');
hold off;


figure(2);  
numPaths = 50000;
numSteps = 1;
for i=1:numPaths 
    plot(0:numSteps,Path_1(:,i));
    hold on;
end
title('1-Step MC - Underlying Stock Price Simulations');
hold off;


figure(3);  
numPaths = 50000;
numSteps = 12;
for i=1:numPaths 
    plot(0:numSteps,Path_mul_B(:,i));
    hold on;
end
hold on;

hline = refline([0 110]);
hline.Color = 'k';
hline.LineWidth = 2.5;
hold off;
title('Multi-Step MC with Barrier - Underlying Stock Price Simulations');


figure(4);  
numPaths = 50000;
numSteps = 1;
for i=1:numPaths 
    plot(0:numSteps,Path_1_B(:,i));
    hold on;
end
hold on;
hline = refline([0 110]);
hline.Color = 'k';
hline.LineWidth = 2.5;
hold off;
title('1-Step MC with Barrier - Underlying Stock Price Simulations');

%% Q1 Black-Scholes pricing model
function [call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma)
    t = 0;
    d1=(log(S0/K)+(r+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
    d2=d1-sigma*sqrt(T-t);
    call_BS_European_Price = normcdf(d1)*S0-normcdf(d2)*K*exp(-r*(T-t));
    putBS_European_Price = normcdf(-d2)*K*exp(-r*(T-t))-normcdf(-d1)*S0;
end
%% Q1 Monte Carlo Geometric Random Walk Path
function paths = GRWPaths(initPrice, mu, sigma, T, numSteps, numPaths)
    % Computes numPaths random paths for a geometric random walk
    % mu    - annual drift
    % sigma - annual volatility
    % T     - total length of time for the path (in years)
    % dT    - time increment (in years) 
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

%% Q1 Monte Carlo pricing model
function [callMC_European_Price_multi_step, putMC_European_Price_multi_step, S] = ...
    MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths)
    
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

%% Q1 Monte Carlo pricing model for knock-in option
function [callMC_Barrier_Knockin_Price_step, putMC_Barrier_Knockin_Price_step, S] = ...
     MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)
    
    % the option becomes a standard European option if the barrier was 
    % crossed some time before expiration.
    % Simulate asset paths for the geometric random walk
    S = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
    determine_barrier=zeros(1,numPaths);
    
    for i=1:numPaths
        if find(S(:,i)>Sb)
            determine_barrier(i)=1;
        else
            determine_barrier(i)=0;
        end
    end
    
    % Calculate the payoff for each path for a Put
    PutPayoff = determine_barrier .*max(K-S(numSteps+1,:),0);
    % Calculate the payoff for each path for a Call
    CallPayoff = determine_barrier .* max(S(numSteps+1,:)-K,0);
    % Discount back
    putMC_Barrier_Knockin_Price_step = mean(PutPayoff)*exp(-r*T);
    callMC_Barrier_Knockin_Price_step = mean(CallPayoff)*exp(-r*T);
end
