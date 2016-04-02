function [Option_Price,lower_bound, upper_bound] = MC_Option_Price_Direct(S0,K,r,sigma,T,optionType,n_sims)
% Returns the bound of the option price calculated using the given PayOff
% function and the initial value of stock price
randn('seed',0);
S_T = zeros;
OptionPayOff = zeros;
for i=1:n_sims
    S_T(i) = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*randn);
    OptionPayOff(i) = get_PayOff(S_T(i),K,optionType)*exp(-r*T);
end
Option_Price = mean(OptionPayOff);
lower_bound = Option_Price - 1.96*std(OptionPayOff)/sqrt(n_sims);
upper_bound = Option_Price + 1.96*std(OptionPayOff)/sqrt(n_sims);
end

