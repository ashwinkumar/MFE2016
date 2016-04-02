function [Option_Price,lower_bound, upper_bound] = MC_Option_Price(S_T,K,r,T,optionType,Within_Barrier,n_sims)
% Returns the bound of the option price calculated using the given PayOff
% function and the simulated stock price path
OptionPayOff = zeros(n_sims,1);
for i=1:n_sims
    OptionPayOff(i,1) = Within_Barrier(i)*get_PayOff(S_T(i),K,optionType).*exp(-r(i)*T);
end
Option_Price = mean(OptionPayOff);
lower_bound = Option_Price - 1.96*std(OptionPayOff)/sqrt(n_sims);
upper_bound = Option_Price + 1.96*std(OptionPayOff)/sqrt(n_sims);
end

