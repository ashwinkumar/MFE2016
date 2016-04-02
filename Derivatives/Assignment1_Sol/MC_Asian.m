function [OptionPrice, StdDeviation] = MC_Asian(S0, K, r, sigma, T, N, n)
% Returns the Asian Option Price and StdDeviation using Monte Carlo Simulations with the
% following parameters
% S0  -  Initial Stock Price
% K   - Strike Price
% T   - Time 
% N   - Steps
% n   - simulations
% sigma - volatality

    randn ('seed',0);
    h = T/(N);
    drift = (r - sigma*sigma/2);
    for i =1:n
        %rng default
        S(:,i) = S0*[1;cumprod(exp(drift*h+sigma*sqrt(h)*randn([N 1])))];
    end
    
    
    discount_factor = exp(-r*T);
    PayOff=discount_factor*max(mean(S(2:N+1,:),1)-K,0);
    OptionPrice =mean(PayOff);
    
    StdDeviation = std(PayOff);
end