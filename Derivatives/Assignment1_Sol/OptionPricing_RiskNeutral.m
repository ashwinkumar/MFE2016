function [ OptionPrice] = OptionPricing_RiskNeutral(S0,K,u,d,N,r,h,discrete_div, optionType)
StockPrice_Tree = (get_price_tree(S0,discrete_div,u,d,N));
PayOff_T = get_PayOffFn(StockPrice_Tree,K, N,optionType);
p_star = (exp(r*h)-d)/(u-d);
q_star = 1-p_star;
OptionPrice = zeros(N+1);
OptionPrice(:,N+1) = PayOff_T;
discount_factor = exp(-r*h);
for i = 1: N
    option_price = zeros(N+1,1);
    for j = 1:N
        option_price(j,1) = discount_factor*(p_star*OptionPrice(j,N-i+2)+ q_star*(OptionPrice(j+1,N-i+2)));
    end
    OptionPrice(:,N+1-i)= max(option_price, get_PayOffFn(StockPrice_Tree,K,N-i,optionType));
end
end

