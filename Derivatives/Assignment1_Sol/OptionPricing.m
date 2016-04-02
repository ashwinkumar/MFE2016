function [ OptionPrice,Delta_Option, B ,Exercise] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div, optionType)
% Returns the option price using the replicating portfolio with the
% following parameters
% S0    - Initial Stock Price
% K     - Strike Price
% u,d   - Upward and downward probability
% N     - Steps
% r     - Risk free rate
% delta - continuous dividend yield
% h     - incremental time period (dt)
% discrete_div - discrete dividend array
% optionType - (American/European Call/Put/Straddle) (A/E C/P/S)
%%
StockPrice_Tree = (get_price_tree(S0,discrete_div,u,d,N));
[PayOff_T,Exercise] = get_PayOffFn(StockPrice_Tree,K, N,optionType);
% Exercise array will contain the N or -1 ( Whether it is in/out of money
% at end of the period
Delta_Option = zeros(N+1,N);
B = zeros(N+1,N);
OptionPrice = zeros(N+1);
OptionPrice(:,N+1) = PayOff_T;
for i = 1: N
    Delta_Option(:,N-i+1) = get_Delta_t(StockPrice_Tree, OptionPrice,delta, h, N-i+1,N);
    B(:,N-i+1) = get_Bond_t(StockPrice_Tree, OptionPrice, Delta_Option(:,N-i+1), r, delta, h, N-i+1,N);
    [OptionPrice(:,N+1-i),Exercise]= get_Option_Price(B(:,N-i+1),Delta_Option(:,N-i+1), StockPrice_Tree(:,N-i+1),K, Exercise, N-i,optionType);
end


end

