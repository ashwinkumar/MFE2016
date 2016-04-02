function [ B ] = get_Bond_t(StockPrice_Tree, Option_Price, Delta, rate, small_delta ,h, t,T)
% Returns the Bond value at given step t, with the
% following parameters
% StockPrice_Tree  -  Stock Price Generated
% OptionPrice   - OptionPrice array (use only for t+1)
% Delta - Delta at t (calulated from get_Delta_t function
% T     - Steps
% r - risk free rate
% h     - incremental time period (dt)
    B = zeros(1,T+1);
    PayOff_t = Option_Price(:,t+1);    
    for i = 1:t
     B(i) = exp(-rate*h)*(PayOff_t(i) - Delta(i)*StockPrice_Tree(i,t+1)*exp(small_delta*h));
    end
    B = transpose(B);
end

