function [ Delta ] = get_Delta_t(StockPrice_Tree, Option_Price,small_delta, h, t,T)
% Returns the Delta of Portfolio at given step t, with the
% following parameters
% StockPrice_Tree  -  Stock Price Generated
% OptionPrice   - OptionPrice array (use only for t+1)
% T     - Steps
% small_delta - continuous dividend yield
% h     - incremental time period (dt)
    Delta = zeros(1,T+1);
    PayOff_t = Option_Price(1:t+1,t+1);    
    for i = 1:t
     Delta(i) = exp(-small_delta*h)*(PayOff_t(i) - PayOff_t(i+1))/ (StockPrice_Tree(i,t+1) - StockPrice_Tree(i+1,t+1));
    end
    Delta = transpose(Delta);
end


