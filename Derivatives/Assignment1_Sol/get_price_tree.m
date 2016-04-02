function [ price_tree ] = get_price_tree(S0,discrete_div,u,d,T)
% Returns the stock price tree till T steps
% S0    - Initial Stock Price
% u,d   - Upward and downward probability
% T     - current step
% discrete_div - discrete dividend array
price_tree = zeros(T+1,T+1);
price_tree(1,1) = S0*(1 - discrete_div(1));
%price_tree(2,1) = S0*u*(1- discrete_div(2));
%price_tree(2,2) = S0*d*(1- discrete_div(2));
for i =2:T+1
    price_tree(i,1) = price_tree(i-1,1)*u*(1 - discrete_div(i));
    for j =2:i
        price_tree(i,j) = price_tree(i-1,j-1)*d*(1 - discrete_div(i));
    end
end
price_tree = transpose(price_tree);
end


