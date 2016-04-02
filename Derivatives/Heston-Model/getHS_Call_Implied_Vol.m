function [ implied_vol ] = getHS_Call_Implied_Vol( S,K,r,T,yield,CallPrice,guess)
% Return the implied vol of BS Option Price using Heston Model's price

result = @(x) bs_call_price(S,K,r,T,x,yield) - CallPrice;
implied_vol = fsolve(result,guess);
end

