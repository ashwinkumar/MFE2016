function [ d1 ] = getd1( Stock_price, Strike_price, rate,T, sigma)
%Returns the d1 of B-S option pricing formula
d1 = (log(Stock_price/Strike_price) + (rate + (sigma*sigma)*0.5).*T)./(sigma.*sqrt(T));
end

