function [ d2 ] = getd2( Stock_price, Strike_price, rate,T, sigma)
%Returns the d2 of B-S option pricing formula
d2 = getd1(Stock_price, Strike_price, rate,T, sigma)- sigma.*sqrt(T);
end