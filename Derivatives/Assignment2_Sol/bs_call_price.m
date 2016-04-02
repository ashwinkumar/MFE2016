function [ bs_call_price] = bs_call_price(Stock_price,Strike_price, rate,T,sigma,yield)
%Returns the BS Call Option price using the given paremeters
d1 = getd1(Stock_price, Strike_price, rate,T, sigma);
d2 = getd2(Stock_price, Strike_price, rate,T, sigma);
bs_call_price = Stock_price.*exp(-yield.*T).*normcdf(d1) -Strike_price.*exp(-rate.*T).*normcdf(d2);
end

